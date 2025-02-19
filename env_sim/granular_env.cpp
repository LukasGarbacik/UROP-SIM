//  Copyright (c) 2021, SBEL GPU Development Team
//  Copyright (c) 2021, University of Wisconsin - Madison
//
//  SPDX-License-Identifier: BSD-3-Clause

#include <core/ApiVersion.h>
#include <core/utils/ThreadManager.h>
#include <DEM/API.h>
#include <DEM/HostSideHelpers.hpp>
#include <DEM/utils/Samplers.hpp>

#include <cstdio>
#include <chrono>
#include <filesystem>
#include <random>

using namespace deme;
using namespace std::filesystem;

int main() {
    // Simulation parameters, could be changed
    double terrain_rad = 0.0025 / 2.0; // radius of granular material
    int run_num = 0;

    // Random number generator initialized once
    std::random_device rd;
    std::mt19937 gen(rd());

    // Solver setup
    DEMSolver DEMSim;
    DEMSim.SetVerbosity("ERROR");
    DEMSim.SetOutputFormat("CSV");
    DEMSim.SetOutputContent({"ABSV"});
    DEMSim.SetMeshOutputFormat("VTK");

    // Output directory
    path out_dir = current_path() / "GranularEnvOutput";
    create_directory(out_dir);

    // Define material properties
    /*
    E (Young's Modulus): 70 MPa, measures stiffness.
    nu (Poisson's Ratio): 0.24, describes the material's deformation behavior.
    CoR (Coefficient of Restitution): 0.9, indicates how elastic collisions are.
    mu (Friction Coefficient): 0.3.
    Crr (Rolling Resistance): 0.
    */
    auto mat_type_rover = DEMSim.LoadMaterial({{"E", 7e7}, {"nu", 0.24}, {"CoR", 0.9}, {"mu", 0.3}, {"Crr", 0.0}});
    auto mat_type_terrain = DEMSim.LoadMaterial({{"E", 7e7}, {"nu", 0.24}, {"CoR", 0.9}, {"mu", 0.3}, {"Crr", 0.0}});
    auto mat_type_terrain_sim = DEMSim.LoadMaterial({{"E", 7e7}, {"nu", 0.24}, {"CoR", 0.9}, {"mu", 0.3}, {"Crr", 0.0}});

    float step_size = 2e-6;
    double world_size = 0.2;
    DEMSim.InstructBoxDomainDimension({-world_size / 2., world_size / 2.}, {-world_size / 2., world_size / 2.}, {0, 10 * world_size});
    DEMSim.InstructBoxDomainBoundingBC("top_open", mat_type_terrain);

    // Load SolidWorks-designed rover model
    auto rover = DEMSim.AddWavefrontMeshObject((GET_DATA_PATH() / "mesh/rover.obj").string(), mat_type_rover);
    rover->Scale(1.0);  // Adjust scale as per model dimensions
    rover->SetInitPos(make_float3(0, 0, 0.1));  // Start just above the terrain
    rover->SetMass(1.0);  // Total mass in kg
    rover->SetMOI(make_float3(0.1, 0.1, 0.1));  // Adjusted MOI for typical rover
    auto rover_tracker = DEMSim.Track(rover);

    std::vector<std::shared_ptr<DEMClumpTemplate>> templates_terrain;
    for (int i = 0; i < 11; i++) {
        templates_terrain.push_back(DEMSim.LoadSphereType(
            terrain_rad * terrain_rad * terrain_rad * 2.5e3 * 4 / 3 * PI, terrain_rad, mat_type_terrain));
        terrain_rad += 0.0001 / 2.;
    }

    unsigned int num_particle = 0;
    float sample_z = 1.5 * terrain_rad;
    float fullheight = world_size * 2.;
    float sample_halfwidth = world_size / 2 - 2 * terrain_rad;

    if (run_num > 0) {
        char cp_filename[200];
        sprintf(cp_filename, "%s/bed.csv", out_dir.c_str());

        auto clump_xyz = DEMSim.ReadClumpXyzFromCsv(std::string(cp_filename));
        auto clump_quaternion = DEMSim.ReadClumpQuatFromCsv(std::string(cp_filename));
        for (int i = 0; i < templates_terrain.size(); i++) {
            char t_name[20];
            sprintf(t_name, "%04d", i);

            auto this_xyz = clump_xyz[std::string(t_name)];
            auto this_quaternion = clump_quaternion[std::string(t_name)];
            auto batch = DEMSim.AddClumps(templates_terrain[i], this_xyz);
            batch->SetOriQ(this_quaternion);
            num_particle += this_quaternion.size();
        }
    } else {
        std::uniform_int_distribution<> dist(0, templates_terrain.size() - 1);

        PDSampler sampler(2.01 * terrain_rad);
        while (sample_z < fullheight) {
            float3 sample_center = make_float3(0, 0, sample_z);
            auto input_xyz = sampler.SampleBox(sample_center, make_float3(sample_halfwidth, sample_halfwidth, 0.000001));
            std::vector<std::shared_ptr<DEMClumpTemplate>> template_to_use(input_xyz.size());
            for (unsigned int i = 0; i < input_xyz.size(); i++) {
                template_to_use[i] = templates_terrain[dist(gen)];
            }
            DEMSim.AddClumps(template_to_use, input_xyz);
            num_particle += input_xyz.size();
            sample_z += 2.01 * terrain_rad;
        }
    }

    std::cout << "Total num of particles: " << num_particle << std::endl;

    // Now add a plane to compress the sample
    auto compressor = DEMSim.AddExternalObject();
    compressor->AddPlane(make_float3(0, 0, 0), make_float3(0, 0, -1), mat_type_terrain);
    compressor->SetFamily(10);
    DEMSim.SetFamilyFixed(10);
    DEMSim.DisableContactBetweenFamilies(0, 10);
    auto compressor_tracker = DEMSim.Track(compressor);

    auto max_z_finder = DEMSim.CreateInspector("clump_max_z");
    auto total_mass_finder = DEMSim.CreateInspector("clump_mass");

    DEMSim.SetInitTimeStep(step_size);
    DEMSim.SetMaxVelocity(30.);
    DEMSim.SetGravitationalAcceleration(make_float3(0, 0, -9.81));

    DEMSim.Initialize();

    float sim_time = 3.0;
    float settle_time = 1.0;
    unsigned int fps = 10;
    float frame_time = 1.0 / fps;
    unsigned int currframe = 0;

    if (run_num == 0) {
        // We can let it settle first
        for (float t = 0; t < settle_time; t += frame_time) {
            std::cout << "Frame: " << currframe << std::endl;
            char filename[200], meshfilename[200];
            sprintf(filename, "%s/DEMdemo_output_%04d.csv", out_dir.c_str(), currframe);
            sprintf(meshfilename, "%s/DEMdemo_mesh_%04d.vtk", out_dir.c_str(), currframe);
            DEMSim.WriteSphereFile(std::string(filename));
            DEMSim.WriteMeshFile(std::string(meshfilename));
            currframe++;

            DEMSim.DoDynamicsThenSync(frame_time);
            DEMSim.ShowThreadCollaborationStats();
        }

        char cp_filename[200];
        sprintf(cp_filename, "%s/bed.csv", out_dir.c_str());
        DEMSim.WriteClumpFile(std::string(cp_filename));
    }

    // Simulate rover rotation
    std::cout << "Starting rover rotation simulation..." << std::endl;
    float rotation_sim_time = 5.0;
    float3 angular_velocity = make_float3(0.0, 0.0, 5.0); // Rotating around Z-axis
    rover->SetAngVel(angular_velocity);

    for (float t = 0; t < rotation_sim_time; t += frame_time) {
        char filename[200], meshfilename[200];
        sprintf(filename, "%s/RoverRotation_output_%04d.csv", out_dir.c_str(), currframe);
        sprintf(meshfilename, "%s/RoverRotation_mesh_%04d.vtk", out_dir.c_str(), currframe);
        DEMSim.WriteSphereFile(std::string(filename));
        DEMSim.WriteMeshFile(std::string(meshfilename));
        currframe++;

        DEMSim.DoDynamics(frame_time);
        DEMSim.ShowThreadCollaborationStats();
    }

    float3 final_pos = rover_tracker->Pos();
    float3 final_ori = rover_tracker->AngVel();
    std::cout << "Final rover position: (" << final_pos.x << ", " << final_pos.y << ", " << final_pos.z << ")" << std::endl;
    std::cout << "Final rover angular velocity: (" << final_ori.x << ", " << final_ori.y << ", " << final_ori.z << ")" << std::endl;

    std::cout << "Granular environment and rover rotation simulation completed." << std::endl;
    return 0;
}

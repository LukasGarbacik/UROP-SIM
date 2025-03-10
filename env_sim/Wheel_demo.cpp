#include <DEM/API.h>
#include <DEM/HostSideHelpers.hpp>
#include <DEM/utils/Samplers.hpp>

#include <chrono>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <iostream>
#include <memory>
#include <random>

using namespace deme;

const double math_PI = 3.1415927;

// === USER CONFIGURABLE ===
const std::string custom_wheel_obj = "path/to/your_wheel.obj";  // Replace with your wheel OBJ file path!
const std::string out_folder = "DemoOutput_CustomWheel";

// === MAIN FUNCTION ===
int main() {
    // Output directory setup
    std::filesystem::path out_dir = std::filesystem::current_path() / out_folder;
    std::filesystem::create_directory(out_dir);

    // Create DEM solver
    DEMSolver DEMSim;
    DEMSim.SetVerbosity(INFO);
    DEMSim.SetOutputFormat(OUTPUT_FORMAT::CSV);
    DEMSim.SetOutputContent(OUTPUT_CONTENT::ABSV);
    DEMSim.SetMeshOutputFormat(MESH_FORMAT::VTK);
    DEMSim.SetContactOutputContent({"OWNER", "FORCE", "POINT"});

    // Optional: If you don't need individual contact forces for each contact pair
    DEMSim.SetNoForceRecord();

    // Define materials for wheel
    auto mat_wheel = DEMSim.LoadMaterial({{"E", 1e9}, {"nu", 0.3}, {"CoR", 0.5}, {"mu", 0.5}, {"Crr", 0.01}});
    auto mat_terrain = DEMSim.LoadMaterial({{"E", 1e9}, {"nu", 0.3}, {"CoR", 0.4}, {"mu", 0.5}, {"Crr", 0.01}});
    DEMSim.SetMaterialPropertyPair("mu", mat_wheel, mat_terrain, 0.8);
    DEMSim.SetMaterialPropertyPair("CoR", mat_wheel, mat_terrain, 0.6);

    // World setup
    float G = 9.81;
    float bottom = -0.5;
    DEMSim.InstructBoxDomainDimension(2.0, 1.0, 2.0);
    DEMSim.InstructBoxDomainBoundingBC("top_open", mat_terrain);
    DEMSim.AddBCPlane(make_float3(0, 0, bottom), make_float3(0, 0, 1), mat_terrain);

    // Load and set up the wheel from your OBJ file
    float wheel_radius = 0.25;
    float wheel_width = 0.2;
    float wheel_weight = 100.0;  // Newtons
    float wheel_mass = wheel_weight / G;
    float wheel_IYY = wheel_mass * wheel_radius * wheel_radius / 2;
    float wheel_IXX = (wheel_mass / 12.0f) * (3 * wheel_radius * wheel_radius + wheel_width * wheel_width);

    auto wheel = DEMSim.AddWavefrontMeshObject(custom_wheel_obj, mat_wheel);
    wheel->SetMass(wheel_mass);
    wheel->SetMOI(make_float3(wheel_IXX, wheel_IYY, wheel_IXX));
    wheel->SetFamily(1);  // Used for prescribing motion
    auto wheel_tracker = DEMSim.Track(wheel);

    // Terrain particle setup (spheres only)
    float terrain_density = 2600;             // kg/m^3
    float sphere_radius = 0.01;               // 1 cm radius spheres
    auto terrain_sphere_template = DEMSim.LoadSphereType(terrain_density, sphere_radius, mat_terrain);

    // Use HCP Sampler for terrain particle positions
    HCPSampler sampler(sphere_radius * 2.1);  // Packing factor (~2 * radius)
    float sample_halfheight = 0.25;
    float sample_halfwidth_x = (2.0 * 0.95) / 2;
    float sample_halfwidth_y = (1.0 * 0.95) / 2;
    float offset_z = bottom + sample_halfheight + 0.03;

    float3 sample_center = make_float3(0, 0, offset_z);
    auto terrain_positions = sampler.SampleBox(sample_center, make_float3(sample_halfwidth_x, sample_halfwidth_y, sample_halfheight));

    std::vector<std::shared_ptr<DEMSphereTemplate>> terrain_spheres(terrain_positions.size(), terrain_sphere_template);
    auto terrain_particles = DEMSim.AddSpheres(terrain_spheres, terrain_positions);

    // Optional: give particles a slight downward velocity to settle them quickly
    terrain_particles->SetVel(make_float3(0.0, 0.0, -0.05));

    std::cout << "Terrain particles added: " << terrain_positions.size() << std::endl;

    // Prescribed motion setup for the wheel family
    float w_r = math_PI / 4;  // Angular velocity (rad/s)
    float v_ref = w_r * wheel_radius;  // Linear velocity reference based on no-slip condition
    float total_pressure = 200.0;  // N
    float added_pressure = total_pressure - wheel_weight;

    DEMSim.SetFamilyPrescribedAngVel(1, "0", std::to_string(w_r), "0", false);
    DEMSim.AddFamilyPrescribedAcc(1, "none", "none", std::to_string(-added_pressure / wheel_mass));

    // Simulation configuration
    float step_size = 5e-6;
    DEMSim.SetInitTimeStep(step_size);
    DEMSim.SetGravitationalAcceleration(make_float3(0, 0, -G));
    DEMSim.SetMaxVelocity(20.);
    DEMSim.SetErrorOutVelocity(35.);
    DEMSim.SetExpandSafetyMultiplier(1.);
    DEMSim.SetCDUpdateFreq(40);
    DEMSim.DisableAdaptiveUpdateFreq();
    DEMSim.Initialize();

    // Initial wheel placement
    auto max_z_finder = DEMSim.CreateInspector("clump_max_z");
    float max_z = max_z_finder->GetValue();
    wheel_tracker->SetPos(make_float3(-0.45, 0, max_z + 0.03 + wheel_radius));

    // Sink-in phase (settling)
    unsigned int fps = 10;
    unsigned int out_steps = static_cast<unsigned int>(1.0 / (fps * step_size));
    unsigned int currframe = 0;
    double frame_time = 1.0 / fps;

    for (double t = 0; t < 1.0; t += frame_time) {
        std::cout << "Outputting frame (sink phase): " << currframe << std::endl;
        DEMSim.WriteSphereFile(out_dir / ("DEMdemo_output_" + std::to_string(currframe) + ".csv"));
        DEMSim.WriteMeshFile(out_dir / ("DEMdemo_mesh_" + std::to_string(currframe) + ".vtk"));
        DEMSim.DoDynamics(frame_time);
        currframe++;
    }

    // Switch to drawbar pull simulation phase
    DEMSim.DoDynamicsThenSync(0);
    DEMSim.ChangeFamily(1, 2);

    // Prescribe rolling with slip for the new family
    DEMSim.SetFamilyPrescribedAngVel(2, "0", std::to_string(w_r), "0", false);
    DEMSim.SetFamilyPrescribedLinVel(2, std::to_string(v_ref * 0.5), "0", "none", false);  // Slip ratio 0.5
    DEMSim.AddFamilyPrescribedAcc(2, "none", "none", std::to_string(-added_pressure / wheel_mass));

    // Main simulation loop
    double sim_end = 6.0;
    unsigned int curr_step = 0;
    unsigned int report_ps = 100;
    unsigned int report_steps = static_cast<unsigned int>(1.0 / (report_ps * step_size));

    auto max_v_finder = DEMSim.CreateInspector("clump_max_absv");

    auto start_time = std::chrono::high_resolution_clock::now();

    for (double t = 0; t < sim_end; t += step_size, curr_step++) {
        if (curr_step % out_steps == 0) {
            std::cout << "Outputting frame: " << currframe << std::endl;
            DEMSim.WriteSphereFile(out_dir / ("DEMdemo_output_" + std::to_string(currframe) + ".csv"));
            DEMSim.WriteMeshFile(out_dir / ("DEMdemo_mesh_" + std::to_string(currframe) + ".vtk"));
            currframe++;
        }

        if (curr_step % report_steps == 0) {
            float3 force = wheel_tracker->ContactAcc() * wheel_mass;
            std::cout << "Time: " << t
                      << " | Force on wheel: [" << force.x << ", " << force.y << ", " << force.z << "]"
                      << " | Drawbar pull coeff: " << force.x / total_pressure
                      << " | Max system velocity: " << max_v_finder->GetValue() << std::endl;
        }

        DEMSim.DoDynamics(step_size);
    }

    auto end_time = std::chrono::high_resolution_clock::now();
    auto elapsed_sec = std::chrono::duration<double>(end_time - start_time).count();

    std::cout << "Simulation finished in " << elapsed_sec << " seconds (wall time)." << std::endl;

    DEMSim.ShowTimingStats();
    DEMSim.ShowAnomalies();

    std::cout << "Custom wheel simulation complete." << std::endl;
    return 0;
}
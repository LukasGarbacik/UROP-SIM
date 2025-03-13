#convert all chrono output files to .vtk, so paraview can process whole visual
import pyvista as pv
import glob
import os
import pandas as pd
import numpy as np

# Create subdirectories for organization
converted_dir = 'converted'
particle_dir = os.path.join(converted_dir, 'particles')
mesh_dir = os.path.join(converted_dir, 'mesh')

# Create directories if they don't exist
for directory in [converted_dir, particle_dir, mesh_dir]:
    if not os.path.exists(directory):
        os.makedirs(directory)

# Get all CSV files in order (particle data)
csv_files = sorted(glob.glob('wheelOutput/*.csv'),
                  key=lambda x: int(x.split('_')[-1].split('.')[0]) if 'output' in x else -1)

# Get all VTK files in order (mesh data)
vtk_files = sorted(glob.glob('wheelOutput/*.vtk'),
                  key=lambda x: int(x.split('_')[-1].split('.')[0]))

# Convert each CSV file to VTK and organize in particle directory
for csv_file in csv_files:
    # Read particle data from CSV
    particles = pd.read_csv(csv_file)
    points = np.column_stack((particles['X'], particles['Y'], particles['Z']))
    particle_cloud = pv.PolyData(points)
    
    # point data
    particle_cloud['velocity'] = particles['absv']
    particle_cloud['radius'] = particles['r']
    
    base_name = os.path.basename(csv_file)
    vtk_filename = os.path.join(particle_dir, base_name.replace('.csv', '_particles.vtk'))
    
    particle_cloud.save(vtk_filename)

# Copy mesh VTK files to mesh directory
for vtk_file in vtk_files:
    # Generate destination path
    dest_file = os.path.join(mesh_dir, os.path.basename(vtk_file))
    
    # Copy the file
    with open(vtk_file, 'rb') as src, open(dest_file, 'wb') as dst:
        dst.write(src.read())

print(f"Conversion complete. Files organized in {converted_dir}:")
print(f"- Particle data: {particle_dir}")
print(f"- Mesh data: {mesh_dir}")
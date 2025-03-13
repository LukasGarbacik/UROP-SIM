import pyvista as pv
import glob
import os
import pandas as pd
import numpy as np

# Set up the plotter for generating HTML
pv.set_plot_theme('document')
pl = pv.Plotter(notebook=True, window_size=[1024, 768], off_screen=True)
pl.set_background('white')

# Enable anti-aliasing for better quality
pl.enable_anti_aliasing()

# Get all files in order
#object directory for wheel/ball over time
vtk_files = sorted(glob.glob('wheelOutput/*.vtk'), 
                  key=lambda x: int(x.split('_')[-1].split('.')[0]))
#particle directory for granular material over time
csv_files = sorted(glob.glob('wheelOutput/*.csv'), 
                  key=lambda x: int(x.split('_')[-1].split('.')[0]) if 'output' in x else -1)

# Process bed -> possible file for our own implimentation, but there will be no deined box right now
#bed_data = pd.read_csv(csv_files[0])
#bed_points = np.column_stack((bed_data['X'], bed_data['Y'], bed_data['Z']))
#bed_cloud = pv.PolyData(bed_points)
#csv_files.pop(0)

# debug comment for file length
print(f"DEBUG: Found {len(vtk_files)} VTK files and {len(csv_files)} CSV files")

# Set up angled camera position and lighting
pl.camera_position = [(2, 1.5, 1),  # Camera position at an angle
                     (0, 0, -0.5),   # Focal point moved down to align with box bottom
                     (0, 0, 1)]   # Up vector unchanged
pl.camera.zoom(0.84)  # Adjusted zoom for better view

# Adjust lighting for better depth perception at the new angle
pl.add_light(pv.Light(position=(2, 2, 2), focal_point=(0, 0, -0.5), intensity=0.7))
pl.add_light(pv.Light(position=(-2, -2, 2), focal_point=(0, 0, -0.5), intensity=0.3))

# Add a wireframe box for reference (matching world dimensions from wheel.cpp)
box = pv.Box(bounds=(-1, 1, -0.5, 0.5, -0.5, 1))
pl.add_mesh(box, style='wireframe', color='black', line_width=1)

# Add axis labels
pl.add_axes(xlabel='X (m)', ylabel='Y (m)', zlabel='Z (m)', 
           line_width=2, labels_off=False)

# Gif initialization
pl.open_gif('wheel_simulation.gif')


#global max velocity
maxV = 0.4 #auto rescales to max seen if this is surpassed by any particle

# Add each frame
for frame_num, (vtk_file, csv_file) in enumerate(zip(vtk_files, csv_files)):
    # Read mesh data for the wheel
    mesh = pv.read(vtk_file)
    
    # Read particle data from CSV
    particles = pd.read_csv(csv_file)
    
    # Create points for particles and set up velocity-based coloring
    points = np.column_stack((particles['X'], particles['Y'], particles['Z']))
    velocities = particles['absv']
    particle_cloud = pv.PolyData(points)
    particle_cloud['velocity'] = velocities
    if frame_num > 4 and maxV < velocities.max():
        maxV = velocities.max()
    
    pl.clear()
    
    # Add the box back after clearing
    pl.add_mesh(box, style='wireframe', color='black', line_width=1)
    
    # Add the particles with velocity-based coloring
    pl.add_mesh(particle_cloud, 
                scalars='velocity',
                cmap='coolwarm',
                point_size=12,
                opacity=0.9,
                render_points_as_spheres=True, 
                style='surface',
                edge_color='black',
                show_edges=True,
                line_width=2,
                lighting=True,
                ambient=0.3,
                diffuse=0.8,
                specular=0.5,
                clim=[0, maxV]) 
    
    # Add the wheel mesh with enhanced visibility
    pl.add_mesh(mesh, 
                color='steelblue',
                opacity=0.95,
                show_edges=True,
                edge_color='black',
                line_width=1,
                ambient=0.4,
                diffuse=0.6,
                specular=0.3,
                )
    
    # Gif add frame
    pl.write_frame()

    # Save frame as PNG
    #pl.screenshot(f'{output_dir}/frame_{frame_num:04d}.png', transparent_background=False)

pl.close()
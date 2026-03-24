import cupy as cp
import pyvista as pv
import numpy as np
import sys

# --- 1. Parameters ---
N = 1000
dt = 0.005
sub_steps = 4

# Lennard-Jones parameters
epsilon = 10.0      # Depth of the potential well (attractive strength)
sigma = 4.0         # Finite distance at which inter-particle potential is zero
mass = 1.0
softening = 1.5     # Prevents 1/r^12 explosion at extremely close ranges
max_force = 500.0   # Hard cap on repulsion

# Dynamic Globals
box_size = 30.0
target_temp = 50.0  # Initial target temperature
gravity = 10.0      # Downward acceleration (Y-axis)
tau = 0.1           # Thermostat coupling time constant

# --- 2. Initialization ---
def init_gas():
    # CuPy PRNG fixed! Generate uniformly directly on the GPU.
    pos = cp.random.uniform(-box_size * 0.8, box_size * 0.8, (N, 3), dtype=cp.float32)
    vel = cp.random.normal(0, cp.sqrt(target_temp), (N, 3), dtype=cp.float32)
    return pos, vel

pos, vel = init_gas()

# --- 3. PyVista Setup ---
plotter = pv.Plotter(title="VDW Thermodynamics (CuPy 11.6 Clean)")
plotter.set_background('#050510')

# Draw the bounding box
base_box = pv.Box(bounds=[-1, 1, -1, 1, -1, 1])
box_mesh = base_box.copy()
box_mesh.points *= box_size
box_actor = plotter.add_mesh(box_mesh, style='wireframe', color='white', opacity=0.2)

point_cloud = pv.PolyData(pos.get())

speeds = cp.linalg.norm(vel, axis=1).get()
point_cloud['Speed'] = speeds

plotter.add_mesh(
    point_cloud,
    scalars='Speed',
    cmap='plasma',
    point_size=8.0,
    render_points_as_spheres=True,
    show_scalar_bar=True,
    clim=[0, 40]
)
plotter.camera_position = [(0, 0, 100.0), (0, 0, 0), (0, 1, 0)]

# --- 4. Controls (VTK Sliders) ---
def set_temp(value):
    global target_temp
    target_temp = value

def set_volume(value):
    global box_size
    box_size = value
    box_mesh.points = base_box.points * box_size

plotter.add_slider_widget(set_temp, [0.1, 150.0], value=target_temp, title="Target Temperature", 
                          pointa=(0.025, 0.1), pointb=(0.3, 0.1), style='modern')
plotter.add_slider_widget(set_volume, [15.0, 50.0], value=box_size, title="Box Size (Volume)", 
                          pointa=(0.675, 0.1), pointb=(0.95, 0.1), style='modern')

# --- 5. The Hardware-Safe Engine (Explicit UI Loop) ---

def update_physics():
    global pos, vel
    
    for _ in range(sub_steps):
        x = pos[:, 0]
        y = pos[:, 1]
        z = pos[:, 2]
        
        diff_x = x[None, :] - x[:, None]
        diff_y = y[None, :] - y[:, None]
        diff_z = z[None, :] - z[:, None]
        
        dist_sq = diff_x**2 + diff_y**2 + diff_z**2
        dist_sq_safe = cp.maximum(dist_sq, softening**2)
        
        w = (sigma**2) / dist_sq_safe
        w3 = w * w * w
        w6 = w3 * w3
        
        coeff = (24.0 * epsilon / dist_sq_safe) * (2.0 * w6 - w3)
        coeff = cp.clip(coeff, -max_force, max_force)
        
        # Math is restored! Correct axis=1 summation.
        accel_x = -cp.sum(diff_x * coeff, axis=1) / mass
        accel_y = -cp.sum(diff_y * coeff, axis=1) / mass
        accel_z = -cp.sum(diff_z * coeff, axis=1) / mass
        
        # Apply downward Gravity to Y axis
        accel_y -= gravity
        
        accel = cp.stack([accel_x, accel_y, accel_z], axis=1)
        
        # Integration
        vel += accel * dt
        
        # Berendsen Thermostat heavily modified for Gravity Free-Fall:
        current_temp_xz = cp.mean(vel[:, 0]**2 + vel[:, 2]**2) / 2.0
        
        if current_temp_xz > 0:
            lmbda = cp.sqrt(target_temp / current_temp_xz)
            vel_scaling = 1.0 + (lmbda - 1.0) * (dt / tau)
            vel_scaling = cp.clip(vel_scaling, 0.8, 1.2)
            
            # Scale ONLY X and Z thermal velocities
            vel[:, 0] *= vel_scaling
            vel[:, 2] *= vel_scaling
            
        pos += vel * dt
        
        # Global Speed Hard-Cap to prevent overlapping LJ crystal explosions
        speed = cp.linalg.norm(vel, axis=1, keepdims=True)
        speed_safe = cp.maximum(speed, 1e-7)
        vel = cp.where(speed > 60.0, (vel / speed_safe) * 60.0, vel)
        
        # Wall Collisions (Inelastic bounce)
        for i in range(3):
            mask_max = pos[:, i] > box_size
            vel[mask_max, i] *= -0.8
            pos[mask_max, i] = box_size
            
            mask_min = pos[:, i] < -box_size
            vel[mask_min, i] *= -0.8
            pos[mask_min, i] = -box_size

    # Update visual points
    cpu_pos = pos.get()
    point_cloud.points = cpu_pos
    point_cloud['Speed'] = cp.linalg.norm(vel, axis=1).get()


# THE FIX: Tell PyVista you are taking over the animation loop manually!
plotter.show(interactive_update=True)
print("Van der Waals gas running. Explicit UI Loop active!")

while True:
    try:
        # 1. Math block (GPU runs physics)
        update_physics()
        
        # 2. UI block (GUI processes mouse/keyboard interactions, sliders, then draws the screen)
        plotter.update()
        
    except Exception as e:
        # Only reached if the window is closed or an error occurs (or Ctrl+C!)
        print("\nSimulation ended safely.")
        break

sys.exit(0)

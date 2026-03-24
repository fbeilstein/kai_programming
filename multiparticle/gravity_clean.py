import cupy as cp
import pyvista as pv
import numpy as np
import sys

# --- 1. Parameters ---
N_stars = 1500
G = 5.0            
softening = 3.0    
dt = 0.01          
sub_steps = 2      

M_blackhole = 2500.0   # The Supermassive Black Hole holding the galaxy together!

def create_galaxy(center, velocity, color_val):
    # --- The Stars ---
    r = cp.random.uniform(5, 30, N_stars, dtype=cp.float32)
    theta = cp.random.uniform(0, 2 * cp.pi, N_stars, dtype=cp.float32)
    z = cp.random.normal(0, 1.5, N_stars, dtype=cp.float32)
    
    # Calculate perfectly balanced velocity
    M_enclosed = M_blackhole + N_stars * cp.maximum((r - 5.0) / 25.0, 0.0)
    v_orbit = cp.sqrt(G * M_enclosed / r)
    
    pos_stars = cp.stack([
        r * cp.cos(theta) + center[0],
        r * cp.sin(theta) + center[1],
        z + center[2]
    ], axis=1)
    
    vel_stars = cp.stack([
        -cp.sin(theta) * v_orbit + velocity[0],
        cp.cos(theta) * v_orbit + velocity[1],
        cp.full(N_stars, velocity[2], dtype=cp.float32)
    ], axis=1)
    
    # --- The Central Supermassive Black Hole ---
    pos_bh = cp.array([center], dtype=cp.float32)
    vel_bh = cp.array([velocity], dtype=cp.float32)
    
    # Combine everything
    pos = cp.vstack([pos_bh, pos_stars])
    vel = cp.vstack([vel_bh, vel_stars])
    
    # Assign Mass (BH is extremely heavy, stars are 1.0)
    masses = cp.ones(N_stars + 1, dtype=cp.float32)
    masses[0] = M_blackhole
    
    # Assign Colors (Make the Black Holes bright yellow/white!)
    cols = cp.full(N_stars + 1, color_val, dtype=cp.float32)
    cols[0] = 5.0 
    
    return pos, vel, masses, cols

pos_a, vel_a, mass_a, col_a = create_galaxy([-35, -8, 0], [1.2, 0.2, 0], 0.2)
pos_b, vel_b, mass_b, col_b = create_galaxy([35, 8, 0], [-1.2, -0.2, 0], 0.8)

pos = cp.vstack([pos_a, pos_b])
vel = cp.vstack([vel_a, vel_b])
masses = cp.concatenate([mass_a, mass_b])
color_scalars = cp.concatenate([col_a, col_b]).get()

# --- 2. PyVista Setup ---
plotter = pv.Plotter(title="GPU Galaxy Collision (SMBH Engine)")
plotter.set_background('#020205') 

point_cloud = pv.PolyData(pos.get())
point_cloud['ID'] = color_scalars

plotter.add_mesh(
    point_cloud, 
    scalars='ID', 
    cmap='plasma', # Swapped to plasma to make the BHs glow beautifully
    point_size=5.0, 
    render_points_as_spheres=True, 
    show_scalar_bar=False
)
plotter.camera_position = [(0, 0, 180), (0, 0, 0), (0, 1, 0)]


# --- 3. The Hardware-Safe Engine ---
def update_physics():
    global pos, vel
    
    for _ in range(sub_steps):
        x = pos[:, 0]
        y = pos[:, 1]
        z = pos[:, 2]
        
        diff_x = x[None, :] - x[:, None]
        diff_y = y[None, :] - y[:, None]
        diff_z = z[None, :] - z[:, None]
        
        dist_sq = diff_x**2 + diff_y**2 + diff_z**2 + (softening**2)
        inv_dist3 = 1.0 / (dist_sq * cp.sqrt(dist_sq))
        
        # Multiply the gravitational pull by the *mass of the particle doing the pulling*!
        # This makes the black holes pull 2500x harder than a regular star.
        pull_strength = inv_dist3 * masses[None, :]
        
        accel_x = G * cp.sum(diff_x * pull_strength, axis=1)
        accel_y = G * cp.sum(diff_y * pull_strength, axis=1)
        accel_z = G * cp.sum(diff_z * pull_strength, axis=1)
        
        accel = cp.stack([accel_x, accel_y, accel_z], axis=1)
        
        vel += accel * dt
        pos += vel * dt
        
    point_cloud.points = pos.get()

plotter.show(interactive_update=True)
print("Black Hole physics injected! Starting simulation...")

while True:
    try:
        update_physics()
        plotter.update()
    except Exception:
        break

sys.exit(0)

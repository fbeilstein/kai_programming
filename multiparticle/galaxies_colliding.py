import numpy as np
import pyvista as pv
import time
import sys

# Instruct students without a GPU to change the next line to: import numpy as cp
import cupy as cp  

# =========================================================================
# CONFIGURATION & PHYSICS CONSTANTS
# =========================================================================
# DEBUG = True: Tiny galaxies (50 stars). Safe for testing Python loops.
# DEBUG = False: Massive galaxies (1500 stars). Python loops will grind to a 
# halt, proving the extreme necessity of GPU broadcasting!
DEBUG = True

N_stars = 25 if DEBUG else 1500
G = 5.0            
softening = 3.0    
dt = 0.01          
sub_steps = 2      

M_blackhole = 2500.0   

# =========================================================================
# SUBPROBLEM 1: COMPUTE ACCELERATIONS (The O(N^2) Bottleneck)
# Optimize this function removing cycles
# =========================================================================
def compute_accelerations(pos, masses):
    """
    Calculates the gravitational pull every star exerts on every other star.
    Input and Output must be standard arrays (NumPy for CPU, CuPy for GPU).
    """
    N = len(pos)
    
    # --- LOOPY PYTHON VERSION (SLOW) ---
    accel = np.zeros_like(pos)
    for i in range(N):
        for j in range(N):
            dx = pos[j, 0] - pos[i, 0]
            dy = pos[j, 1] - pos[i, 1]
            dz = pos[j, 2] - pos[i, 2]
            
            dist_sq = dx**2 + dy**2 + dz**2 + softening**2
            inv_dist3 = 1.0 / (dist_sq * np.sqrt(dist_sq))
            
            pull_strength = inv_dist3 * masses[j]
            
            accel[i, 0] += G * dx * pull_strength
            accel[i, 1] += G * dy * pull_strength
            accel[i, 2] += G * dz * pull_strength
            
    return accel

    # --- CUPY / VECTORIZED VERSION (FAST) ---
    # TODO: Comment out the slow loop above, and uncomment this!
    # Hint: Use matrix broadcasting to avoid loops entirely.
    # 
    # x, y, z = pos[:, 0], pos[:, 1], pos[:, 2]
    # 
    # # diff_x[i, j] represents the x-distance between particle i and j
    # diff_x = x[None, :] - x[:, None]
    # diff_y = y[None, :] - y[:, None]
    # diff_z = z[None, :] - z[:, None]
    # 
    # dist_sq = diff_x**2 + diff_y**2 + diff_z**2 + (softening**2)
    # inv_dist3 = 1.0 / (dist_sq * cp.sqrt(dist_sq))
    # 
    # pull_strength = inv_dist3 * masses[None, :]
    # 
    # accel_x = G * cp.sum(diff_x * pull_strength, axis=1)
    # accel_y = G * cp.sum(diff_y * pull_strength, axis=1)
    # accel_z = G * cp.sum(diff_z * pull_strength, axis=1)
    # 
    # return cp.stack([accel_x, accel_y, accel_z], axis=1)


# =========================================================================
# SUBPROBLEM 2: UPDATE KINEMATICS
# Optimize this function removing cycles
# =========================================================================
def update_kinematics(pos, vel, accel):
    """
    Applies the acceleration to the velocities, and velocities to positions.
    """
    N = len(pos)
    
    # --- LOOPY PYTHON VERSION (SLOW) ---
    new_vel = np.zeros_like(vel)
    new_pos = np.zeros_like(pos)
    for i in range(N):
        new_vel[i, 0] = vel[i, 0] + accel[i, 0] * dt
        new_vel[i, 1] = vel[i, 1] + accel[i, 1] * dt
        new_vel[i, 2] = vel[i, 2] + accel[i, 2] * dt
        
        new_pos[i, 0] = pos[i, 0] + new_vel[i, 0] * dt
        new_pos[i, 1] = pos[i, 1] + new_vel[i, 1] * dt
        new_pos[i, 2] = pos[i, 2] + new_vel[i, 2] * dt
        
    return new_pos, new_vel

    # --- CUPY / VECTORIZED VERSION (FAST) ---
    # TODO: Comment out the slow loop above, and uncomment this!
    # 
    # new_vel = vel + accel * dt
    # new_pos = pos + new_vel * dt
    # return new_pos, new_vel


# ==============================================================================
# SUBPROBLEM 3: THE MEMORY PIPELINE (CPU -> GPU -> CPU)
# This is the ONLY function using update_kinematics and compute_accelerations.
# Gain some more optimization by transferring data to GPU here,
# performing update_kinematics and compute_accelerations on GPU data, 
# and copying to CPU before return.
# NOTE: leave the cycle here UNOPTIMIZED
# ==============================================================================
def update_step(pos_cpu, vel_cpu, masses_cpu):
    """
    Orchestrates the physics pipeline over multiple sub-steps.
    """
    # --- CPU PIPELINE (SLOW) ---
    p, v = pos_cpu, vel_cpu
    for _ in range(sub_steps):
        a = compute_accelerations(p, masses_cpu)
        p, v = update_kinematics(p, v, a)
    return p, v

    # --- GPU PIPELINE (FAST) ---
    # TODO: Comment out the CPU pipeline above, and uncomment this!
    # By looping on the GPU, we do 2x the physics math but only pay the 
    # penalty of copying the data across the motherboard ONCE per frame!
    #
    # # 1. Push arrays to GPU
    # d_pos = cp.asarray(pos_cpu)
    # d_vel = cp.asarray(vel_cpu)
    # d_masses = cp.asarray(masses_cpu)
    # 
    # # 2. Run the simulation loop entirely on the device
    # for _ in range(sub_steps):
    #     d_accel = compute_accelerations(d_pos, d_masses)
    #     d_pos, d_vel = update_kinematics(d_pos, d_vel, d_accel)
    # 
    # # 3. Pull results back to CPU for rendering
    # return cp.asnumpy(d_pos), cp.asnumpy(d_vel)


# =========================================================================
# GALAXY INITIALIZATION (PURE NUMPY)
# =========================================================================
def create_galaxy(center, velocity, color_val):
    # Pure NumPy initialization to keep the CPU/GPU boundary clean
    r = np.random.uniform(5, 30, N_stars).astype(np.float32)
    theta = np.random.uniform(0, 2 * np.pi, N_stars).astype(np.float32)
    z = np.random.normal(0, 1.5, N_stars).astype(np.float32)
    
    M_enclosed = M_blackhole + N_stars * np.maximum((r - 5.0) / 25.0, 0.0)
    v_orbit = np.sqrt(G * M_enclosed / r)
    
    pos_stars = np.column_stack([
        r * np.cos(theta) + center[0],
        r * np.sin(theta) + center[1],
        z + center[2]
    ])
    
    vel_stars = np.column_stack([
        -np.sin(theta) * v_orbit + velocity[0],
        np.cos(theta) * v_orbit + velocity[1],
        np.full(N_stars, velocity[2], dtype=np.float32)
    ])
    
    pos_bh = np.array([center], dtype=np.float32)
    vel_bh = np.array([velocity], dtype=np.float32)
    
    pos = np.vstack([pos_bh, pos_stars])
    vel = np.vstack([vel_bh, vel_stars])
    
    masses = np.ones(N_stars + 1, dtype=np.float32)
    masses[0] = M_blackhole
    
    cols = np.full(N_stars + 1, color_val, dtype=np.float32)
    cols[0] = 0.5 
    
    return pos, vel, masses, cols


def main():
    print(f"Initializing {N_stars * 2 + 2} bodies...")
    pos_a, vel_a, mass_a, col_a = create_galaxy([-35, -8, 0], [1.2, 0.2, 0], 0.0) # Red
    pos_b, vel_b, mass_b, col_b = create_galaxy([35, 8, 0], [-1.2, -0.2, 0], 1.0) # Green

    pos = np.vstack([pos_a, pos_b])
    vel = np.vstack([vel_a, vel_b])
    masses = np.concatenate([mass_a, mass_b])
    colors = np.concatenate([col_a, col_b])

    plotter = pv.Plotter(title="GPU Galaxy Collision")
    plotter.set_background('#000000')

    point_cloud = pv.PolyData(pos)
    point_cloud['ID'] = colors

    plotter.add_mesh(
        point_cloud, 
        scalars='ID', 
        cmap='RdYlGn', 
        point_size=10.0 if DEBUG else 5.0, 
        render_points_as_spheres=True, 
        show_scalar_bar=False
    )
    plotter.camera_position = [(0, 0, 180), (0, 0, 0), (0, 1, 0)]
    

    plotter.show(interactive_update=True)

    print("\nSimulation Started!")
    if not DEBUG:
        print("WARNING: You are in DEBUG = False mode with loopy Python!")
        print("This will be unbelievably slow until you implement CuPy.")

    while True:
        try:
            start_time = time.time()
            
            # The Magic Pipeline Step
            pos, vel = update_step(pos, vel, masses)
            
            elapsed = time.time() - start_time
            
            # Update the UI
            point_cloud.points = pos
            
            # The modern PyVista way to update text: use the 'name' parameter!
            plotter.add_text(f"Step Time: {elapsed:.4f}s", position='upper_left', 
                             color='white', font_size=12, name='steptime')
            
            plotter.update()
            
        except Exception as e:
            # Never fly blind again! This will print the error if you click the 'X' 
            # or if the math crashes.
            print(f"Simulation stopped: {e}")
            break

if __name__ == '__main__':
    main()

import numpy as np
import pyvista as pv
import time
import sys

# Instruct students without a GPU to change the next line to: import numpy as cp
import cupy as cp  

# =========================================================================
# CONFIGURATION & PHYSICS CONSTANTS
# =========================================================================
# DEBUG = True: 100 particles. Safe for testing the intense Python loops.
# DEBUG = False: 1500 particles. The Python loop will freeze, demanding CuPy!
DEBUG = True

N = 10 if DEBUG else 1500
dt = 0.005
sub_steps = 4

# Lennard-Jones parameters
epsilon = 10.0      
sigma = 4.0         
mass = 1.0
softening = 1.5     
max_force = 500.0   

gravity = 10.0      
tau = 0.1           

# Dynamic Globals (Modified by UI Sliders)
current_box_size = 30.0
current_target_temp = 50.0  

# =========================================================================
# SUBPROBLEM 1: COMPUTE LENNARD-JONES FORCES (The O(N^2) Bottleneck)
# Optimize this function removing cycles
# =========================================================================
def compute_accelerations(pos):
    """
    Calculates the Lennard-Jones Van der Waals forces and applies gravity.
    Input and Output must be standard arrays (NumPy for CPU, CuPy for GPU).
    """
    # --- LOOPY PYTHON VERSION (SLOW) ---
    accel = np.zeros_like(pos)
    for i in range(N):
        for j in range(N):
            if i == j: continue
            
            dx = pos[j, 0] - pos[i, 0]
            dy = pos[j, 1] - pos[i, 1]
            dz = pos[j, 2] - pos[i, 2]
            
            dist_sq = dx**2 + dy**2 + dz**2
            dist_sq_safe = max(dist_sq, softening**2)
            
            w = (sigma**2) / dist_sq_safe
            w3 = w * w * w
            w6 = w3 * w3
            
            coeff = (24.0 * epsilon / dist_sq_safe) * (2.0 * w6 - w3)
            coeff = max(min(coeff, max_force), -max_force) # Clip force
            
            accel[i, 0] += -dx * coeff / mass
            accel[i, 1] += -dy * coeff / mass
            accel[i, 2] += -dz * coeff / mass
            
        # Apply downward gravity to the Y-axis
        accel[i, 1] -= gravity
        
    return accel

    # --- CUPY / VECTORIZED VERSION (FAST) ---
    # TODO: Comment out the slow loop above, and uncomment this!
    # 
    # x, y, z = pos[:, 0], pos[:, 1], pos[:, 2]
    # 
    # diff_x = x[None, :] - x[:, None]
    # diff_y = y[None, :] - y[:, None]
    # diff_z = z[None, :] - z[:, None]
    # 
    # dist_sq = diff_x**2 + diff_y**2 + diff_z**2
    # dist_sq_safe = cp.maximum(dist_sq, softening**2)
    # 
    # w = (sigma**2) / dist_sq_safe
    # w3 = w * w * w
    # w6 = w3 * w3
    # 
    # coeff = (24.0 * epsilon / dist_sq_safe) * (2.0 * w6 - w3)
    # coeff = cp.clip(coeff, -max_force, max_force)
    # 
    # accel_x = -cp.sum(diff_x * coeff, axis=1) / mass
    # accel_y = -cp.sum(diff_y * coeff, axis=1) / mass
    # accel_z = -cp.sum(diff_z * coeff, axis=1) / mass
    # 
    # accel_y -= gravity
    # 
    # return cp.stack([accel_x, accel_y, accel_z], axis=1)


# =========================================================================
# SUBPROBLEM 2: KINEMATICS, THERMOSTAT & BOUNDARIES
# Optimize this function removing cycles up to N
# NOTE: you can leave "small" cycles that are O(1)
# =========================================================================
def update_kinematics(pos, vel, accel, target_temp, box_size):
    """
    Integrates motion, scales velocity to match target temperature, 
    and bounces particles off the walls.
    """
    # --- LOOPY PYTHON VERSION (SLOW) ---
    new_pos = np.zeros_like(pos)
    new_vel = np.zeros_like(vel)
    
    # 1. Integrate Velocity
    for i in range(N):
        new_vel[i, 0] = vel[i, 0] + accel[i, 0] * dt
        new_vel[i, 1] = vel[i, 1] + accel[i, 1] * dt
        new_vel[i, 2] = vel[i, 2] + accel[i, 2] * dt
        
    # 2. Berendsen Thermostat (Only tracking X and Z thermal energy)
    current_temp_xz = 0.0
    for i in range(N):
        current_temp_xz += new_vel[i, 0]**2 + new_vel[i, 2]**2
    current_temp_xz = (current_temp_xz / N) / 2.0
    
    if current_temp_xz > 0:
        lmbda = np.sqrt(target_temp / current_temp_xz)
        vel_scaling = 1.0 + (lmbda - 1.0) * (dt / tau)
        vel_scaling = max(min(vel_scaling, 1.2), 0.8)
        for i in range(N):
            new_vel[i, 0] *= vel_scaling
            new_vel[i, 2] *= vel_scaling
            
    # 3. Integrate Position, Enforce Speed Limit, and Wall Bounce
    for i in range(N):
        speed = np.sqrt(new_vel[i, 0]**2 + new_vel[i, 1]**2 + new_vel[i, 2]**2)
        if speed > 60.0:
            new_vel[i] = (new_vel[i] / speed) * 60.0
            
        new_pos[i, 0] = pos[i, 0] + new_vel[i, 0] * dt
        new_pos[i, 1] = pos[i, 1] + new_vel[i, 1] * dt
        new_pos[i, 2] = pos[i, 2] + new_vel[i, 2] * dt
        
        for j in range(3):
            if new_pos[i, j] > box_size:
                new_vel[i, j] *= -0.8
                new_pos[i, j] = box_size
            elif new_pos[i, j] < -box_size:
                new_vel[i, j] *= -0.8
                new_pos[i, j] = -box_size
                
    return new_pos, new_vel

    # --- CUPY / VECTORIZED VERSION (FAST) ---
    # TODO: Comment out the slow loop above, and uncomment this!
    # 
    # vel += accel * dt
    # 
    # current_temp_xz = cp.mean(vel[:, 0]**2 + vel[:, 2]**2) / 2.0
    # if current_temp_xz > 0:
    #     lmbda = cp.sqrt(target_temp / current_temp_xz)
    #     vel_scaling = 1.0 + (lmbda - 1.0) * (dt / tau)
    #     vel_scaling = cp.clip(vel_scaling, 0.8, 1.2)
    #     vel[:, 0] *= vel_scaling
    #     vel[:, 2] *= vel_scaling
    # 
    # pos += vel * dt
    # 
    # speed = cp.linalg.norm(vel, axis=1, keepdims=True)
    # speed_safe = cp.maximum(speed, 1e-7)
    # vel = cp.where(speed > 60.0, (vel / speed_safe) * 60.0, vel)
    # 
    # for i in range(3):
    #     mask_max = pos[:, i] > box_size
    #     vel[mask_max, i] *= -0.8
    #     pos[mask_max, i] = box_size
    #     mask_min = pos[:, i] < -box_size
    #     vel[mask_min, i] *= -0.8
    #     pos[mask_min, i] = -box_size
    # 
    # return pos, vel


# =========================================================================
# SUBPROBLEM 3: THE MEMORY PIPELINE (CPU -> GPU -> CPU)
# This is the ONLY function using compute_accelerations and update_kinematics
# Gain some more optimization by transferring data to GPU here,
# performing compute_accelerations and update_kinematics on GPU data, 
# and copying to CPU before return
# NOTE: leave the cycle here UNOPTIMIZED
# =========================================================================
def update_step(pos_cpu, vel_cpu, target_temp, box_size):
    """
    Orchestrates the thermodynamics pipeline over multiple sub-steps.
    """
    # --- CPU PIPELINE (SLOW) ---
    p, v = pos_cpu, vel_cpu
    for _ in range(sub_steps):
        a = compute_accelerations(p)
        p, v = update_kinematics(p, v, a, target_temp, box_size)
    return p, v

    # --- GPU PIPELINE (FAST) ---
    # TODO: Comment out the CPU pipeline above, and uncomment this!
    #
    # # 1. Push arrays to GPU
    # d_pos = cp.asarray(pos_cpu)
    # d_vel = cp.asarray(vel_cpu)
    # 
    # # 2. Run the simulation loop entirely on the device
    # for _ in range(sub_steps):
    #     d_accel = compute_accelerations(d_pos)
    #     d_pos, d_vel = update_kinematics(d_pos, d_vel, d_accel, target_temp, box_size)
    # 
    # # 3. Pull results back to CPU for rendering
    # return cp.asnumpy(d_pos), cp.asnumpy(d_vel)


# =========================================================================
# UI ENGINE & INITIALIZATION (PURE NUMPY)
# =========================================================================
def init_gas():
    # Keep initialization in NumPy to maintain the CPU/GPU boundary
    pos = np.random.uniform(-current_box_size * 0.8, current_box_size * 0.8, (N, 3)).astype(np.float32)
    vel = np.random.normal(0, np.sqrt(current_target_temp), (N, 3)).astype(np.float32)
    return pos, vel

def main():
    global current_target_temp, current_box_size
    
    print(f"Initializing {N} gas molecules...")
    pos, vel = init_gas()

    plotter = pv.Plotter(title="VDW Thermodynamics")
    plotter.set_background('#050510')

    # Draw the bounding box
    base_box = pv.Box(bounds=[-1, 1, -1, 1, -1, 1])
    box_mesh = base_box.copy()
    box_mesh.points *= current_box_size
    plotter.add_mesh(box_mesh, style='wireframe', color='white', opacity=0.2)

    point_cloud = pv.PolyData(pos)
    speeds = np.linalg.norm(vel, axis=1)
    point_cloud['Speed'] = speeds

    plotter.add_mesh(
        point_cloud,
        scalars='Speed',
        cmap='plasma',
        point_size=12.0 if DEBUG else 6.0,
        render_points_as_spheres=True,
        show_scalar_bar=True,
        clim=[0, 40]
    )
    plotter.camera_position = [(0, 0, 100.0), (0, 0, 0), (0, 1, 0)]
    
    text_actor = plotter.add_text("Step Time: 0.000s", position='upper_left', color='white', font_size=12, name='steptime')

    # --- Controls (VTK Sliders) ---
    def set_temp(value):
        global current_target_temp
        current_target_temp = value

    def set_volume(value):
        global current_box_size
        current_box_size = value
        box_mesh.points = base_box.points * current_box_size

    plotter.add_slider_widget(set_temp, [0.1, 150.0], value=current_target_temp, title="Target Temperature", 
                              pointa=(0.025, 0.1), pointb=(0.3, 0.1), style='modern')
    plotter.add_slider_widget(set_volume, [15.0, 50.0], value=current_box_size, title="Box Size (Volume)", 
                              pointa=(0.675, 0.1), pointb=(0.95, 0.1), style='modern')

    plotter.show(interactive_update=True)

    print("\nSimulation Started!")
    if not DEBUG:
        print("WARNING: You are in DEBUG = False mode with loopy Python!")
        print("This will be unbelievably slow until you implement CuPy.")

    while True:
        try:
            start_time = time.time()
            
            # The Magic Pipeline Step
            pos, vel = update_step(pos, vel, current_target_temp, current_box_size)
            
            elapsed = time.time() - start_time
            
            # Update the UI
            point_cloud.points = pos
            point_cloud['Speed'] = np.linalg.norm(vel, axis=1)
            plotter.add_text(f"Step Time: {elapsed:.4f}s", position='upper_left', 
                             color='white', font_size=12, name='steptime')
            
            plotter.update()
            
        except Exception as e:
            print(f"\nSimulation stopped: {e}")
            break

if __name__ == '__main__':
    main()

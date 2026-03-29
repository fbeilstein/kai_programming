# pip install vispy
# pip install PyQt5

import numpy as np
import time
import sys

# Instruct students without a GPU to change the next line to: import numpy as cp
import cupy as cp

import vispy
vispy.use('PyQt5')

# VisPy for raw OpenGL rendering
from vispy import app, scene
from vispy.color import Colormap

# =========================================================================
# CONFIGURATION
# =========================================================================
# DEBUG = True: 100x100 grid. Safe for testing the slow Python FDTD loops.
# DEBUG = False: 500x500 grid. The Python loops will crawl, demanding CuPy!
DEBUG = True

GRID_SIZE = 100 if DEBUG else 500
sub_steps = 4

alpha = 0.25      
damping = 0.995 

# =========================================================================
# SUBPROBLEM 1: THE PHYSICS MATH
# Optimize this function removing cycles
# =========================================================================
def compute_laplacian(z):
    rows, cols = z.shape
    
    lap = np.zeros_like(z)
    for i in range(1, rows - 1):
        for j in range(1, cols - 1):
            lap[i, j] = (z[i-1, j] + z[i+1, j] + 
                         z[i, j-1] + z[i, j+1] - 
                         4.0 * z[i, j])
    return lap

    # --- CUPY / VECTORIZED VERSION (FAST) ---
    # TODO: Comment out the slow loop above, and uncomment this!
    # 
    # lap = cp.zeros_like(z)
    # lap[1:-1, 1:-1] = (z[:-2, 1:-1] + z[2:, 1:-1] + 
    #                    z[1:-1, :-2] + z[1:-1, 2:] - 
    #                    4.0 * z[1:-1, 1:-1])
    # return lap

# =========================================================================
# SUBPROBLEM 2: THE PHYSICS MATH
# Optimize this function removing cycles
# =========================================================================

def integrate_wave(z_prev, z_curr, laplacian):
    rows, cols = z_curr.shape
    
    # --- LOOPY PYTHON VERSION (SLOW) ---
    z_next = np.zeros_like(z_curr)
    for i in range(1, rows - 1):
        for j in range(1, cols - 1):
            z_next[i, j] = (2.0 * z_curr[i, j] - z_prev[i, j] + 
                            alpha * laplacian[i, j]) * damping
            
    return z_curr, z_next

    # --- CUPY / VECTORIZED VERSION (FAST) ---
    # TODO: Comment out the slow loop above, and uncomment this!
    # 
    # z_next = cp.zeros_like(z_curr)
    # z_next[1:-1, 1:-1] = (2.0 * z_curr[1:-1, 1:-1] - z_prev[1:-1, 1:-1] + 
    #                       alpha * laplacian[1:-1, 1:-1]) * damping
    # return z_curr, z_next


# =========================================================================
# SUBPROBLEM 3: THE MEMORY PIPELINE (CPU -> GPU -> CPU)
# This is the ONLY function using compute_laplacian and integrate_wave.
# Gain some more optimization by transferring data to GPU here,
# performing compute_laplacian and integrate_wave on GPU data, 
# and copying to CPU before return.
# NOTE: leave the cycle here UNOPTIMIZED
# =========================================================================
def update_step(z_prev_cpu, z_curr_cpu):
    # --- CPU PIPELINE (SLOW) ---
    p, c = z_prev_cpu, z_curr_cpu
    for _ in range(sub_steps):
        lap = compute_laplacian(c)
        p, c = integrate_wave(p, c, lap)
    return p, c

    # --- GPU PIPELINE (FAST) ---
    # TODO: Comment out the CPU pipeline above, and uncomment this!
    #
    # d_p = cp.asarray(z_prev_cpu)
    # d_c = cp.asarray(z_curr_cpu)
    # 
    # for _ in range(sub_steps):
    #     d_lap = compute_laplacian(d_c)
    #     d_p, d_c = integrate_wave(d_p, d_c, d_lap)
    # 
    # return cp.asnumpy(d_p), cp.asnumpy(d_c)


# =========================================================================
# UI ENGINE & INITIALIZATION (VISPY OPENGL)
# =========================================================================
def add_drop(z, cx, cy, radius, height):
    rows, cols = z.shape
    y, x = np.ogrid[-cy:rows-cy, -cx:cols-cx]
    mask = x**2 + y**2 <= radius**2
    dist = np.sqrt(x**2 + y**2)
    z[mask] += height * np.cos(np.pi / 2.0 * dist[mask] / radius)
    return z

def main():
    print(f"Initializing {GRID_SIZE}x{GRID_SIZE} OpenGL wave simulation...")
    
    z_prev = np.zeros((GRID_SIZE, GRID_SIZE), dtype=np.float32)
    z_curr = np.zeros((GRID_SIZE, GRID_SIZE), dtype=np.float32)

    # Prevent explosion on launch
    rad = 10 if DEBUG else 40
    z_curr = add_drop(z_curr, GRID_SIZE//2, GRID_SIZE//2, rad, 2.0)
    z_prev = np.copy(z_curr) 

    # --- Setup VisPy Canvas ---
    canvas = scene.SceneCanvas(keys='interactive', show=True, bgcolor='#050510', size=(1000, 800))
    canvas.title = "GPU 3D Wave Equation (VisPy)"
    view = canvas.central_widget.add_view()
    
    # Setup 3D Turntable Camera
    view.camera = scene.TurntableCamera(up='z', fov=45, distance=GRID_SIZE * 1.5, 
                                        elevation=35, azimuth=45)
    
    # Center the camera on the grid
    view.camera.center = (GRID_SIZE / 2, GRID_SIZE / 2, 0)

    # Create the 3D Surface
    # Notice how clean this is: it just takes the 2D array natively!
    surface = scene.visuals.SurfacePlot(z=z_curr, color=(0.3, 0.5, 1, 1), shading='smooth')
    view.add(surface)

    # Add text overlay
    # Add text overlay (Anchors set inside the constructor!)
    text = scene.visuals.Text("Initializing...", parent=canvas.scene, color='white', pos=(20, 30), font_size=14, anchor_x='left', anchor_y='center')

    # --- Interaction State ---
    raining = False

    @canvas.events.key_press.connect
    def on_key_press(event):
        nonlocal z_curr, z_prev, raining
        if event.text == ' ':
            drop_rad = 10 if DEBUG else 40
            z_curr = add_drop(z_curr, GRID_SIZE//2, GRID_SIZE//2, drop_rad, 2.0)
            z_prev = np.copy(z_curr)
        elif event.text.lower() == 'r':
            raining = not raining

    # --- The Physics Animation Loop ---
    def update(ev):
        nonlocal z_curr, z_prev
        start_time = time.time()
        
        if raining and np.random.rand() < 0.2:
            rx = np.random.randint(5, GRID_SIZE - 5)
            ry = np.random.randint(5, GRID_SIZE - 5)
            rad = 5 if DEBUG else 15
            z_curr = add_drop(z_curr, rx, ry, rad, 1.0)
            z_prev = add_drop(z_prev, rx, ry, rad, 1.0)
        
        z_prev, z_curr = update_step(z_prev, z_curr)
        
        # VisPy Magic: Just pass the 2D array. OpenGL handles the rest instantly.
        surface.set_data(z=z_curr)
        
        elapsed = time.time() - start_time
        fps = 1.0 / elapsed if elapsed > 0 else 0
        text.text = f"Grid: {GRID_SIZE}x{GRID_SIZE} | Time: {elapsed:.4f}s ({fps:.0f} FPS)\nSpace: Drop | R: Toggle Rain"

    # Start the hardware timer
    timer = app.Timer(interval='auto', connect=update, start=True)

    if not DEBUG:
        print("WARNING: You are in DEBUG = False mode with loopy Python!")
        print("The framerate will be terrible until you implement CuPy.")

    app.run()

if __name__ == '__main__':
    # You might need to 'pip install vispy' if you haven't already!
    main()

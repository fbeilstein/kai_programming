import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import time

# Instruct students without a GPU to change the next line to: import numpy as cp
import cupy as cp  

# =========================================================================
# CONFIGURATION
# =========================================================================
# DEBUG = True: 20x20 grid, starts empty. Perfect for step-by-step logic checking.
# DEBUG = False: 1024x1024 grid, starts random. Perfect for GPU benchmarking.
DEBUG = True

GRID_SIZE = 20 if DEBUG else 1024

# =========================================================================
# SUBPROBLEM 1: COUNT NEIGHBORS
# Optimize this function removing cycles
# =========================================================================
def count_neighbors(grid):
    """
    Calculates the number of alive neighbors for each cell in the grid.
    """
    rows, cols = grid.shape
    
    # --- LOOPY PYTHON VERSION (SLOW) ---
    neighbors = np.zeros_like(grid)
    for i in range(1, rows - 1):
        for j in range(1, cols - 1):
            neighbors[i, j] = (grid[i-1, j-1] + grid[i-1, j] + grid[i-1, j+1] +
                               grid[i, j-1]                  + grid[i, j+1] +
                               grid[i+1, j-1] + grid[i+1, j] + grid[i+1, j+1])
    return neighbors

    # --- CUPY / VECTORIZED VERSION (FAST) ---
    # TODO: Comment out the slow loop above, and uncomment this!
    # (Note: When using the GPU Pipeline in Subproblem 3, 'grid' will 
    # already be a CuPy array when it enters this function!)
    # 
    # neighbors = cp.zeros_like(grid)
    # N = (grid[:-2,:-2] + grid[:-2,1:-1] + grid[:-2,2:] +
    #      grid[1:-1,:-2]                 + grid[1:-1,2:] +
    #      grid[2:,:-2]  + grid[2:,1:-1]  + grid[2:,2:])
    # neighbors[1:-1, 1:-1] = N
    # return neighbors

# =========================================================================
# SUBPROBLEM 2: APPLY RULES
# Optimize this function removing cycles
# =========================================================================
def apply_rules(grid, neighbors):
    """
    Applies Conway's rules to determine the next state of the grid.
    """
    rows, cols = grid.shape
    
    # --- LOOPY PYTHON VERSION (SLOW) ---
    new_grid = np.zeros_like(grid)
    for i in range(1, rows - 1):
        for j in range(1, cols - 1):
            if grid[i, j] == 1:
                if neighbors[i, j] == 2 or neighbors[i, j] == 3:
                    new_grid[i, j] = 1
            else:
                if neighbors[i, j] == 3:
                    new_grid[i, j] = 1
    return new_grid

    # --- CUPY / VECTORIZED VERSION (FAST) ---
    # TODO: Comment out the slow loop above, and uncomment this!
    # 
    # new_grid = cp.zeros_like(grid)
    # inner_grid = grid[1:-1, 1:-1]
    # inner_neighbors = neighbors[1:-1, 1:-1]
    # 
    # birth = (inner_neighbors == 3) & (inner_grid == 0)
    # survive = ((inner_neighbors == 2) | (inner_neighbors == 3)) & (inner_grid == 1)
    # new_grid[1:-1, 1:-1][birth | survive] = 1
    # return new_grid

# =========================================================================
# SUBPROBLEM 3: THE MEMORY PIPELINE (CPU -> GPU -> CPU)
# This is the ONLY function using count_neighbors and apply_rules
# Gain some more optimization by transferring data to GPU here,
# performing count_neighbors and apply_rules on GPU data, and copying 
# to CPU before return
# =========================================================================
def update_step(grid):
    """
    Orchestrates one full step of the simulation.
    Input and Output must be standard NumPy arrays for the UI to render.
    """
    # --- CPU PIPELINE (SLOW) ---
    neighbors = count_neighbors(grid)
    new_grid = apply_rules(grid, neighbors)
    return new_grid

    # --- GPU PIPELINE (FAST) ---
    # TODO: Comment out the CPU pipeline above, and uncomment this!
    # This teaches a core rule of GPU programming: move the data across the 
    # PCI-e bus ONCE, do all your math on the device, and bring it back ONCE.
    #
    # # 1. Move data from CPU (NumPy) to GPU (CuPy)
    # d_grid = cp.asarray(grid)
    # 
    # # 2. Perform computations entirely on the GPU
    # d_neighbors = count_neighbors(d_grid)
    # d_new_grid = apply_rules(d_grid, d_neighbors)
    # 
    # # 3. Move data from GPU (CuPy) back to CPU (NumPy)
    # return cp.asnumpy(d_new_grid)


# =========================================================================
# THE UI ENGINE (PURE NUMPY)
# =========================================================================
def main():
    if DEBUG:
        grid = np.zeros((GRID_SIZE, GRID_SIZE), dtype=np.uint8)
    else:
        grid = np.random.choice([0, 1], size=(GRID_SIZE, GRID_SIZE), p=[0.8, 0.2]).astype(np.uint8)
        
    paused = True

    fig, ax = plt.subplots(figsize=(8, 8), facecolor='#1e1e1e')
    fig.canvas.manager.set_window_title('Interactive Game of Life')
    ax.axis('off')
    
    img = ax.imshow(grid, cmap='inferno', interpolation='nearest', vmin=0, vmax=1)
    
    def update_title(time_taken=None):
        if paused:
            ax.set_title(f"PAUSED [{GRID_SIZE}x{GRID_SIZE}] | Space: Play | Right Arrow: Step | C: Clear", color='white', pad=20)
        else:
            time_str = f" - Step Time: {time_taken:.4f}s" if time_taken else ""
            ax.set_title(f"RUNNING [{GRID_SIZE}x{GRID_SIZE}]{time_str} | Space: Pause", color='white', pad=20)
            
    update_title()

    def on_click(event):
        if event.inaxes != ax or event.xdata is None or event.ydata is None: return
        col, row = int(event.xdata + 0.5), int(event.ydata + 0.5)
        if 0 <= row < GRID_SIZE and 0 <= col < GRID_SIZE:
            grid[row, col] = 1 - grid[row, col]
            img.set_data(grid)
            fig.canvas.draw_idle()

    def on_drag(event):
        if event.inaxes != ax or not event.button or event.xdata is None or event.ydata is None: return
        col, row = int(event.xdata + 0.5), int(event.ydata + 0.5)
        if 0 <= row < GRID_SIZE and 0 <= col < GRID_SIZE:
            grid[row, col] = 1 if event.button == 1 else 0
            img.set_data(grid)
            fig.canvas.draw_idle()

    def on_key(event):
        nonlocal paused, grid
        if event.key == ' ':
            paused = not paused
            update_title()
        elif event.key == 'right': 
            if paused:
                grid = update_step(grid)
                img.set_data(grid)
        elif event.key == 'c':
            paused = True
            grid[:] = 0
            img.set_data(grid)
            update_title()
        elif event.key == 'r':
            grid[:] = np.random.choice([0, 1], size=(GRID_SIZE, GRID_SIZE), p=[0.8, 0.2]).astype(np.uint8)
            img.set_data(grid)
        fig.canvas.draw_idle()

    fig.canvas.mpl_connect('button_press_event', on_click)
    fig.canvas.mpl_connect('motion_notify_event', on_drag)
    fig.canvas.mpl_connect('key_press_event', on_key)

    def animate(frameNum):
        nonlocal grid
        if not paused:
            start_time = time.time()
            grid = update_step(grid)
            elapsed = time.time() - start_time
            update_title(elapsed)
            img.set_data(grid)
        return [img, ax.title]

    ani = animation.FuncAnimation(fig, animate, interval=50, blit=False, cache_frame_data=False)
    plt.tight_layout()
    plt.show()

if __name__ == '__main__':
    main()

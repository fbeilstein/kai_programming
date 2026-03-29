import numpy as np
import matplotlib.pyplot as plt
import time
import sys

# Instruct students without a GPU to change the next line to: import numpy as cp
import cupy as cp

# =========================================================================
# CONFIGURATION
# =========================================================================
# DEBUG = True: 100x100 resolution. The Python loops will take ~0.5s per frame.
# DEBUG = False: 1000x1000 resolution. The Python loop will completely freeze 
# your computer, proving why GPUs dominate graphics and fractals!
DEBUG = True

WIDTH = 100 if DEBUG else 1000
HEIGHT = 100 if DEBUG else 1000
INITIAL_MAX_ITER = 50 if DEBUG else 150

# =========================================================================
# SUBPROBLEM 1: COMPUTE MANDELBROT (The Escape Time Bottleneck)
# Optimize this function removing cycles over WIDTH and HEIGHT
# NOTE: cycle over max_iter may stay
# =========================================================================
def compute_mandelbrot(x_min, x_max, y_min, y_max, max_iter):
    """
    Calculates the Mandelbrot set using pure Python loops.
    """
    N = np.zeros((HEIGHT, WIDTH), dtype=np.int32)
    dx = (x_max - x_min) / (WIDTH - 1)
    dy = (y_max - y_min) / (HEIGHT - 1)

    for j in range(HEIGHT):
        y = y_min + j * dy
        for i in range(WIDTH):
            x = x_min + i * dx
            c = complex(x, y)
            z = 0j
            
            for k in range(max_iter):
                # Standard escape radius is 2.0
                if abs(z) > 2.0:
                    break
                z = z*z + c
                N[j, i] += 1
                
    return N

    # --- CUPY / VECTORIZED VERSION (FAST) ---
    # TODO: Implement the vectorized GPU version!
    # Hint 1: Use cp.linspace and cp.meshgrid to create a grid of X and Y coordinates.
    # Hint 2: Combine them into a complex array: C = X + 1j * Y
    # Hint 3: Use boolean masking (cp.abs(Z) <= 2.0) to only update points that haven't escaped.
    #
    # x = cp.linspace(x_min, x_max, WIDTH, dtype=cp.float32)
    # y = cp.linspace(y_min, y_max, HEIGHT, dtype=cp.float32)
    # X, Y = cp.meshgrid(x, y)
    # C = X + 1j * Y
    # 
    # Z = cp.zeros_like(C)
    # N = cp.zeros(C.shape, dtype=cp.int32)
    # 
    # for _ in range(max_iter):
    #     mask = cp.abs(Z) <= 2.0
    #     Z[mask] = Z[mask] * Z[mask] + C[mask]
    #     N[mask] += 1
    # 
    # return cp.asnumpy(N)



# =========================================================================
# UI ENGINE & INTERACTION (PURE NUMPY)
# =========================================================================
def update_step(bounds, max_iter):
    x_min, x_max, y_min, y_max = bounds
    return compute_mandelbrot(x_min, x_max, y_min, y_max, max_iter)


def main():
    # Initial physical boundaries of the complex plane
    bounds = [-2.0, 0.5, -1.25, 1.25] # x_min, x_max, y_min, y_max
    max_iter = INITIAL_MAX_ITER

    fig, ax = plt.subplots(figsize=(10, 10), facecolor='#1e1e1e')
    fig.canvas.manager.set_window_title('GPU Mandelbrot Zoomer')
    ax.axis('off')
    
    title = ax.set_title("Initializing Engine...", color='white', pad=20)
    
    # We use extent=bounds so matplotlib maps the mouse clicks to actual math coordinates
    img_display = ax.imshow(np.zeros((HEIGHT, WIDTH)), cmap='magma', extent=bounds, origin='lower', vmin=0)

    def update_view():
        try:
            t1 = time.perf_counter()
            mandel_data = update_step(bounds, max_iter)
            t2 = time.perf_counter()
            
            img_display.set_data(mandel_data)
            img_display.set_clim(0, max_iter)
            img_display.set_extent(bounds)
            
            title.set_text(f"Render Time: {(t2-t1)*1000:.1f}ms | Iter: {max_iter} | Left-Click: Zoom In | Right-Click: Zoom Out")
            fig.canvas.draw_idle()
            
        except Exception as e:
            print(f"Render Error: {e}")
            sys.exit(1)

    def on_click(event):
        nonlocal max_iter
        if event.inaxes != ax or event.xdata is None or event.ydata is None: return
        
        x_center, y_center = event.xdata, event.ydata
        dx = (bounds[1] - bounds[0])
        dy = (bounds[3] - bounds[2])
        
        # Left click zoom in (0.2x), Right click zoom out (5.0x)
        zoom_factor = 0.2 if event.button == 1 else 5.0 
        
        bounds[0] = x_center - (dx * zoom_factor / 2)
        bounds[1] = x_center + (dx * zoom_factor / 2)
        bounds[2] = y_center - (dy * zoom_factor / 2)
        bounds[3] = y_center + (dy * zoom_factor / 2)
        
        # Increase iteration depth to keep detail sharp as we zoom
        if event.button == 1:
            max_iter = int(max_iter * 1.2)
        else:
            max_iter = max(INITIAL_MAX_ITER, int(max_iter / 1.2))

        update_view()

    fig.canvas.mpl_connect('button_press_event', on_click)
    
    print("\nMandelbrot initialized! Click on the plot to zoom.")
    if not DEBUG:
        print("WARNING: DEBUG = False. Loopy Python will be brutally slow until you implement CuPy.")

    # Trigger first draw
    update_view()
    
    plt.tight_layout()
    plt.show()

if __name__ == '__main__':
    main()

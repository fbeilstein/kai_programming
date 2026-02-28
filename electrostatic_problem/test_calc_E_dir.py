import numpy as np
import matplotlib.pyplot as plt

# Import the student's function from their assignment file
from implementation_tasks import calc_E_dir

def plot_field(ax, charges_positions, charges_values, title):
    # Create a 2D grid of coordinates (Z=0)
    grid_res = 15
    x = np.linspace(-3, 3, grid_res)
    y = np.linspace(-3, 3, grid_res)
    X, Y = np.meshgrid(x, y)
    
    Ex = np.zeros_like(X)
    Ey = np.zeros_like(Y)
    
    # Format the charge data for the Numba function
    pos_arr = np.array(charges_positions, dtype=np.float64)
    q_arr = np.array(charges_values, dtype=np.float64)
    
    # Calculate the normalized electric field direction at each point
    for i in range(X.shape[0]):
        for j in range(X.shape[1]):
            # The student's function expects a 3D coordinate array
            pt = np.array([X[i, j], Y[i, j], 0.0], dtype=np.float64)
            
            # Extract the resulting 3D direction vector
            E_dir = calc_E_dir(pt, pos_arr, q_arr)
            
            # We only care about X and Y for the 2D plot
            Ex[i, j] = E_dir[0]
            Ey[i, j] = E_dir[1]
            
    # Draw the normalized directional arrows
    ax.quiver(X, Y, Ex, Ey, color='gray', pivot='middle')
    
    # Draw the charge sources as colored dots
    for pos, q in zip(charges_positions, charges_values):
        color = 'red' if q > 0 else 'blue'
        label = '+' if q > 0 else '-'
        ax.scatter(pos[0], pos[1], color=color, s=200, zorder=5)
        ax.text(pos[0], pos[1], label, color='white', ha='center', va='center', 
                fontweight='bold', fontsize=12, zorder=6)
        
    # Format the subplot
    ax.set_title(title, pad=10)
    ax.set_aspect('equal')
    ax.set_xlim(-3.5, 3.5)
    ax.set_ylim(-3.5, 3.5)
    ax.set_xticks([])
    ax.set_yticks([])

def run_tests():
    print("Running Field Direction Tests...")
    fig, axs = plt.subplots(2, 2, figsize=(10, 10))
    fig.canvas.manager.set_window_title('Task 1 Debugger: Electric Field Direction')
    
    # 1. Single Positive
    plot_field(axs[0, 0], [[0, 0, 0]], [1.0], "Single Positive Charge")
    
    # 2. Single Negative
    plot_field(axs[0, 1], [[0, 0, 0]], [-1.0], "Single Negative Charge")
    
    # 3. Two Positives (Repulsion)
    plot_field(axs[1, 0], [[-1.5, 0, 0], [1.5, 0, 0]], [1.0, 1.0], "Two Positive Charges")
    
    # 4. Dipole (Attraction)
    plot_field(axs[1, 1], [[-1.5, 0, 0], [1.5, 0, 0]], [1.0, -1.0], "Electric Dipole")
    
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    run_tests()

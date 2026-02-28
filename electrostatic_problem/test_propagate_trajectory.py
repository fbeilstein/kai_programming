import numpy as np
import matplotlib.pyplot as plt

# Import the student's function from their assignment file
from electrostatics_solver import propagate_trajectory

def test_rk2_integration():
    print("Running RK2 Trajectory Tests...")
    fig, ax = plt.subplots(figsize=(8, 8))
    fig.canvas.manager.set_window_title('Task 2 Debugger: RK2 Trajectory')
    
    # 1. Define an Electric Dipole
    charges_positions = np.array([[-1.5, 0.0, 0.0], [1.5, 0.0, 0.0]], dtype=np.float64)
    charges_values = np.array([1.0, -1.0], dtype=np.float64)
    
    max_steps = 200
    step_size = 0.1
    
    # Generate 16 seed points in a circle
    num_seeds = 16
    angles = np.linspace(0, 2 * np.pi, num_seeds, endpoint=False)
    radius = 0.4
    
    # 2. Plot trajectories from the Positive Charge (Trace Forward)
    for angle in angles:
        pt = np.array([-1.5 + radius * np.cos(angle), radius * np.sin(angle), 0.0], dtype=np.float64)
        
        # Pre-allocate memory for the student's function
        trajectory = np.zeros((max_steps, 3), dtype=np.float64)
        trajectory[0] = pt
        
        # Call the student's loop (current_charge_sign = 1.0)
        count = propagate_trajectory(
            trajectory, 1.0, charges_positions, charges_values, max_steps, step_size
        )
        
        # Slice the array to only plot the valid steps, then plot X vs Y
        valid_path = trajectory[:count]
        ax.plot(valid_path[:, 0], valid_path[:, 1], color='gold', linewidth=1.5)

    # 3. Plot trajectories from the Negative Charge (Trace Backward)
    for angle in angles:
        pt = np.array([1.5 + radius * np.cos(angle), radius * np.sin(angle), 0.0], dtype=np.float64)
        
        trajectory = np.zeros((max_steps, 3), dtype=np.float64)
        trajectory[0] = pt
        
        # Call the student's loop (current_charge_sign = -1.0)
        count = propagate_trajectory(
            trajectory, -1.0, charges_positions, charges_values, max_steps, step_size
        )
        
        valid_path = trajectory[:count]
        ax.plot(valid_path[:, 0], valid_path[:, 1], color='darkorange', linewidth=1.5, linestyle='--')

    # 4. Draw the source charges
    ax.scatter([-1.5], [0.0], color='red', s=300, zorder=5)
    ax.text(-1.5, 0.0, '+', color='white', ha='center', va='center', fontweight='bold', fontsize=14, zorder=6)
    
    ax.scatter([1.5], [0.0], color='blue', s=300, zorder=5)
    ax.text(1.5, 0.0, '-', color='white', ha='center', va='center', fontweight='bold', fontsize=14, zorder=6)
    
    # Format the plot 
    ax.set_title("RK2 Streamlines: Electric Dipole", pad=10)
    ax.set_aspect('equal')
    ax.set_xlim(-4, 4)
    ax.set_ylim(-4, 4)
    ax.grid(True, linestyle=':', alpha=0.6)
    plt.show()

if __name__ == "__main__":
    test_rk2_integration()

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# --- Configuration & Hyperparameters ---
G = -9.81         # Acceleration due to gravity (m/s^2)
DT = 0.05         # Time step (seconds)
BOUNCE_COEFF = 1.0  # Energy retained after a bounce (0 to 1)

# Initial State: [y_position, y_velocity]
initial_state = np.array([10.0, 0.0]) 

def update_state(state):
    """
    Students will implement this logic.
    Current state is a NumPy array: [position, velocity]
    Returns: Updated NumPy array: [new_position, new_velocity]
    """
    pos, vel = state

    # 1. Physics: Update velocity (v = v + a*dt)
    new_vel = vel + G * DT

    # 2. Physics: Update position (y = y + v*dt)
    new_pos = pos + new_vel * DT

    # 3. Collision Logic: If it hits the floor (y <= 0)
    if new_pos <= 0:
        new_pos = 0            # Reset to floor level
        new_vel = -new_vel * BOUNCE_COEFF  # Reverse and dampen velocity

    return np.array([new_pos, new_vel])

# --- Animation Boilerplate ---
fig, ax = plt.subplots()
ax.set_xlim(-1, 1)
ax.set_ylim(0, 12)
ax.set_title("Simple Bouncing Ball (Euler Method)")
ball, = ax.plot([], [], 'ro', ms=15)

current_state = initial_state.copy()

def init():
    ball.set_data([], [])
    return ball,

def animate(frame):
    global current_state
    current_state = update_state(current_state)
    
    # We only care about y-position for this 1D model, x is fixed at 0
    ball.set_data([0], [current_state[0]])
    return ball,

# Create animation
ani = FuncAnimation(fig, animate, frames=200, init_func=init, blit=True, interval=20)

plt.show()

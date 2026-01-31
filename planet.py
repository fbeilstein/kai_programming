import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# --- Simulation Constants ---
G = 1.0           # Gravitational constant (normalized for the toy model)
M_SUN = 100.0     # Mass of the sun
DT = 0.01         # Time step

# Initial State: [x, y, vx, vy]
# Placing planet at x=10, y=0, with a velocity primarily in the y-direction
initial_state = np.array([10.0, 0.0, 0.0, 3.2]) 

def update_state(state):
    """
    Logic for students to implement:
    1. Calculate distance (r) and the acceleration vector.
    2. Update velocity based on acceleration.
    3. Update position based on velocity.
    """
    x, y, vx, vy = state
    
    # 1. Calculate distance from Sun (origin)
    r = np.sqrt(x**2 + y**2)
    
    # 2. Calculate Acceleration (Inverse Square Law)
    # The direction is always toward the origin (-x, -y)
    accel_mag = G * M_SUN / r**2
    ax = -accel_mag * (x / r)
    ay = -accel_mag * (y / r)
    
    # 3. Update State using Euler Method
    new_vx = vx + ax * DT
    new_vy = vy + ay * DT
    new_x = x + new_vx * DT
    new_y = y + new_vy * DT
    
    return np.array([new_x, new_y, new_vx, new_vy])

# --- Animation Boilerplate ---
fig, ax = plt.subplots(figsize=(6, 6))
ax.set_xlim(-15, 15)
ax.set_ylim(-15, 15)
ax.set_aspect('equal') # Crucial so the orbit doesn't look stretched

# Draw the Sun
ax.plot(0, 0, 'yo', ms=15, label="Sun")
# Draw the Planet and its path
planet, = ax.plot([], [], 'bo', ms=8, label="Planet")
trail, = ax.plot([], [], 'b-', alpha=0.3) # To see the orbital path

state_history = [initial_state[0:2]]
current_state = initial_state.copy()

def animate(frame):
    global current_state
    current_state = update_state(current_state)
    state_history.append(current_state[0:2])
    
    # Update visual elements
    planet.set_data([current_state[0]], [current_state[1]])
    
    # Draw trail (limited to last 500 points for performance)
    history_arr = np.array(state_history[-500:])
    trail.set_data(history_arr[:, 0], history_arr[:, 1])
    
    return planet, trail

ani = FuncAnimation(fig, animate, frames=1000, interval=10, blit=True)
plt.legend()
plt.show()

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider

class CausticSimulation:
    def __init__(self, ax, radius=3.0):
        self.ax = ax
        self.radius = radius
        self.center = np.array([0.0, 0.0])
        
        # UI State
        self.dragging = False
        self.num_rays = 60
        
        # Visual Elements
        self.circle_patch = plt.Circle(self.center, self.radius, color='cyan', fill=False, lw=2, picker=True)
        self.ax.add_patch(self.circle_patch)
        self.ray_lines = [self.ax.plot([], [], 'r-', lw=0.5, alpha=0.5)[0] for _ in range(self.num_rays)]
        
        # Connect Events
        self.cid_press = ax.figure.canvas.mpl_connect('button_press_event', self.on_press)
        self.cid_move = ax.figure.canvas.mpl_connect('motion_notify_event', self.on_motion)
        self.cid_release = ax.figure.canvas.mpl_connect('button_release_event', self.on_release)

    def calculate_rays(self):
        """
        Isolated physics logic. 
        Rays start at x = -10 (infinity) and travel horizontally.
        """
        cx, cy = self.center
        # Generate rays covering the full height of the view
        y_starts = np.linspace(-10, 10, self.num_rays)
        ray_data = []

        for y in y_starts:
            # Check if ray hits the sphere: (y - cy)^2 must be < r^2
            dy = y - cy
            if abs(dy) < self.radius:
                # Solve for x-intercept: (x-cx)^2 + dy^2 = r^2
                dx_hit = -np.sqrt(self.radius**2 - dy**2)
                hit_x = cx + dx_hit
                hit_y = y
                
                # Normal vector at hit point
                nx, ny = (hit_x - cx) / self.radius, (hit_y - cy) / self.radius
                
                # Reflection: I = [1, 0]. R = I - 2(I·N)N
                dot = 1.0 * nx  # since I is [1, 0], dot is just nx
                rx = 1.0 - 2 * dot * nx
                ry = -2 * dot * ny
                
                # Ray path: Start -> Hit -> Reflected (extended to screen edge)
                ray_data.append(([ -10, hit_x, hit_x + rx * 20], [y, hit_y, hit_y + ry * 20]))
            else:
                # Ray misses: just a straight line across the screen
                ray_data.append(([-10, 10], [y, y]))
        return ray_data

    def update(self):
        self.circle_patch.center = self.center
        self.circle_patch.set_radius(self.radius)
        
        rays = self.calculate_rays()
        for line, data in zip(self.ray_lines, rays):
            line.set_data(data[0], data[1])
        self.ax.figure.canvas.draw_idle()

    # --- Event Handlers ---
    def on_press(self, event):
        if event.inaxes != self.ax: return
        contains, _ = self.circle_patch.contains(event)
        if contains:
            self.dragging = True
            self.offset = self.center - [event.xdata, event.ydata]

    def on_motion(self, event):
        if self.dragging and event.inaxes == self.ax:
            self.center = np.array([event.xdata, event.ydata]) + self.offset
            self.update()

    def on_release(self, event):
        self.dragging = False

# --- Setup Figure ---
fig, ax = plt.subplots(figsize=(8, 8))
plt.subplots_adjust(bottom=0.2)
ax.set_xlim(-10, 10)
ax.set_ylim(-10, 10)
ax.set_aspect('equal')
ax.set_title("Interactive Caustics: Drag Sphere / Use Slider")

sim = CausticSimulation(ax)

# Add Slider for Radius
ax_slider = plt.axes([0.2, 0.05, 0.6, 0.03])
r_slider = Slider(ax_slider, 'Radius', 0.5, 7.0, valinit=3.0)

def slide_update(val):
    sim.radius = val
    sim.update()

r_slider.on_changed(slide_update)

sim.update()
plt.show()

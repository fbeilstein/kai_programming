import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider

class RefractionSimulation:
    def __init__(self, ax, radius=3.0, n=1.33):
        self.ax = ax
        self.radius = radius
        self.n = n
        self.center = np.array([0.0, 0.0])
        self.num_rays = 60
        self.dragging = False
        
        # We use ONE line object for ALL rays. 
        # We separate them using np.nan so they don't connect visually.
        self.rays_plot, = self.ax.plot([], [], 'r-', lw=0.8, alpha=0.6)
        
        self.circle_patch = plt.Circle(self.center, self.radius, color='cyan', fill=False, lw=2)
        self.ax.add_patch(self.circle_patch)
        
        fig.canvas.mpl_connect('button_press_event', self.on_press)
        fig.canvas.mpl_connect('motion_notify_event', self.on_motion)
        fig.canvas.mpl_connect('button_release_event', self.on_release)

    def refract(self, I, N, n1, n2):
        eta = n1 / n2
        cos_i = -np.dot(N, I)
        sin2_t = eta**2 * (1 - cos_i**2)
        if sin2_t > 1.0: return None 
        cos_t = np.sqrt(1.0 - sin2_t)
        return eta * I + (eta * cos_i - cos_t) * N


    def calculate_rays(self):
            cx, cy = self.center
            y_starts = np.linspace(-9.5, 9.5, self.num_rays)
            
            all_x = []
            all_y = []

            for y in y_starts:
                # Default: Line from left to right
                px = [-10.0, 10.0]
                py = [float(y), float(y)]
                
                dy = y - cy
                if abs(dy) < self.radius:
                    dx_entry = -np.sqrt(self.radius**2 - dy**2)
                    p1 = np.array([cx + dx_entry, y])
                    n_entry = (p1 - self.center) / self.radius
                    v_in = np.array([1.0, 0.0])
                    v_mid = self.refract(v_in, n_entry, 1.0, self.n)
                    
                    if v_mid is not None:
                        cos_theta = -np.dot(v_mid, n_entry)
                        dist = 2 * cos_theta * self.radius
                        p2 = p1 + v_mid * dist
                        n_exit = (p2 - self.center) / self.radius
                        v_out = self.refract(v_mid, n_exit, self.n, 1.0)
                        
                        if v_out is not None:
                            # We use .tolist() to ensure we aren't nesting arrays
                            px = [-10.0, float(p1[0]), float(p2[0]), float(p2[0] + v_out[0] * 20)]
                            py = [float(y), float(p1[1]), float(p2[1]), float(p2[1] + v_out[1] * 20)]

                # Extend the master list with flattened values
                all_x.extend(px)
                all_x.append(np.nan)
                all_y.extend(py)
                all_y.append(np.nan)
                
            return np.array(all_x, dtype=float), np.array(all_y, dtype=float)


    def update(self):
        self.circle_patch.center = self.center
        self.circle_patch.set_radius(self.radius)
        
        x_data, y_data = self.calculate_rays()
        self.rays_plot.set_data(x_data, y_data)
        
        self.ax.figure.canvas.draw_idle()

    def on_press(self, event):
        if event.inaxes == self.ax and self.circle_patch.contains(event)[0]:
            self.dragging = True
            self.offset = self.center - [event.xdata, event.ydata]
    def on_motion(self, event):
        if self.dragging and event.inaxes == self.ax:
            self.center = np.array([event.xdata, event.ydata]) + self.offset
            self.update()
    def on_release(self, event): self.dragging = False

# --- Setup ---
fig, ax = plt.subplots(figsize=(8, 8))
plt.subplots_adjust(bottom=0.25)
ax.set_xlim(-10, 10); ax.set_ylim(-10, 10); ax.set_aspect('equal')

sim = RefractionSimulation(ax)

ax_n = plt.axes([0.2, 0.05, 0.6, 0.03])
ax_r = plt.axes([0.2, 0.10, 0.6, 0.03])
n_slider = Slider(ax_n, 'Index (n)', 1.0, 2.5, valinit=1.33)
r_slider = Slider(ax_r, 'Radius', 0.5, 7.0, valinit=3.0)

def ui_update(val):
    sim.n = n_slider.val
    sim.radius = r_slider.val
    sim.update()

n_slider.on_changed(ui_update)
r_slider.on_changed(ui_update)

sim.update()
plt.show()

import tkinter as tk
from tkinter import ttk, messagebox
import unittest
import numpy as np
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.patches import Arc, Circle

import implementation_tasks as tasks
import test_suite

class OpticsDebugger(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Optics Engine Debugger")
        self.geometry("1100x750")
        self.configure(bg="#1e1e1e")

        self.setup_sidebar()
        self.setup_main_area()
        self.refresh_tests()

    def setup_sidebar(self):
        # Using standard tk.Frame: removed 'padding', replaced with container padx/pady
        self.sidebar = tk.Frame(self, bg="#252526", highlightbackground="#333333", highlightthickness=1)
        self.sidebar.pack(side="left", fill="y", padx=0, pady=0)
        
        tk.Label(self.sidebar, text="Tasks", font=("Arial", 12, "bold"), 
                 bg="#252526", fg="white").pack(pady=15, padx=20)
        
        self.tasks = [
            (1, "Infinite Line"), (2, "Segment Bounds"), (3, "Circle Math"),
            (4, "Arc Sector"), (5, "Line Normals"), (6, "Arc Normals"), (7, "Refraction")
        ]
        
        self.status_indicators = {}
        for num, name in self.tasks:
            frame = tk.Frame(self.sidebar, bg="#252526")
            frame.pack(fill="x", pady=4, padx=10)
            
            ind = tk.Canvas(frame, width=15, height=15, highlightthickness=0, bg="#252526")
            ind.pack(side="left", padx=5)
            light = ind.create_oval(2, 2, 13, 13, fill="gray")
            self.status_indicators[num] = (ind, light)
            
            tk.Button(frame, text=f"L{num} {name}", width=20, font=("Arial", 9),
                       command=lambda n=num: self.switch_sandbox(n)).pack(side="left")

        tk.Button(self.sidebar, text="🔄 Reload & Retest", bg="#3e3e42", fg="white",
                  command=self.refresh_tests).pack(pady=30, padx=20, fill="x")

    def setup_main_area(self):
        self.main_container = tk.Frame(self, bg="#1e1e1e")
        self.main_container.pack(side="right", fill="both", expand=True)

        self.fig = Figure(figsize=(7, 7), facecolor="#252526")
        self.ax = self.fig.add_subplot(111)
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.main_container)
        self.canvas.get_tk_widget().pack(fill="both", expand=True)
        
        self.current_lvl = 1
        self.handles = {
            'ray_o': np.array([-30.0, 0.0]), 'ray_t': np.array([-10.0, 10.0]),
            'p1': np.array([10.0, -20.0]), 'p2': np.array([10.0, 20.0]),
            'center': np.array([10.0, 0.0]), 'radius_h': np.array([25.0, 5.0])
        }
        self.dragging = None
        
        self.canvas.mpl_connect('button_press_event', self.on_press)
        self.canvas.mpl_connect('button_release_event', self.on_release)
        self.canvas.mpl_connect('motion_notify_event', self.on_motion)

    def refresh_tests(self):
        print("\n>>> Running Optics Test Suite...")
        loader = unittest.TestLoader()
        suite = loader.loadTestsFromTestCase(test_suite.TestOpticsMath)
        
        for num, _ in self.tasks:
            level_suite = unittest.TestSuite([t for t in suite if f"test_l{num}" in t._testMethodName])
            res = unittest.TextTestRunner(verbosity=0).run(level_suite)
            
            canvas, light = self.status_indicators[num]
            if res.wasSuccessful() and res.testsRun > 0:
                canvas.itemconfig(light, fill="#4ec9b0") 
            elif res.testsRun > 0:
                canvas.itemconfig(light, fill="#f44747")
            else:
                canvas.itemconfig(light, fill="gray")
        self.update_plot()

    def switch_sandbox(self, level):
        self.current_lvl = level
        self.update_plot()

    def update_plot(self):
        self.ax.clear()
        self.ax.set_facecolor("#1e1e1e")
        self.ax.set_xlim(-60, 60); self.ax.set_ylim(-50, 50)
        self.ax.grid(True, color="#333333", linestyle='--')
        
        r_o, r_t = self.handles['ray_o'], self.handles['ray_t']
        p1, p2 = self.handles['p1'], self.handles['p2']
        center = self.handles['center']
        ray_dir = (r_t - r_o) / np.linalg.norm(r_t - r_o)
        
        SURF_COLOR = "#00ffff" # Cyan
        RAY_COLOR = "#ff8c00"  # Orange
        HANDLE_COLOR = "#ffff00" # Yellow
        
        t = float('inf')
        curve = None

        # Level 1, 2, 5, 7: Line-based surfaces
        if self.current_lvl in [1, 2, 5, 7]:
            slope_vec = p2 - p1
            line_pts = np.array([p1 - 100*slope_vec, p1 + 100*slope_vec])
            self.ax.plot(line_pts[:,0], line_pts[:,1], color="#444444", linestyle=":", alpha=0.6)
            self.ax.plot([p1[0], p2[0]], [p1[1], p2[1]], color=SURF_COLOR, lw=3, zorder=3)
            
            if self.current_lvl == 1:
                t, _ = tasks.intersect_line_infinite(r_o, ray_dir, p1, p2)
            else:
                t = tasks.intersect_segment(r_o, ray_dir, p1, p2)
            curve = {'type': 'line', 'p1': p1, 'p2': p2}
            
        # Level 3, 4, 6: Curve-based surfaces
        elif self.current_lvl in [3, 4, 6]:
            radius = np.linalg.norm(self.handles['radius_h'] - center)
            axis = np.array([-1, 0]); cos_half = 0.7 
            
            if self.current_lvl == 3:
                ts = tasks.intersect_circle_infinite(r_o, ray_dir, center, radius)
                t = min([val for val in ts if val > 1e-4], default=float('inf'))
                self.ax.add_patch(Circle(center, radius, color=SURF_COLOR, fill=False, lw=2, alpha=0.4))
            else:
                t = tasks.intersect_arc(r_o, ray_dir, center, radius, axis, cos_half)
                angle = np.degrees(np.arctan2(axis[1], axis[0]))
                span = np.degrees(np.arccos(cos_half)) * 2
                self.ax.add_patch(Arc(center, radius*2, radius*2, angle=angle, 
                                      theta1=-span/2, theta2=span/2, color=SURF_COLOR, lw=3))
            curve = {'type': 'arc', 'center': center}

        # Draw Ray
        hit_dist = t if (t != float('inf') and t > 0) else 150
        hit_pt = r_o + hit_dist * ray_dir
        self.ax.plot([r_o[0], hit_pt[0]], [r_o[1], hit_pt[1]], color=RAY_COLOR, lw=2, zorder=2)
        
        # Level 7 Refraction Demo
        if self.current_lvl == 7 and t != float('inf') and t > 0:
            norm = tasks.calculate_normal(hit_pt, ray_dir, curve)
            refr_dir = tasks.refract_vector(ray_dir, norm, 1.0, 1.5)
            if refr_dir is not None:
                refr_end = hit_pt + refr_dir * 40
                self.ax.plot([hit_pt[0], refr_end[0]], [hit_pt[1], refr_end[1]], color="#00ff00", lw=2)

        # Draw relevant Handles
        visible_handles = ['ray_o', 'ray_t']
        if self.current_lvl in [1, 2, 5, 7]: visible_handles += ['p1', 'p2']
        if self.current_lvl in [3, 4, 6]: visible_handles += ['center', 'radius_h']
        
        for k in visible_handles:
            v = self.handles[k]
            self.ax.scatter(v[0], v[1], c=HANDLE_COLOR, s=60, edgecolors='black', zorder=5)

        # Draw Intersection and Normal
        if t != float('inf') and t > 0:
            self.ax.scatter(hit_pt[0], hit_pt[1], c="#4ec9b0", s=120, edgecolors='white', zorder=10)
            if self.current_lvl in [5, 6, 7]:
                norm = tasks.calculate_normal(hit_pt, ray_dir, curve)
                self.ax.quiver(hit_pt[0], hit_pt[1], norm[0], norm[1], color="#c586c0", 
                               scale=12, width=0.012, pivot='tail', zorder=15)

        self.canvas.draw()

    def on_press(self, event):
        if event.xdata is None: return
        visible_keys = ['ray_o', 'ray_t']
        if self.current_lvl in [1, 2, 5, 7]: visible_keys += ['p1', 'p2']
        if self.current_lvl in [3, 4, 6]: visible_keys += ['center', 'radius_h']
        
        for k in visible_keys:
            v = self.handles[k]
            if np.hypot(event.xdata - v[0], event.ydata - v[1]) < 4:
                self.dragging = k; break

    def on_motion(self, event):
        if self.dragging and event.xdata is not None:
            self.handles[self.dragging] = np.array([event.xdata, event.ydata])
            self.update_plot()

    def on_release(self, event): self.dragging = None

if __name__ == "__main__":
    OpticsDebugger().mainloop()

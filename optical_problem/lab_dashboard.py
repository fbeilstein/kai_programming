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
        self.title("Optics Engine Debugger Pro")
        self.geometry("1150x850")
        self.configure(bg="#1e1e1e")

        self.setup_sidebar()
        self.setup_main_area()
        self.refresh_tests()

    def setup_sidebar(self):
        self.sidebar = tk.Frame(self, bg="#252526", highlightbackground="#333333", highlightthickness=1)
        self.sidebar.pack(side="left", fill="y", padx=0, pady=0)
        
        tk.Label(self.sidebar, text="Implementation Tasks", font=("Arial", 12, "bold"), 
                 bg="#252526", fg="white").pack(pady=15, padx=20)
        
        self.tasks = [
            (1, "Infinite Line"), (2, "Segment Bounds"), (3, "Circle Math"),
            (4, "Arc Sector"), (5, "Segment Normals"), (6, "Arc Normals"), (7, "Refraction")
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
                  command=self.refresh_tests).pack(pady=(30, 10), padx=20, fill="x")

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
            'center': np.array([15.0, 0.0]), 'arc_p1': np.array([0.0, 15.0]),
            'arc_p2': np.array([0.0, -15.0]), 'radius_h': np.array([30.0, 0.0])
        }
        self.dragging = None
        
        self.canvas.mpl_connect('button_press_event', self.on_press)
        self.canvas.mpl_connect('button_release_event', self.on_release)
        self.canvas.mpl_connect('motion_notify_event', self.on_motion)

    def refresh_tests(self):
        loader = unittest.TestLoader()
        suite = loader.loadTestsFromTestCase(test_suite.TestOpticsMath)
        for num, _ in self.tasks:
            level_suite = unittest.TestSuite([t for t in suite if f"test_l{num}" in t._testMethodName])
            res = unittest.TextTestRunner(verbosity=0).run(level_suite)
            canvas, light = self.status_indicators[num]
            color = "#4ec9b0" if (res.wasSuccessful() and res.testsRun > 0) else "#f44747" if res.testsRun > 0 else "gray"
            canvas.itemconfig(light, fill=color)
        self.update_plot()

    def update_plot(self):
        self.ax.clear(); self.ax.set_facecolor("#1e1e1e"); self.ax.set_xlim(-60, 60); self.ax.set_ylim(-50, 50); self.ax.set_aspect('equal'); self.ax.grid(True, color="#333333", linestyle='--')
        r_o, r_t = self.handles['ray_o'], self.handles['ray_t']
        ray_dir = (r_t - r_o) / np.linalg.norm(r_t - r_o)
        GUIDE_COLOR, SURF_COLOR, RAY_COLOR, HANDLE_COLOR = "#ff0000", "#00ffff", "#ff8c00", "#ffff00"
        
        t, curve = float('inf'), None
        try:
            if self.current_lvl in [1, 2, 5, 7]:
                p1, p2 = self.handles['p1'], self.handles['p2']
                slope = p2 - p1
                self.ax.plot([p1[0]-100*slope[0], p1[0]+100*slope[0]], [p1[1]-100*slope[1], p1[1]+100*slope[1]], color=GUIDE_COLOR, linestyle="--", alpha=0.4)
                self.ax.plot([p1[0], p2[0]], [p1[1], p2[1]], color=SURF_COLOR, lw=3)
                t = tasks.intersect_curve(r_o, ray_dir, {'type': 'line', 'p1': p1, 'p2': p2, 'medium_id': 1})
                curve = {'type': 'line', 'p1': p1, 'p2': p2}
            elif self.current_lvl in [3, 4, 6]:
                center = self.handles['center']; radius = np.linalg.norm(self.handles['arc_p1'] - center)
                self.ax.add_patch(Circle(center, radius, color=GUIDE_COLOR, fill=False, lw=1, linestyle="--", alpha=0.3))
                if self.current_lvl == 3:
                    self.ax.add_patch(Circle(center, radius, color=SURF_COLOR, fill=False, lw=2, alpha=0.3))
                else:
                    mid = (self.handles['arc_p1'] + self.handles['arc_p2']) / 2.0; axis = (mid - center); axis /= np.linalg.norm(axis)
                    cos_half = np.dot((self.handles['arc_p1'] - center) / radius, axis)
                    span = np.degrees(np.arccos(np.clip(cos_half, -1, 1))) * 2
                    self.ax.add_patch(Arc(center, radius*2, radius*2, angle=np.degrees(np.arctan2(axis[1], axis[0])), theta1=-span/2, theta2=span/2, color=SURF_COLOR, lw=3))
                t = tasks.intersect_curve(r_o, ray_dir, {'type': 'arc', 'center': center, 'radius': radius, 'axis': np.array([-1, 0]), 'cos_half_angle': 0.7, 'medium_id': 1})
                curve = {'type': 'arc', 'center': center}
        except Exception: t = float('inf') # Silent fail for plotting

        # Ray drawing (Stops at hit, or goes through if t is infinity)
        hit_dist = t if (isinstance(t, (float, int)) and t != float('inf') and t > 0) else 150
        hit_pt = r_o + hit_dist * ray_dir
        self.ax.plot([r_o[0], hit_pt[0]], [r_o[1], hit_pt[1]], color=RAY_COLOR, lw=2)

        if self.current_lvl == 7 and t != float('inf') and t > 0:
            norm = tasks.calculate_normal(hit_pt, ray_dir, curve)
            if np.linalg.norm(norm) > 1e-3:
                refr_dir = tasks.refract_vector(ray_dir, norm, 1.0, 1.5)
                if refr_dir is not None: self.ax.plot([hit_pt[0], hit_pt[0] + refr_dir[0]*40], [hit_pt[1], hit_pt[1] + refr_dir[1]*40], color="#00ff00", lw=2)

        # Draw Intersection Point & Normal Arrow
        if isinstance(t, (float, int)) and t != float('inf') and t > 0:
            self.ax.scatter(hit_pt[0], hit_pt[1], c="#4ec9b0", s=120, edgecolors='white', zorder=10)
            if self.current_lvl in [5, 6, 7]:
                norm = tasks.calculate_normal(hit_pt, ray_dir, curve)
                if np.linalg.norm(norm) > 1e-3: self.ax.quiver(hit_pt[0], hit_pt[1], norm[0], norm[1], color="#c586c0", scale=12, pivot='tail')

        # Visible Handles
        vis = ['ray_o', 'ray_t']
        if self.current_lvl in [1, 2, 5, 7]: vis += ['p1', 'p2']
        elif self.current_lvl == 3: vis += ['center', 'arc_p1'] # radius controlled by p1
        elif self.current_lvl in [4, 6]: vis += ['center', 'arc_p1', 'arc_p2']
        for k in vis: self.ax.scatter(self.handles[k][0], self.handles[k][1], c=HANDLE_COLOR, s=60, edgecolors='black', zorder=5)
        self.canvas.draw()

    # Interaction methods (on_motion, on_press, on_release) same as previous...
    def switch_sandbox(self, l): self.current_lvl = l; self.update_plot()
    def on_press(self, e):
        if e.xdata is None: return
        for k, v in self.handles.items():
            if np.hypot(e.xdata - v[0], e.ydata - v[1]) < 4: self.dragging = k; break
    def on_release(self, e): self.dragging = None
    def on_motion(self, e):
        if self.dragging and e.xdata is not None:
            new = np.array([e.xdata, e.ydata])
            if self.dragging == 'center':
                d1, d2 = self.handles['arc_p1'] - self.handles['center'], self.handles['arc_p2'] - self.handles['center']
                self.handles['center'] = new; self.handles['arc_p1'], self.handles['arc_p2'] = new + d1, new + d2
            elif self.dragging == 'arc_p1':
                self.handles['arc_p1'] = new
                r = np.linalg.norm(new - self.handles['center'])
                v2 = self.handles['arc_p2'] - self.handles['center']
                self.handles['arc_p2'] = self.handles['center'] + r * (v2 / np.linalg.norm(v2))
            elif self.dragging == 'arc_p2':
                r = np.linalg.norm(self.handles['arc_p1'] - self.handles['center'])
                v2 = new - self.handles['center']
                self.handles['arc_p2'] = self.handles['center'] + r * (v2 / np.linalg.norm(v2))
            else: self.handles[self.dragging] = new
            self.update_plot()

if __name__ == "__main__":
    OpticsDebugger().mainloop()
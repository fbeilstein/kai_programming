import tkinter as tk
from tkinter import ttk, messagebox
import unittest
import numpy as np
import subprocess
import sys
import importlib
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.patches import Arc, Circle

# Ensure these modules exist in your directory
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
                  command=self.refresh_tests).pack(pady=(30, 10), padx=20, fill="x")
        
        tk.Button(self.sidebar, text="🚀 Run Main Simulation", bg="#007acc", fg="white", font=("Arial", 10, "bold"),
                  command=self.run_main_script).pack(pady=10, padx=20, fill="x")

    def run_main_script(self):
        try:
            subprocess.Popen([sys.executable, "optic_bench.py"])
        except Exception as e:
            messagebox.showerror("Execution Error", f"Could not launch optic_bench.py: {e}")

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
            'center': np.array([15.0, 0.0]), 
            'radius_h': np.array([30.0, 0.0]), # Unique to L3
            'arc_p1': np.array([0.0, 15.0]),   # Unique to L4/L6
            'arc_p2': np.array([0.0, -15.0])   # Unique to L4/L6
        }
        self.dragging = None
        
        self.canvas.mpl_connect('button_press_event', self.on_press)
        self.canvas.mpl_connect('button_release_event', self.on_release)
        self.canvas.mpl_connect('motion_notify_event', self.on_motion)

    def refresh_tests(self):
        importlib.reload(tasks)
        importlib.reload(test_suite)
        loader = unittest.TestLoader()
        suite = loader.loadTestsFromTestCase(test_suite.TestOpticsMath)
        
        for num, _ in self.tasks:
            level_suite = unittest.TestSuite([t for t in suite if f"test_l{num}" in t._testMethodName])
            res = unittest.TextTestRunner(verbosity=0).run(level_suite)
            canvas, light = self.status_indicators[num]
            color = "#4ec9b0" if (res.wasSuccessful() and res.testsRun > 0) else "#f44747" if res.testsRun > 0 else "red"
            canvas.itemconfig(light, fill=color)
        self.update_plot()

    def switch_sandbox(self, level):
        self.current_lvl = level
        self.update_plot()

    def update_plot(self):
        self.ax.clear()
        self.ax.set_facecolor("#1e1e1e")
        self.ax.set_xlim(-60, 60); self.ax.set_ylim(-50, 50); self.ax.set_aspect('equal')
        self.ax.grid(True, color="#333333", linestyle='--')
        
        r_o, r_t = self.handles['ray_o'], self.handles['ray_t']
        center = self.handles['center']
        ray_dir = (r_t - r_o) / (np.linalg.norm(r_t - r_o) + 1e-9)
        
        GUIDE_COLOR, SURF_COLOR, RAY_COLOR, HANDLE_COLOR = "#ff0000", "#00ffff", "#ff8c00", "#ffff00"
        t_raw, curve = float('inf'), None

        # --- LINEAR LEVELS (1, 2, 5, 7) ---
        if self.current_lvl in [1, 2, 5, 7]:
            p1, p2 = self.handles['p1'], self.handles['p2']
            slope = p2 - p1
            self.ax.plot([p1[0]-100*slope[0], p1[0]+100*slope[0]], [p1[1]-100*slope[1], p1[1]+100*slope[1]], 
                         color=GUIDE_COLOR, linestyle="--", alpha=0.5, lw=1)
            self.ax.plot([p1[0], p2[0]], [p1[1], p2[1]], color=SURF_COLOR, lw=3, zorder=3)
            
            if self.current_lvl == 1:
                t_raw = tasks.intersect_line_infinite(r_o, ray_dir, p1, p2)
            else:
                t_raw = tasks.intersect_segment(r_o, ray_dir, p1, p2)
            curve = {'type': 'line', 'p1': p1, 'p2': p2}

        # --- CIRCLE LEVEL (3) ---
        elif self.current_lvl == 3:
            radius = np.linalg.norm(self.handles['radius_h'] - center)
            self.ax.add_patch(Circle(center, radius, color=SURF_COLOR, fill=False, lw=3, zorder=3))
            t_raw = tasks.intersect_circle_infinite(r_o, ray_dir, center, radius)
            curve = {'type': 'circle', 'center': center, 'radius': radius}

        # --- ARC LEVELS (4, 6) ---
        elif self.current_lvl in [4, 6]:
            radius = np.linalg.norm(self.handles['arc_p1'] - center)
            self.ax.add_patch(Circle(center, radius, color=GUIDE_COLOR, fill=False, lw=1, linestyle="--", alpha=0.3))
            
            # Recalculate Arc Geometry from Handles
            mid = (self.handles['arc_p1'] + self.handles['arc_p2']) / 2.0
            axis = (mid - center) / (np.linalg.norm(mid - center) + 1e-9)
            vec_p1 = (self.handles['arc_p1'] - center) / (radius + 1e-9)
            cos_half = np.dot(vec_p1, axis)

            span = np.degrees(np.arccos(np.clip(cos_half, -1, 1))) * 2
            self.ax.add_patch(Arc(center, radius*2, radius*2, angle=np.degrees(np.arctan2(axis[1], axis[0])), 
                                  theta1=-span/2, theta2=span/2, color=SURF_COLOR, lw=4, zorder=3))
            
            t_raw = tasks.intersect_arc(r_o, ray_dir, center, radius, axis, cos_half)
            curve = {'type': 'arc', 'center': center, 'radius': radius, 'axis': axis, 'cos_half_angle': cos_half}

        # Handle Return Types (L1/L3 return lists or tuples)
        t = t_raw[0] if isinstance(t_raw, (list, tuple)) and len(t_raw) > 0 else t_raw
        if isinstance(t, (list, tuple)): t = float('inf') # Safety for double nesting
        
        is_hit = (t is not None) and (t != float('inf')) and (t > 1e-4)
        hit_dist = t if is_hit else 150
        hit_pt = r_o + hit_dist * ray_dir
        self.ax.plot([r_o[0], hit_pt[0]], [r_o[1], hit_pt[1]], color=RAY_COLOR, lw=2)

        if is_hit:
            self.ax.scatter(hit_pt[0], hit_pt[1], c="#4ec9b0", s=120, edgecolors='white', zorder=10)
            if self.current_lvl in [5, 6, 7]:
                norm = tasks.calculate_normal(hit_pt, ray_dir, curve)
                if norm is not None and np.linalg.norm(norm) > 1e-3:
                    self.ax.quiver(hit_pt[0], hit_pt[1], norm[0], norm[1], color="#c586c0", scale=12, pivot='tail')
                    if self.current_lvl == 7:
                        refr_dir = tasks.refract_vector(ray_dir, norm, 1.0, 1.5)
                        if refr_dir is not None:
                            refr_end = hit_pt + refr_dir * 40
                            self.ax.plot([hit_pt[0], refr_end[0]], [hit_pt[1], refr_end[1]], color="#00ff00", lw=2)

        # Selective Handle Visibility
        vis = ['ray_o', 'ray_t']
        if self.current_lvl in [1, 2, 5, 7]: vis += ['p1', 'p2']
        if self.current_lvl == 3: vis += ['center', 'radius_h']
        if self.current_lvl in [4, 6]: vis += ['center', 'arc_p1', 'arc_p2']
        
        for k in vis:
            v = self.handles[k]
            self.ax.scatter(v[0], v[1], c=HANDLE_COLOR, s=60, edgecolors='black', zorder=5)

        self.canvas.draw()

    def on_press(self, event):
        if event.xdata is None: return
        # Filter for visible handles only
        vis = ['ray_o', 'ray_t']
        if self.current_lvl in [1, 2, 5, 7]: vis += ['p1', 'p2']
        if self.current_lvl == 3: vis += ['center', 'radius_h']
        if self.current_lvl in [4, 6]: vis += ['center', 'arc_p1', 'arc_p2']
        
        for k in vis:
            v = self.handles[k]
            if np.hypot(event.xdata - v[0], event.ydata - v[1]) < 4: self.dragging = k; break

    def on_motion(self, event):
        if self.dragging and event.xdata is not None:
            new = np.array([event.xdata, event.ydata])
            if self.dragging == 'center':
                # Sync all geometric handles to the center
                d_rh = self.handles['radius_h'] - self.handles['center']
                d1, d2 = self.handles['arc_p1'] - self.handles['center'], self.handles['arc_p2'] - self.handles['center']
                self.handles['center'] = new
                self.handles['radius_h'], self.handles['arc_p1'], self.handles['arc_p2'] = new + d_rh, new + d1, new + d2
            elif self.dragging == 'arc_p1':
                self.handles['arc_p1'] = new
                r = np.linalg.norm(new - self.handles['center'])
                v2 = self.handles['arc_p2'] - self.handles['center']
                self.handles['arc_p2'] = self.handles['center'] + r * (v2 / (np.linalg.norm(v2) + 1e-9))
            elif self.dragging == 'arc_p2':
                r = np.linalg.norm(self.handles['arc_p1'] - self.handles['center'])
                v2 = new - self.handles['center']
                self.handles['arc_p2'] = self.handles['center'] + r * (v2 / (np.linalg.norm(v2) + 1e-9))
            else:
                self.handles[self.dragging] = new
            self.update_plot()

    def on_release(self, event): self.dragging = None

if __name__ == "__main__":
    app = OpticsDebugger()
    app.mainloop()

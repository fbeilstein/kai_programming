import tkinter as tk
from tkinter import ttk, messagebox
import unittest
import numpy as np
import subprocess
import sys
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

        # --- ACTION BUTTONS ---
        tk.Button(self.sidebar, text="🔄 Reload & Retest", bg="#3e3e42", fg="white",
                  command=self.refresh_tests).pack(pady=(30, 10), padx=20, fill="x")
        
        tk.Button(self.sidebar, text="🚀 Run Main Simulation", bg="#007acc", fg="white", font=("Arial", 10, "bold"),
                  command=self.run_main_script).pack(pady=10, padx=20, fill="x")

    def run_main_script(self):
        try:
            subprocess.Popen([sys.executable, "caustics.py"])
        except Exception as e:
            messagebox.showerror("Execution Error", f"Could not launch caustics.py: {e}")

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
            'arc_p1': np.array([0.0, 15.0]), # Controls Radius + Start Angle
            'arc_p2': np.array([0.0, -15.0]), # Controls End Angle
            'radius_h': np.array([30.0, 0.0]) # Circle-only handle
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
        self.ax.clear()
        self.ax.set_facecolor("#1e1e1e")
        self.ax.set_xlim(-60, 60); self.ax.set_ylim(-50, 50); self.ax.set_aspect('equal')
        self.ax.grid(True, color="#333333", linestyle='--')
        
        r_o, r_t = self.handles['ray_o'], self.handles['ray_t']
        ray_dir = (r_t - r_o) / np.linalg.norm(r_t - r_o)
        
        # BRIGHT VIBRANT COLORS
        GUIDE_COLOR, SURF_COLOR, RAY_COLOR, HANDLE_COLOR = "#ff0000", "#00ffff", "#ff8c00", "#ffff00"
        t = float('inf')
        curve = None

        if self.current_lvl in [1, 2, 5, 7]:
            p1, p2 = self.handles['p1'], self.handles['p2']
            slope = p2 - p1
            self.ax.plot([p1[0]-100*slope[0], p1[0]+100*slope[0]], [p1[1]-100*slope[1], p1[1]+100*slope[1]], 
                         color=GUIDE_COLOR, linestyle="--", alpha=0.4, lw=1)
            self.ax.plot([p1[0], p2[0]], [p1[1], p2[1]], color=SURF_COLOR, lw=3, zorder=3)
            t = tasks.intersect_line_infinite(r_o, ray_dir, p1, p2)[0] if self.current_lvl == 1 else tasks.intersect_segment(r_o, ray_dir, p1, p2)
            curve = {'type': 'line', 'p1': p1, 'p2': p2}
            
        elif self.current_lvl in [3, 4, 6]:
            center = self.handles['center']
            radius = np.linalg.norm(self.handles['arc_p1'] - center) if self.current_lvl != 3 else np.linalg.norm(self.handles['radius_h'] - center)
            
            # Show full ghost circle in bright red
            self.ax.add_patch(Circle(center, radius, color=GUIDE_COLOR, fill=False, lw=1, linestyle="--", alpha=0.3))

            if self.current_lvl == 3:
                ts = tasks.intersect_circle_infinite(r_o, ray_dir, center, radius)
                t = min([val for val in ts if val > 1e-4], default=float('inf'))
                self.ax.add_patch(Circle(center, radius, color=SURF_COLOR, fill=False, lw=2, alpha=0.3))
            else:
                # Arc construction
                midpoint = (self.handles['arc_p1'] + self.handles['arc_p2']) / 2.0
                axis = (midpoint - center); axis /= np.linalg.norm(axis)
                vec_p1 = (self.handles['arc_p1'] - center) / radius
                cos_half = np.dot(vec_p1, axis)
                
                t = tasks.intersect_arc(r_o, ray_dir, center, radius, axis, cos_half)
                angle_axis = np.degrees(np.arctan2(axis[1], axis[0]))
                span = np.degrees(np.arccos(np.clip(cos_half, -1, 1))) * 2
                self.ax.add_patch(Arc(center, radius*2, radius*2, angle=angle_axis, 
                                      theta1=-span/2, theta2=span/2, color=SURF_COLOR, lw=3))
            curve = {'type': 'arc', 'center': center, 'radius': radius, 'axis': axis if self.current_lvl != 3 else None, 'cos_half_angle': cos_half if self.current_lvl != 3 else None}

        # Draw Ray
        hit_dist = t if (t != float('inf') and t > 0) else 150
        hit_pt = r_o + hit_dist * ray_dir
        self.ax.plot([r_o[0], hit_pt[0]], [r_o[1], hit_pt[1]], color=RAY_COLOR, lw=2, zorder=2)
        
        if self.current_lvl == 7 and t != float('inf') and t > 0:
            norm = tasks.calculate_normal(hit_pt, ray_dir, curve)
            refr_dir = tasks.refract_vector(ray_dir, norm, 1.0, 1.5)
            if refr_dir is not None:
                refr_end = hit_pt + refr_dir * 40
                self.ax.plot([hit_pt[0], refr_end[0]], [hit_pt[1], refr_end[1]], color="#00ff00", lw=2)

        # Draw Draggable Handles
        visible = ['ray_o', 'ray_t']
        if self.current_lvl in [1, 2, 5, 7]: visible += ['p1', 'p2']
        if self.current_lvl == 3: visible += ['center', 'radius_h']
        if self.current_lvl in [4, 6]: visible += ['center', 'arc_p1', 'arc_p2']
        
        for k in visible:
            v = self.handles[k]
            self.ax.scatter(v[0], v[1], c=HANDLE_COLOR, s=60, edgecolors='black', zorder=5)

        if t != float('inf') and t > 0:
            self.ax.scatter(hit_pt[0], hit_pt[1], c="#4ec9b0", s=120, edgecolors='white', zorder=10)
            if self.current_lvl in [5, 6, 7]:
                norm = tasks.calculate_normal(hit_pt, ray_dir, curve)
                self.ax.quiver(hit_pt[0], hit_pt[1], norm[0], norm[1], color="#c586c0", scale=12, pivot='tail')

        self.canvas.draw()

    def on_motion(self, event):
        if self.dragging and event.xdata is not None:
            new_pos = np.array([event.xdata, event.ydata])
            if self.dragging == 'center':
                diff1, diff2 = self.handles['arc_p1'] - self.handles['center'], self.handles['arc_p2'] - self.handles['center']
                diff_rh = self.handles['radius_h'] - self.handles['center']
                self.handles['center'] = new_pos
                self.handles['arc_p1'], self.handles['arc_p2'] = new_pos + diff1, new_pos + diff2
                self.handles['radius_h'] = new_pos + diff_rh
            elif self.dragging == 'arc_p1':
                # Master handle sets radius and angle
                self.handles['arc_p1'] = new_pos
                # Snap arc_p2 to the new radius
                r = np.linalg.norm(new_pos - self.handles['center'])
                vec2 = self.handles['arc_p2'] - self.handles['center']
                self.handles['arc_p2'] = self.handles['center'] + r * (vec2 / np.linalg.norm(vec2))
            elif self.dragging == 'arc_p2':
                # Sibling handle only changes angle
                r = np.linalg.norm(self.handles['arc_p1'] - self.handles['center'])
                vec2 = new_pos - self.handles['center']
                self.handles['arc_p2'] = self.handles['center'] + r * (vec2 / np.linalg.norm(vec2))
            else:
                self.handles[self.dragging] = new_pos
            self.update_plot()

    def on_press(self, event):
        if event.xdata is None: return
        for k, v in self.handles.items():
            if np.hypot(event.xdata - v[0], event.ydata - v[1]) < 4: self.dragging = k; break
    def on_release(self, event): self.dragging = None
    def switch_sandbox(self, level): self.current_lvl = level; self.update_plot()

if __name__ == "__main__":
    app = OpticsDebugger()
    app.mainloop()

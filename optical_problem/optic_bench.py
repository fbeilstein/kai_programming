import tkinter as tk
from tkinter import ttk
from matplotlib.figure import Figure 
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.patches import Polygon, Rectangle
from matplotlib.path import Path
import numpy as np
from implementation_tasks import trace_ray_step

from lens_object import LensObject
from lens_architect import LensArchitect


class OpticBenchApp(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Physics Lab: Optic Bench (Simulation Ready)")
        self.geometry("1100x650")
        self.lenses = []; self.selected_lens = None; self.dragging_lens = None; self.dragging_source = False
        self.last_mouse_pos = (0, 0); self.start_mouse_angle = 0; self.start_lens_angle = 0
        self.source_type = "parallel"; self.source_pos = [-100, 20]

        toolbar = ttk.Frame(self, padding=5); toolbar.pack(side=tk.TOP, fill=tk.X)
        ttk.Button(toolbar, text="+ Add Lens", command=self.open_architect).pack(side=tk.LEFT, padx=5)
        ttk.Button(toolbar, text="Delete", command=self.delete_lens).pack(side=tk.LEFT, padx=5)
        self.lbl_status = ttk.Label(toolbar, text="Ready", foreground="green")
        self.lbl_status.pack(side=tk.RIGHT, padx=10)
        ttk.Separator(toolbar, orient=tk.VERTICAL).pack(side=tk.LEFT, padx=10, fill=tk.Y)
        self.btn_source = ttk.Button(toolbar, text="Mode: Parallel", command=self.toggle_source); self.btn_source.pack(side=tk.LEFT, padx=5)
        ttk.Label(toolbar, text="Lens n:").pack(side=tk.LEFT)
        self.n_var = tk.DoubleVar(value=1.52)
        self.spin_n = ttk.Spinbox(toolbar, from_=1.0, to=3.0, increment=0.01, textvariable=self.n_var, width=5, command=self.update_n)
        self.spin_n.pack(side=tk.LEFT); self.spin_n.bind('<Return>', lambda e: self.update_n())

        self.fig = Figure(figsize=(5, 4), dpi=100); self.ax = self.fig.add_subplot(111)
        self.canvas = FigureCanvasTkAgg(self.fig, master=self); self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)
        self.setup_plot(); self.connect_events()

    def setup_plot(self):
        self.ax.set_xlim(-150, 150); self.ax.set_ylim(-80, 80); self.ax.set_aspect('equal')
        self.ax.grid(True, alpha=0.3); self.draw_scene()

    def connect_events(self):
        self.canvas.mpl_connect('button_press_event', self.on_press)
        self.canvas.mpl_connect('motion_notify_event', self.on_drag)
        self.canvas.mpl_connect('button_release_event', self.on_release)

    def open_architect(self): LensArchitect(self, self.add_lens_callback)
    def add_lens_callback(self, geometry_data):
        new_lens = LensObject(geometry_data, x_pos=0, y_pos=0)
        spawn_found = False
        for offset in [0, 20, -20, 40, -40, 60, -60, 80, -80]:
            new_lens.x = offset
            collision = False
            for other in self.lenses:
                if new_lens.intersects(other): collision = True; break
            if not collision: spawn_found = True; break
        if spawn_found: self.lenses.append(new_lens); self.draw_scene()
        else: self.lbl_status.config(text="Error: Bench too crowded!", foreground="red")

    def delete_lens(self):
        if self.selected_lens and self.selected_lens in self.lenses: self.lenses.remove(self.selected_lens); self.selected_lens = None; self.draw_scene()
    def update_n(self):
        if self.selected_lens: self.selected_lens.n = self.n_var.get(); self.draw_scene()

    def on_press(self, event):
        if event.xdata is None: return
        self.last_mouse_pos = (event.xdata, event.ydata)
        if self.source_type == "point":
            if np.hypot(event.xdata - self.source_pos[0], event.ydata - self.source_pos[1]) < 8: self.dragging_source = True; return
        self.selected_lens = None
        for lens in reversed(self.lenses):
            if lens.contains((event.xdata, event.ydata)):
                self.selected_lens = lens; self.dragging_lens = lens; self.n_var.set(lens.n)
                dx, dy = event.xdata - lens.x, event.ydata - lens.y
                self.start_mouse_angle = np.degrees(np.arctan2(dy, dx)); self.start_lens_angle = lens.angle; break
        self.draw_scene()

    def on_drag(self, event):
        if event.xdata is None: return
        mouse_dx = event.xdata - self.last_mouse_pos[0]; mouse_dy = event.ydata - self.last_mouse_pos[1]
        self.last_mouse_pos = (event.xdata, event.ydata)

        if self.dragging_source: self.source_pos = [event.xdata, event.ydata]; self.draw_scene()
        elif self.dragging_lens:
            if np.hypot(event.xdata - self.dragging_lens.x, event.ydata - self.dragging_lens.y) > 60:
                self.dragging_lens = None; self.lbl_status.config(text="Lost grip!", foreground="orange"); self.draw_scene(); return
            old_x, old_y, old_angle = self.dragging_lens.x, self.dragging_lens.y, self.dragging_lens.angle
            if event.button == 1: self.dragging_lens.x += mouse_dx; self.dragging_lens.y += mouse_dy
            elif event.button == 3:
                dx, dy = event.xdata - self.dragging_lens.x, event.ydata - self.dragging_lens.y
                if np.hypot(dx, dy) > 1:
                    delta = np.degrees(np.arctan2(dy, dx)) - self.start_mouse_angle
                    self.dragging_lens.angle = self.start_lens_angle + delta
            collision = False
            for other in self.lenses:
                if other is not self.dragging_lens and self.dragging_lens.intersects(other): collision = True; break
            if collision:
                self.dragging_lens.x = old_x; self.dragging_lens.y = old_y; self.dragging_lens.angle = old_angle
                self.lbl_status.config(text="Blocked!", foreground="red")
            else: self.lbl_status.config(text="Ready", foreground="green")
            self.draw_scene()

    def on_release(self, event): self.dragging_lens = None; self.dragging_source = False; self.draw_scene()
    def toggle_source(self):
        self.source_type = "point" if self.source_type == "parallel" else "parallel"
        self.btn_source.config(text=f"Mode: {self.source_type} Source"); self.draw_scene()

    # --- THE SIMULATION LOOP ---
    def run_simulation(self):
        # 1. Build "Media" and "Curves" lists for the Student Function
        media = {0: {'name': 'Air', 'n': 1.0}}
        all_curves = []
        
        for i, lens in enumerate(self.lenses):
            medium_id = i + 1
            media[medium_id] = {'name': f'Lens_{i}', 'n': lens.n}
            all_curves.extend(lens.get_curves(medium_id))
            
        # 2. Determine Source Rays
        rays = []
        if self.source_type == "parallel":
            # FULL HEIGHT COVERAGE (-75 to 75)
            # Increased density to 50 rays
            for y in np.linspace(-75, 75, 50):
                rays.append({'origin': np.array([-150.0, y]), 'dir': np.array([1.0, 0.0]), 'medium': 0})
        else:
            # 360 DEGREE POINT SOURCE
            src = np.array([self.source_pos[0], self.source_pos[1]])
            
            # Check if source is starting inside a lens
            start_medium = 0
            for i, lens in enumerate(self.lenses):
                if lens.contains(src): start_medium = i + 1; break
            
            # Full circle (0 to 2pi), 120 rays for smooth look
            # We omit the last point to avoid duplicating 0 and 360
            for angle in np.linspace(0, 2*np.pi, 120, endpoint=False):
                d = np.array([np.cos(angle), np.sin(angle)])
                rays.append({'origin': src, 'dir': d, 'medium': start_medium})

        # 3. Trace Rays
        for ray in rays:
            curr_pos = ray['origin']
            curr_dir = ray['dir']
            curr_med = ray['medium']
            
            path_points = [curr_pos]
            
            for step in range(12): # Increased bounces to 12 for complex traps
                hit, new_dir, new_med = trace_ray_step(curr_pos, curr_dir, curr_med, all_curves, media)
                
                if hit is None:
                    # Ray flies to infinity
                    # Extend visualization slightly off-screen
                    path_points.append(curr_pos + curr_dir * 300)
                    break
                else:
                    path_points.append(hit)
                    if new_dir is None: break # Absorbed/TIR
                    
                    # Nudge forward to prevent self-intersection on next step
                    curr_pos = hit + new_dir * 1e-3 
                    curr_dir = new_dir
                    curr_med = new_med
            
            # Draw Path
            pts = np.array(path_points)
            self.ax.plot(pts[:,0], pts[:,1], 'r-', linewidth=0.8, alpha=0.5)

    def draw_scene(self):
        for artist in list(self.ax.patches) + list(self.ax.lines) + list(self.ax.texts): artist.remove()
        self.ax.axhline(0, color='black', linewidth=0.5, alpha=0.5)

        # Draw Lenses
        for lens in self.lenses:
            pts = lens.get_transformed_points()
            color = 'yellow' if lens == self.selected_lens else 'cyan'
            poly = Polygon(pts, closed=True, facecolor=color, edgecolor='blue', alpha=0.3)
            self.ax.add_patch(poly)
            self.ax.text(lens.x, lens.y+40, f"n={lens.n}", ha='center', fontsize=8, color='blue', zorder=20)
            
        # Draw Source Marker
        if self.source_type == 'point':
            self.ax.plot(self.source_pos[0], self.source_pos[1], 'go', markersize=6, zorder=30)

        # RUN SIMULATION
        self.run_simulation()
        
        self.canvas.draw()

if __name__ == "__main__":
    app = OpticBenchApp()
    app.mainloop()

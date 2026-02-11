import tkinter as tk
from tkinter import ttk
from matplotlib.figure import Figure 
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.patches import Polygon, Rectangle
from matplotlib.path import Path
import numpy as np
from implementation_tasks import trace_ray_step

from lens_object import LensObject


class LensArchitect(tk.Toplevel):
    def __init__(self, parent, on_save_callback):
        super().__init__(parent)
        self.title("Lens Designer")
        self.geometry("800x600")
        self.on_save = on_save_callback
        self.transient(parent); self.grab_set()

        self.width = 15.0; self.diameter = 40.0
        self.sag_front = 10.0; self.sag_back = 10.0
        self.geometry_cache = {}; self.dragging = None

        self.fig = Figure(figsize=(6, 5), dpi=100)
        self.ax = self.fig.add_subplot(111)
        self.canvas = FigureCanvasTkAgg(self.fig, master=self)
        self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        btn_frame = ttk.Frame(self)
        btn_frame.pack(side=tk.BOTTOM, fill=tk.X, pady=10)
        self.lbl_warn = ttk.Label(btn_frame, text="Drag handles to shape lens.", foreground="black")
        self.lbl_warn.pack(side=tk.LEFT, padx=10)
        ttk.Button(btn_frame, text="Add to Bench", command=self.finish).pack(side=tk.RIGHT, padx=10)

        self.ax.set_xlim(-60, 60); self.ax.set_ylim(-40, 40); self.ax.set_aspect('equal')
        self.ax.grid(True, alpha=0.3)
        self.ax.set_title("Top View (Drag Handles)", fontsize=10)
        
        self.ref_box = Rectangle((-7.5, -20), 15, 40, edgecolor='gray', facecolor='none', linestyle='--', alpha=0.5)
        self.ax.add_patch(self.ref_box)
        self.lens_poly = Polygon([[0,0]], closed=True, facecolor='cyan', edgecolor='blue', alpha=0.4)
        self.ax.add_patch(self.lens_poly)
        
        self.h_front, = self.ax.plot([], [], 'ro', markersize=8, picker=10, markeredgecolor='k', label="Front")
        self.h_back, = self.ax.plot([], [], 'bo', markersize=8, picker=10, markeredgecolor='k', label="Back")
        self.h_height, = self.ax.plot([], [], 'gs', markersize=8, picker=10, markeredgecolor='k', label="Diam")

        self.canvas.mpl_connect('pick_event', self.on_pick)
        self.canvas.mpl_connect('motion_notify_event', self.on_drag)
        self.canvas.mpl_connect('button_release_event', self.on_release)
        self.update_visuals()

    def calculate_arc(self, sagitta, aperture, is_left_side):
        h = aperture / 2.0; s = sagitta
        if abs(s) < 0.1:
            x_pos = -self.width/2 if is_left_side else self.width/2
            return {'type': 'line', 'x': x_pos, 'y_top': h, 'y_bot': -h}
        R_signed = (s**2 + h**2) / (2 * s); R = abs(R_signed)
        base_x = -self.width/2 if is_left_side else self.width/2
        direction = -1 if is_left_side else 1 
        apex_x = base_x + (direction * s); center_x = apex_x - (direction * R_signed)
        if h > R: h = R 
        alpha = np.degrees(np.arcsin(h / R))
        angle_to_apex = 0 if (apex_x > center_x) else 180
        t1, t2 = angle_to_apex - alpha, angle_to_apex + alpha
        return {'type': 'arc', 'center': (center_x, 0), 'radius': R, 
                'theta1': min(t1, t2), 'theta2': max(t1, t2), 'apex_x': apex_x}

    def get_poly_path(self, front, back):
        def get_pts(geo, n=30):
            if geo['type'] == 'line': return np.array([[geo['x'], geo['y_top']], [geo['x'], geo['y_bot']]])
            thetas = np.radians(np.linspace(geo['theta1'], geo['theta2'], n))
            xs = geo['center'][0] + geo['radius'] * np.cos(thetas)
            ys = geo['center'][1] + geo['radius'] * np.sin(thetas)
            pts = np.column_stack([xs, ys])
            return pts[pts[:,1].argsort()[::-1]] 
        pts_f = get_pts(front); pts_b = get_pts(back)
        return Path(pts_f), Path(pts_b)

    def update_visuals(self):
        self.ref_box.set_width(self.width); self.ref_box.set_height(self.diameter); self.ref_box.set_xy((-self.width/2, -self.diameter/2))
        front = self.calculate_arc(self.sag_front, self.diameter, True)
        back = self.calculate_arc(self.sag_back, self.diameter, False)
        self.geometry_cache = {'front': front, 'back': back}

        fx = front['x'] if front['type'] == 'line' else front['apex_x']
        bx = back['x'] if back['type'] == 'line' else back['apex_x']
        self.h_front.set_data([fx], [0]); self.h_back.set_data([bx], [0]); self.h_height.set_data([0], [self.diameter/2])

        path_f, path_b = self.get_poly_path(front, back)
        poly_pts = np.vstack([path_f.vertices, path_b.vertices[::-1]])
        self.lens_poly.set_xy(poly_pts)
        self.canvas.draw()

    def check_self_intersection(self, test_sag_front, test_sag_back, test_diam):
        f = self.calculate_arc(test_sag_front, test_diam, True)
        b = self.calculate_arc(test_sag_back, test_diam, False)
        path_f, path_b = self.get_poly_path(f, b)
        if path_f.intersects_path(path_b, filled=False): return True
        fx = f['x'] if f['type']=='line' else f['apex_x']
        bx = b['x'] if b['type']=='line' else b['apex_x']
        if fx >= bx - 1.0: return True 
        return False

    def on_pick(self, event): self.dragging = event.artist
    def on_release(self, event): self.dragging = None
    
    def on_drag(self, event):
        if not self.dragging or event.xdata is None: return
        x, y = event.xdata, event.ydata
        new_sag_f = self.sag_front; new_sag_b = self.sag_back; new_diam  = self.diameter
        
        if self.dragging == self.h_front: new_sag_f = (-self.width / 2) - x
        elif self.dragging == self.h_back: new_sag_b = x - (self.width / 2)
        elif self.dragging == self.h_height: new_diam = abs(y) * 2

        if self.check_self_intersection(new_sag_f, new_sag_b, new_diam):
            self.lbl_warn.config(text="Limit Reached!", foreground="red")
        else:
            self.lbl_warn.config(text="Drag handles to shape lens.", foreground="black")
            self.sag_front = new_sag_f; self.sag_back = new_sag_b; self.diameter = new_diam
            self.update_visuals()

    def finish(self):
        if self.on_save: self.on_save([self.geometry_cache['front'], self.geometry_cache['back']])
        self.destroy()


import tkinter as tk
from tkinter import ttk
from matplotlib.figure import Figure 
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.patches import Polygon, Rectangle
from matplotlib.path import Path
import numpy as np

# =============================================================================
#  SECTION 1: THE "STUDENT CODE" (EXEMPLARY IMPLEMENTATION)
#  This is the logic students are expected to write.
# =============================================================================

def trace_ray_step(ray_origin, ray_dir, current_medium_id, curves, media):
    """
    Advances the ray to the next intersection and refracts it.
    """
    
    # --- 1. GEOMETRY: FIND NEAREST INTERSECTION ---
    min_t = float('inf')
    best_hit = None
    hit_curve = None
    
    for curve in curves:
        t = float('inf')
        
        # A. Line Intersection
        if curve['type'] == 'line':
            p1, p2 = curve['p1'], curve['p2']
            v = p2 - p1
            denom = np.cross(v, ray_dir)
            if abs(denom) > 1e-9:
                u = np.cross(ray_origin - p1, ray_dir) / denom
                t_candidate = np.cross(ray_origin - p1, v) / denom
                if 0 <= u <= 1 and t_candidate > 1e-4: # Epsilon for self-intersection
                    t = t_candidate

        # B. Arc Intersection
        elif curve['type'] == 'arc':
            C = curve['center']
            R = curve['radius']
            OC = ray_origin - C
            
            # Quadratic: |O + td - C|^2 = R^2
            # t^2 + 2(d . OC)t + (|OC|^2 - R^2) = 0
            b = 2 * np.dot(ray_dir, OC)
            c = np.dot(OC, OC) - R**2
            discriminant = b**2 - 4*c
            
            if discriminant >= 0:
                sqrt_d = np.sqrt(discriminant)
                t1 = (-b - sqrt_d) / 2
                t2 = (-b + sqrt_d) / 2
                
                # Check both solutions
                for ti in [t1, t2]:
                    if ti > 1e-4:
                        P = ray_origin + ti * ray_dir
                        # Cone Check: Is P within the angular sector?
                        vec_CP = P - C
                        # Dot product > threshold
                        if np.dot(vec_CP, curve['axis']) >= R * curve['cos_half_angle'] - 1e-4:
                            if ti < t: t = ti

        if t < min_t:
            min_t = t
            best_hit = ray_origin + min_t * ray_dir
            hit_curve = curve

    if best_hit is None:
        return None, None, None # Ray escapes to infinity

    # --- 2. PHYSICS: DETERMINE INDICES ---
    n1 = media[current_medium_id]['n']
    
    if hit_curve['medium_id'] == current_medium_id:
        # Exiting current medium -> Air
        new_medium_id = 0
    else:
        # Entering new medium
        new_medium_id = hit_curve['medium_id']
        
    n2 = media[new_medium_id]['n']

    # --- 3. REFRACTION: SNELL'S LAW ---
    # Normal Vector Calculation
    if hit_curve['type'] == 'line':
        # Normal is perpendicular to segment
        p1, p2 = hit_curve['p1'], hit_curve['p2']
        tangent = p2 - p1
        normal = np.array([-tangent[1], tangent[0]])
        normal /= np.linalg.norm(normal)
    else:
        # Normal is Radius vector (Center -> Hit)
        normal = (best_hit - hit_curve['center'])
        normal /= np.linalg.norm(normal)

    # Ensure normal points AGAINST incident ray for standard formula
    if np.dot(ray_dir, normal) > 0:
        normal = -normal

    # Vector Snell's Law
    eta = n1 / n2
    cos_theta1 = -np.dot(ray_dir, normal)
    sin2_theta1 = 1 - cos_theta1**2
    sin2_theta2 = eta**2 * sin2_theta1
    
    if sin2_theta2 > 1.0:
        # Total Internal Reflection (TIR) - Optional
        # For this demo, let's just stop the ray or reflect
        return best_hit, None, None 
    else:
        cos_theta2 = np.sqrt(1 - sin2_theta2)
        new_dir = eta * ray_dir + (eta * cos_theta1 - cos_theta2) * normal
        new_dir /= np.linalg.norm(new_dir) # Re-normalize to be safe

    return best_hit, new_dir, new_medium_id

# =============================================================================
#  SECTION 2: LENS ARCHITECT (Unchanged)
# =============================================================================
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


# =============================================================================
#  SECTION 3: OPTIC BENCH (Now with Real-Time Simulation)
# =============================================================================
class LensObject:
    def __init__(self, geometry_data, x_pos=0, y_pos=0, angle=0, n=1.52):
        self.geo = geometry_data
        self.x = x_pos; self.y = y_pos; self.angle = angle; self.n = n
        self.selected = False; self.collision_state = False

    def get_transformed_points(self):
        # Used for Drawing and Collision (Visuals)
        front, back = self.geo[0], self.geo[1]
        def get_arc_pts(g, n=20):
            if g['type'] == 'line': return np.array([[g['x'], g['y_top']], [g['x'], g['y_bot']]])
            thetas = np.radians(np.linspace(g['theta1'], g['theta2'], n))
            cx, cy = g['center']
            xs = cx + g['radius'] * np.cos(thetas)
            ys = cy + g['radius'] * np.sin(thetas)
            return np.column_stack([xs, ys])

        pts_f = get_arc_pts(front); pts_f = pts_f[pts_f[:,1].argsort()[::-1]] 
        pts_b = get_arc_pts(back); pts_b = pts_b[pts_b[:,1].argsort()] 
        pts = np.vstack([pts_f, pts_b])
        theta = np.radians(self.angle)
        c, s = np.cos(theta), np.sin(theta)
        R = np.array([[c, -s], [s, c]])
        rotated_pts = pts @ R.T
        return rotated_pts + [self.x, self.y]

    def contains(self, point):
        return Path(self.get_transformed_points()).contains_point(point)

    def intersects(self, other_lens):
        pts1 = self.get_transformed_points(); pts2 = other_lens.get_transformed_points()
        if (np.max(pts1[:,0]) < np.min(pts2[:,0]) or np.min(pts1[:,0]) > np.max(pts2[:,0]) or
            np.max(pts1[:,1]) < np.min(pts2[:,1]) or np.min(pts1[:,1]) > np.max(pts2[:,1])): return False
        path1 = Path(pts1); path2 = Path(pts2)
        if path1.intersects_path(path2, filled=True): return True
        if any(path2.contains_points(pts1)): return True
        if any(path1.contains_points(pts2)): return True
        if path2.contains_point(np.mean(pts1, axis=0)): return True
        if path1.contains_point(np.mean(pts2, axis=0)): return True
        return False

    def get_curves(self, medium_id):
        """ EXPORT FOR STUDENT SIMULATOR """
        curves = []
        
        # 1. Transform Helpers
        theta = np.radians(self.angle)
        c, s = np.cos(theta), np.sin(theta)
        R_mat = np.array([[c, -s], [s, c]])
        
        def transform(p): return (np.array(p) @ R_mat.T) + [self.x, self.y]
        def rotate_vec(v): return np.array(v) @ R_mat.T

        # 2. Helper to find endpoints of a geometry part (Top/Bot)
        def get_local_endpoints(part):
            if part['type'] == 'line':
                # Line endpoints are explicit
                return [part['x'], part['y_top']], [part['x'], part['y_bot']]
            elif part['type'] == 'arc':
                # Arc endpoints are at theta1/theta2
                # We calculate them and sort by Y to find Top/Bot
                c = np.array(part['center'])
                r = part['radius']
                p1 = c + r * np.array([np.cos(np.radians(part['theta1'])), np.sin(np.radians(part['theta1']))])
                p2 = c + r * np.array([np.cos(np.radians(part['theta2'])), np.sin(np.radians(part['theta2']))])
                # Return Top (Higher Y), Bot (Lower Y)
                if p1[1] > p2[1]: return p1, p2
                return p2, p1

        # 3. Process Main Surfaces (Front/Back)
        corners = [] # Will store [(F_top, F_bot), (B_top, B_bot)]
        
        for part in self.geo:
            # Save endpoints for the side walls
            top_local, bot_local = get_local_endpoints(part)
            corners.append((transform(top_local), transform(bot_local)))

            if part['type'] == 'line':
                p1 = transform([part['x'], part['y_top']])
                p2 = transform([part['x'], part['y_bot']])
                curves.append({'type': 'line', 'p1': p1, 'p2': p2, 'medium_id': medium_id})
            
            elif part['type'] == 'arc':
                center = transform(part['center'])
                local_apex = np.array([part['apex_x'], 0])
                local_axis = local_apex - np.array(part['center'])
                local_axis /= np.linalg.norm(local_axis)
                axis = rotate_vec(local_axis)
                
                span = abs(part['theta2'] - part['theta1'])
                cos_phi = np.cos(np.radians(span / 2.0))
                
                curves.append({'type': 'arc', 'center': center, 'radius': part['radius'],
                               'axis': axis, 'cos_half_angle': cos_phi, 'medium_id': medium_id})

        # 4. Generate Side Walls (Top Edge and Bottom Edge)
        # corners[0] is Front (Top, Bot), corners[1] is Back (Top, Bot)
        front_top, front_bot = corners[0]
        back_top, back_bot = corners[1]
        
        # Add Top Edge Line
        curves.append({
            'type': 'line', 
            'p1': front_top, 
            'p2': back_top, 
            'medium_id': medium_id
        })
        
        # Add Bottom Edge Line
        curves.append({
            'type': 'line', 
            'p1': front_bot, 
            'p2': back_bot, 
            'medium_id': medium_id
        })
            
        return curves

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

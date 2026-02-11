import numpy as np
from matplotlib.path import Path


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


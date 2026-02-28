import numpy as np
import pyvista as pv
import vtk

def compute_electrostatic_streamlines(scene_objects):
    charges = [obj for obj in scene_objects if obj.get_rasterization_data().get("type") == "charge"]
    if not charges: return None

    charges_data = []
    seed_points = []
    seed_dirs = []

    # 1. EXTRACT TRUE TRANSFORMED POSITIONS
    for obj in charges:
        q = obj.charge
        if q == 0: continue
        
        # Pull the exact World Center from the VTK matrix (Fixes Ghost Lines)
        matrix = vtk.vtkMatrix4x4()
        obj.transform.GetMatrix(matrix)
        m4 = np.zeros((4, 4))
        for i in range(4):
            for j in range(4):
                m4[i, j] = matrix.GetElement(i, j)
        
        b = obj.actors[0].mapper.dataset.bounds
        base_center = np.array([(b[0]+b[1])/2, (b[2]+b[3])/2, (b[4]+b[5])/2, 1.0])
        pos = (m4 @ base_center)[:3]
        
        charges_data.append((q, pos))
        
        # Create seeds just outside the visual sphere
        sphere = pv.Sphere(radius=0.6, center=pos, theta_resolution=6, phi_resolution=6)
        pts = sphere.points.tolist()
        seed_points.extend(pts)
        
        # Positive charges trace forward (+E), negative trace backward (-E)
        direction = 1 if q > 0 else -1
        seed_dirs.extend([direction] * len(pts))

    if not charges_data or not seed_points: return None

    current_pts = np.array(seed_points)
    seed_dirs = np.array(seed_dirs)[:, np.newaxis]
    
    # 2. VECTORIZED RK2 SOLVER (Gridless & Lightning Fast)
    lines = [[pt] for pt in current_pts]
    step_size = 0.2
    max_steps = 150
    
    def compute_E_direction(pts):
        E_total = np.zeros_like(pts)
        for q, p in charges_data:
            r = pts - p
            r_mag = np.linalg.norm(r, axis=1, keepdims=True)
            r_mag[r_mag < 0.2] = 0.2  # Prevent divide-by-zero singularities
            E_total += (q * r) / (r_mag**3) # Superposition of Coulomb's Law
        
        mags = np.linalg.norm(E_total, axis=1, keepdims=True)
        mags[mags < 1e-9] = 1e-9
        return E_total / mags

    # Tracking array to stop lines if they hit another charge
    active = np.ones(len(current_pts), dtype=bool)

    for _ in range(max_steps):
        if not np.any(active): break
            
        active_pts = current_pts[active]
        active_dirs = seed_dirs[active]
        
        # Standard Runge-Kutta 2 (Midpoint) Integration
        k1 = compute_E_direction(active_pts) * active_dirs
        k2 = compute_E_direction(active_pts + k1 * step_size * 0.5) * active_dirs
        
        new_pts = active_pts + k2 * step_size
        current_pts[active] = new_pts
        
        # Update lists and check for collisions
        active_indices = np.where(active)[0]
        for idx, pt in zip(active_indices, new_pts):
            lines[idx].append(pt.copy())
            
            # Stop tracing if it sinks into a charge
            for q, p in charges_data:
                if np.linalg.norm(pt - p) < 0.3:
                    active[idx] = False

    # 3. CONVERT TO PYVISTA SPLINES
    poly_lines = pv.PolyData()
    all_pts = []
    lines_cells = []
    
    pt_offset = 0
    for line_pts in lines:
        if len(line_pts) > 1:
            all_pts.extend(line_pts)
            lines_cells.append(len(line_pts))
            lines_cells.extend(range(pt_offset, pt_offset + len(line_pts)))
            pt_offset += len(line_pts)
            
    if len(all_pts) == 0: return None
    
    poly_lines.points = np.array(all_pts)
    poly_lines.lines = np.array(lines_cells)

    return poly_lines

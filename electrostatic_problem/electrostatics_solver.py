import numpy as np
import pyvista as pv
import vtk
from numba import njit

# ---------------------------------------------------------
# THE NUMBA COMPILED CORE (Runs at C++ speed)
# ---------------------------------------------------------
import numpy as np
import pyvista as pv
import vtk
from numba import njit

# --- 1. THE VECTORIZED E-FIELD HELPER ---
@njit(cache=True)
def calc_E_dir(pt, c_pos, c_q):
    """Calculates the normalized E-field direction using NumPy broadcasting."""
    # pt is (3,), c_pos is (N, 3) -> r broadcasts to (N, 3)
    r = pt - c_pos 
    
    # Calculate squared distance for all charges at once
    r2 = np.sum(r**2, axis=1)
    # Prevent singularities
    r2[r2 < 0.04] = 0.04
    r_mag = np.sqrt(r2)
    
    # Coulomb's law scalar factors: q / r^3
    factor = c_q / (r_mag * r2)
    
    # Broadcast factors to (N, 1), multiply by r vectors, and sum into a final (3,) vector
    E = np.sum(r * factor[:, np.newaxis], axis=0)
    
    # Normalize to get pure direction
    E_mag = np.linalg.norm(E)
    if E_mag < 1e-9: 
        E_mag = 1e-9
    
    return (E / E_mag)

# --- 2. THE COMPACT RK2 LOOP ---
@njit(cache=True)
def integrate_rk2(seed_pts, seed_dirs, c_pos, c_q, max_steps, step_size):
    num_lines = seed_pts.shape[0]
    trajectories = np.zeros((num_lines, max_steps, 3), dtype=np.float64)
    counts = np.zeros(num_lines, dtype=np.int32)
    
    for i in range(num_lines):
        pt = seed_pts[i].copy()
        dir_val = seed_dirs[i]
        
        trajectories[i, 0] = pt
        counts[i] = 1
        
        for step in range(1, max_steps):
            # 1. E-field at current point
            E1_dir = calc_E_dir(pt, c_pos, c_q)
            
            # 2. Scout the midpoint
            mid_pt = pt + E1_dir * step_size * 0.5 * dir_val
            
            # 3. E-field at midpoint
            E2_dir = calc_E_dir(mid_pt, c_pos, c_q)
            
            # 4. Take final step
            pt += E2_dir * step_size * dir_val
            
            trajectories[i, step] = pt
            counts[i] += 1
            
            # 5. Vectorized Collision Check: stop if we hit any charge
            dist_sq = np.sum((pt - c_pos)**2, axis=1)
            if np.any(dist_sq < 0.09):  # 0.3 squared
                break
                
    return trajectories, counts

# ---------------------------------------------------------
# THE PYTHON WRAPPER
# ---------------------------------------------------------
def compute_electrostatic_streamlines(scene_objects):
    charges = [obj for obj in scene_objects if obj.get_rasterization_data().get("type") == "charge"]
    if not charges: return None

    charge_pos = []
    charge_q = []
    seed_points = []
    seed_dirs = []

    # Extract coordinates
    for obj in charges:
        q = obj.charge
        if q == 0: continue
        
        matrix = vtk.vtkMatrix4x4()
        obj.transform.GetMatrix(matrix)
        m4 = np.zeros((4, 4))
        for i in range(4):
            for j in range(4):
                m4[i, j] = matrix.GetElement(i, j)
        
        b = obj.actors[0].mapper.dataset.bounds
        base_center = np.array([(b[0]+b[1])/2, (b[2]+b[3])/2, (b[4]+b[5])/2, 1.0])
        pos = (m4 @ base_center)[:3]
        
        charge_pos.append(pos)
        charge_q.append(q)
        
        # Lower resolution (4 instead of 6) for cleaner visuals
        sphere = pv.Sphere(radius=0.6, center=pos, theta_resolution=4, phi_resolution=4)
        pts = sphere.points.tolist()
        seed_points.extend(pts)
        
        direction = 1.0 if q > 0 else -1.0
        seed_dirs.extend([direction] * len(pts))

    if not charge_pos or not seed_points: return None

    # Convert to pure contiguous NumPy arrays for Numba
    c_pos_arr = np.array(charge_pos, dtype=np.float64)
    c_q_arr = np.array(charge_q, dtype=np.float64)
    seed_pts_arr = np.array(seed_points, dtype=np.float64)
    seed_dirs_arr = np.array(seed_dirs, dtype=np.float64)
    
    # Execute compiled C++ speed integration
    trajectories, counts = integrate_rk2(
        seed_pts_arr, seed_dirs_arr, c_pos_arr, c_q_arr, 
        max_steps=150, step_size=0.2
    )

    # Convert the pre-allocated block back into PyVista splines
    all_pts = []
    lines_cells = []
    pt_offset = 0
    
    for i in range(len(counts)):
        c = counts[i]
        if c > 1:
            # Only slice out the array steps that were actually used
            all_pts.extend(trajectories[i, :c])
            lines_cells.append(c)
            lines_cells.extend(range(pt_offset, pt_offset + c))
            pt_offset += c
            
    if not all_pts: return None
    
    poly_lines = pv.PolyData()
    poly_lines.points = np.array(all_pts)
    poly_lines.lines = np.array(lines_cells)

    return poly_lines

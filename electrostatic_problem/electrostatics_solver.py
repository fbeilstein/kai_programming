import numpy as np
import pyvista as pv
import vtk
from numba import njit
from implementation_tasks import propagate_trajectory



# =========================================================
# PART 2: THE ENGINE WRAPPERS (Students ignore this)
# =========================================================

@njit(cache=True)
def integrate_rk2(seed_pts, seed_dirs, charges_positions, charges_values, max_steps, step_size):
    """Iterates over all seed points and calls the student's trajectory function."""
    num_lines = seed_pts.shape[0]
    
    # PRE-ALLOCATION: Giant block of memory for all lines
    trajectories = np.zeros((num_lines, max_steps, 3), dtype=np.float64)
    counts = np.zeros(num_lines, dtype=np.int32)
    
    for i in range(num_lines):
        # Set the starting point for this specific line
        trajectories[i, 0] = seed_pts[i]
        current_charge_sign = seed_dirs[i]
        
        # --- Call the isolated student logic, passing the memory slice in-place ---
        count = propagate_trajectory(
            trajectories[i], 
            current_charge_sign, 
            charges_positions, 
            charges_values, 
            max_steps, 
            step_size
        )
        
        counts[i] = count
            
    return trajectories, counts

def compute_electrostatic_streamlines(scene_objects):
    """Extracts UI coordinates and converts numerical arrays to PyVista lines."""
    charges = [obj for obj in scene_objects if obj.get_rasterization_data().get("type") == "charge"]
    if not charges: return None

    charge_pos = []
    charge_q = []
    seed_points = []
    seed_dirs = []

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
        
        res = int(np.clip(abs(q) * 4, 3, 15))
        sphere = pv.Sphere(radius=0.6, center=pos, theta_resolution=res, phi_resolution=res)
        pts = sphere.points.tolist()
        seed_points.extend(pts)
        
        direction = 1.0 if q > 0 else -1.0
        seed_dirs.extend([direction] * len(pts))

    if not charge_pos or not seed_points: return None

    c_pos_arr = np.array(charge_pos, dtype=np.float64)
    c_q_arr = np.array(charge_q, dtype=np.float64)
    seed_pts_arr = np.array(seed_points, dtype=np.float64)
    seed_dirs_arr = np.array(seed_dirs, dtype=np.float64)
    
    trajectories, counts = integrate_rk2(
        seed_pts_arr, seed_dirs_arr, c_pos_arr, c_q_arr, 
        max_steps=150, step_size=0.2
    )

    all_pts = []
    lines_cells = []
    pt_offset = 0
    
    for i in range(len(counts)):
        c = counts[i]
        if c > 1:
            # Slicing the exact number of populated steps
            all_pts.extend(trajectories[i, :c])
            lines_cells.append(c)
            lines_cells.extend(range(pt_offset, pt_offset + c))
            pt_offset += c
            
    if not all_pts: return None
    
    poly_lines = pv.PolyData()
    poly_lines.points = np.array(all_pts)
    poly_lines.lines = np.array(lines_cells)

    return poly_lines

import numpy as np
import pyvista as pv

def compute_electrostatic_streamlines(scene_objects, bounds=(-10, 10, -10, 10, -10, 10), res=30):
    """
    Calculates the exact electric field using Coulomb's superposition principle
    and returns PyVista streamlines.
    """
    # Filter out everything except charges
    charges = [obj for obj in scene_objects if obj.get_rasterization_data().get("type") == "charge"]
    
    if not charges:
        return None

    # 1. Create a lightweight, fixed evaluation grid
    x = np.linspace(bounds[0], bounds[1], res)
    y = np.linspace(bounds[2], bounds[3], res)
    z = np.linspace(bounds[4], bounds[5], res)
    X, Y, Z = np.meshgrid(x, y, z, indexing='ij')
    
    grid = pv.StructuredGrid(X, Y, Z)
    pts = grid.points
    
    E_total = np.zeros_like(pts, dtype=np.float64)
    seed_points = []

    # 2. Vectorized Coulomb's Law (Superposition)
    for charge_obj in charges:
        pos = np.array(charge_obj.actors[0].GetCenter())
        q = charge_obj.charge
        
        if q == 0:
            continue

        # r vector from the charge to every point in space
        r = pts - pos
        r_mag = np.linalg.norm(r, axis=1)[:, np.newaxis]
        
        # Prevent division-by-zero explosions if a grid point is exactly inside the charge
        r_mag[r_mag < 0.2] = 0.2 
        
        # E = q * r / |r|^3
        E_total += (q * r) / (r_mag**3)
        
        # Drop streamline seeds directly around the charge's surface
        sphere = pv.Sphere(radius=0.7, center=pos, theta_resolution=8, phi_resolution=8)
        seed_points.extend(sphere.points.tolist())

    if len(seed_points) == 0:
        return None

    grid["E_field"] = E_total
    grid.set_active_vectors("E_field")

    # 3. Trace Streamlines
    seed_poly = pv.PolyData(np.array(seed_points))
    streamlines = grid.streamlines_from_source(
        seed_poly,
        vectors="E_field",
        max_length=40.0,
        integration_direction="both"
    )
    
    return streamlines

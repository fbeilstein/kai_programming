import numpy as np
import vtk
import pyvista as pv

def rasterize_scene(scene_objects, resolution=200, padding=4.0):
    print("1a. Auto-calculating grid bounds...")
    
    # 1. Find the absolute physical limits of all relevant objects
    mins = []
    maxs = []
    has_solvable_objects = False
    
    for obj in scene_objects:
        if type(obj).__name__ == "Cannon": 
            continue # Ignore cannon for magnetic bounds
            
        has_solvable_objects = True
        # Ensure we capture all parts of the object (like both poles of a magnet)
        for actor in obj.actors:
            b = actor.GetBounds() # VTK returns (xmin, xmax, ymin, ymax, zmin, zmax)
            mins.append([b[0], b[2], b[4]])
            maxs.append([b[1], b[3], b[5]])
        
    if not has_solvable_objects:
        center = np.array([0.0, 0.0, 0.0])
        bounds_half = 5.0
    else:
        mins = np.min(mins, axis=0)
        maxs = np.max(maxs, axis=0)
        center = (maxs + mins) / 2.0
        
        # 2. Force a cubic grid by finding the largest dimension
        max_span = np.max(maxs - mins) / 2.0
        bounds_half = max_span + padding 
        
    print(f"    -> Center: [{center[0]:.1f}, {center[1]:.1f}, {center[2]:.1f}]")
    print(f"    -> Box Width: {bounds_half * 2:.1f} units")

    # 3. Build the strictly cubic coordinates
    x = np.linspace(center[0] - bounds_half, center[0] + bounds_half, resolution)
    y = np.linspace(center[1] - bounds_half, center[1] + bounds_half, resolution)
    z = np.linspace(center[2] - bounds_half, center[2] + bounds_half, resolution)
    
    X, Y, Z = np.meshgrid(x, y, z, indexing='ij')
    
    # --- THE FIX: Create the flattened homogeneous points array ---
    points = np.column_stack((X.ravel(), Y.ravel(), Z.ravel()))
    points_h = np.column_stack((points, np.ones(points.shape[0])))
    
    # Initialize FLAT grids for easy boolean masking
    n_points = resolution ** 3
    mu_grid = np.ones(n_points, dtype=np.float32)
    M_grid = np.zeros((n_points, 3), dtype=np.float32)

    print("1b. Sampling object properties onto grid...")
    # 4. Populate the grids
    for obj in scene_objects:
        if type(obj).__name__ == "Cannon":
            continue 
            
        data = obj.get_rasterization_data()

        # Extract the 4x4 transform matrix to NumPy
        matrix = vtk.vtkMatrix4x4()
        data["transform"].GetMatrix(matrix)
        m4 = np.zeros((4, 4))
        for i in range(4):
            for j in range(4):
                m4[i, j] = matrix.GetElement(i, j)
        
        # Invert the matrix to map world points -> local object space
        m4_inv = np.linalg.inv(m4)

        # Transform the entire grid at once using matrix multiplication
        local_pts_h = points_h @ m4_inv.T
        local_pts = local_pts_h[:, :3] / local_pts_h[:, 3:] # Perspective divide

        if data["type"] == "material":
            # Get original untransformed bounds of the mesh
            bounds = obj.actors[0].mapper.dataset.bounds
            
            # Boolean mask: Which grid points fall inside these local bounds?
            mask = (
                (local_pts[:, 0] >= bounds[0]) & (local_pts[:, 0] <= bounds[1]) &
                (local_pts[:, 1] >= bounds[2]) & (local_pts[:, 1] <= bounds[3]) &
                (local_pts[:, 2] >= bounds[4]) & (local_pts[:, 2] <= bounds[5])
            )
            mu_grid[mask] = data["mu"]

        elif data["type"] == "magnet":
            # For magnets, mask both the North (red) and South (blue) poles
            bounds_n = obj.actors[0].mapper.dataset.bounds
            bounds_s = obj.actors[1].mapper.dataset.bounds
            
            mask_n = (
                (local_pts[:, 0] >= bounds_n[0]) & (local_pts[:, 0] <= bounds_n[1]) &
                (local_pts[:, 1] >= bounds_n[2]) & (local_pts[:, 1] <= bounds_n[3]) &
                (local_pts[:, 2] >= bounds_n[4]) & (local_pts[:, 2] <= bounds_n[5])
            )
            mask_s = (
                (local_pts[:, 0] >= bounds_s[0]) & (local_pts[:, 0] <= bounds_s[1]) &
                (local_pts[:, 1] >= bounds_s[2]) & (local_pts[:, 1] <= bounds_s[3]) &
                (local_pts[:, 2] >= bounds_s[4]) & (local_pts[:, 2] <= bounds_s[5])
            )
            full_mask = mask_n | mask_s

            # Calculate the global Magnetization vector M
            local_M = np.array([data["strength"], 0.0, 0.0, 0.0])
            
            # Rotate M to global space
            global_M = (m4 @ local_M)[:3]
            
            M_grid[full_mask] = global_M

    # 5. Reshape back to 3D spatial dimensions for the solver
    mu_grid = mu_grid.reshape((resolution, resolution, resolution))
    M_grid = M_grid.reshape((resolution, resolution, resolution, 3))

    print("Rasterization complete!")
    return X, Y, Z, mu_grid, M_grid

def debug_visualize_grid(X, Y, Z, mu_grid, M_grid):
    """Creates a quick point cloud to visualize the rasterized grid data."""
    points = np.column_stack((X.ravel(), Y.ravel(), Z.ravel()))
    mu_flat = mu_grid.ravel()
    
    M_mag_flat = np.linalg.norm(M_grid, axis=-1).ravel()
    
    cloud = pv.PolyData(points)
    cloud["Permeability"] = mu_flat
    cloud["Magnetization"] = M_mag_flat
    
    cloud["Solid"] = ((mu_flat > 1.01) | (M_mag_flat > 0.0)).astype(float)
    
    solid_objects = cloud.threshold(0.5, scalars="Solid")
    
    if solid_objects.n_points == 0:
        print("WARNING: Grid is empty! Objects might be outside the grid bounds.")
        return
        
    p = pv.Plotter()
    p.add_mesh(solid_objects, scalars="Permeability", cmap="plasma", point_size=6, render_points_as_spheres=True)
    p.add_axes()
    p.show(title="Rasterization Debug View")

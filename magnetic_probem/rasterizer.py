import numpy as np
import vtk
import pyvista as pv

def rasterize_scene(scene_objects, grid_bounds=(-8, 8, -8, 8, -8, 8), resolution=100):
    """
    Converts the analytical scene into 3D NumPy arrays for physics simulation.
    Returns: X, Y, Z, mu_grid (permeability), M_grid (magnetization vectors)
    """
    print(f"Rasterizing grid at {resolution}^3 resolution...")
    
    # 1. Generate the 3D grid
    x = np.linspace(grid_bounds[0], grid_bounds[1], resolution)
    y = np.linspace(grid_bounds[2], grid_bounds[3], resolution)
    z = np.linspace(grid_bounds[4], grid_bounds[5], resolution)
    X, Y, Z = np.meshgrid(x, y, z, indexing='ij')

    # Flatten into a list of points: shape (N, 3)
    points = np.column_stack((X.ravel(), Y.ravel(), Z.ravel()))
    N = len(points)

    # Convert to homogeneous coordinates for matrix multiplication: shape (N, 4)
    points_h = np.column_stack((points, np.ones(N)))

    # 2. Initialize physical property arrays
    mu_grid = np.ones(N, dtype=np.float32)       # Air has relative permeability of 1.0
    M_grid = np.zeros((N, 3), dtype=np.float32)  # Magnetization vectors (Mx, My, Mz)

    # 3. Populate the grids
    for obj in scene_objects:
        data = obj.get_rasterization_data()
        
        # We don't map the Cannon into the static magnetic field properties
        if data["type"] == "cannon":
            continue 

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
            # In our local space setup, North is at +X and South is at -X.
            # Therefore, M points along the positive local X axis.
            local_M = np.array([data["strength"], 0.0, 0.0, 0.0])
            
            # Rotate M to global space (ignore translation, hence the 0.0)
            global_M = (m4 @ local_M)[:3]
            
            M_grid[full_mask] = global_M

    # 4. Reshape back to 3D spatial dimensions
    mu_grid = mu_grid.reshape((resolution, resolution, resolution))
    M_grid = M_grid.reshape((resolution, resolution, resolution, 3))

    print("Rasterization complete!")
    return X, Y, Z, mu_grid, M_grid

def debug_visualize_grid(X, Y, Z, mu_grid, M_grid):
    """Creates a quick point cloud to visualize the rasterized grid data."""
    # Flatten the 3D arrays back into 1D lists for the point cloud
    points = np.column_stack((X.ravel(), Y.ravel(), Z.ravel()))
    mu_flat = mu_grid.ravel()
    
    # Calculate the magnitude of the magnetization vectors
    M_mag_flat = np.linalg.norm(M_grid, axis=-1).ravel()
    
    # Create the point cloud
    cloud = pv.PolyData(points)
    cloud["Permeability"] = mu_flat
    cloud["Magnetization"] = M_mag_flat
    
    # Define a "Solid": high permeability OR actively magnetized
    cloud["Solid"] = ((mu_flat > 1.01) | (M_mag_flat > 0.0)).astype(float)
    
    # Extract only the solid points
    solid_objects = cloud.threshold(0.5, scalars="Solid")
    
    # SAFETY CATCH: Prevent the PyVista zero-point crash
    if solid_objects.n_points == 0:
        print("WARNING: Grid is empty! Objects might be outside the grid bounds (-4 to 4).")
        return
        
    p = pv.Plotter()
    # Color by permeability, but magnets will also be visible
    p.add_mesh(solid_objects, scalars="Permeability", cmap="plasma", point_size=6, render_points_as_spheres=True)
    p.add_axes()
    p.show(title="Rasterization Debug View")

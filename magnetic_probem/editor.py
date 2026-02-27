import pyvista as pv
import vtk
from magnet import Magnet
from material import Material
from cannon import Cannon
from rasterizer import rasterize_scene, debug_visualize_grid
from solver import solve_magnetic_field
import numpy as np

# -----------------------------
# Main Setup
# -----------------------------
plotter = pv.Plotter()
plotter.enable_depth_peeling()

scene_objects = []
selected = {"current": None}

# -----------------------------
# UI Callbacks
# -----------------------------
def on_strength_change(value):
    if selected["current"] and isinstance(selected["current"], Magnet):
        selected["current"].strength = value

def on_mu_change(value):
    if selected["current"] and isinstance(selected["current"], Material):
        selected["current"].mu = value

def toggle_cannon(state):
    for obj in scene_objects:
        if isinstance(obj, Cannon):
            obj.is_firing = state
            if state: obj.actors[1].GetProperty().SetColor(1.0, 1.0, 0.0) 
            else: obj.actors[1].GetProperty().SetColor(0.0, 0.39, 0.0)
            plotter.render()

def spawn_magnet(state):
    offset = len(scene_objects) * 0.1
    new_mag = Magnet(plotter, center=(offset, offset, 0))
    scene_objects.append(new_mag)
    plotter.reset_camera_clipping_range()
    plotter.render()

def spawn_material(state):
    offset = len(scene_objects) * 0.1
    new_mat = Material(plotter, center=(offset, offset, 0))
    scene_objects.append(new_mat)
    plotter.reset_camera_clipping_range()
    plotter.render()

def compute_and_show_grid(state=None):
    if len(scene_objects) == 0:
        print("Scene is empty. Add some objects first!")
        return
        
    print("Extracting physical properties to grid...")
    X, Y, Z, mu_grid, M_grid = rasterize_scene(scene_objects)
    
    # THE FIX: Pass M_grid to the visualizer!
    debug_visualize_grid(X, Y, Z, mu_grid, M_grid)

def simulate_and_draw_lines(state=None):
    if len(scene_objects) == 0:
        print("Scene is empty!")
        return
        
    print("1. Rasterizing grid...")
    X, Y, Z, mu_grid, M_grid = rasterize_scene(scene_objects)
    
    print("2. Solving magnetic field on GPU...")
    B_grid = solve_magnetic_field(X, mu_grid, M_grid)
    
    print("3. Tracing streamlines...")
    # Rebuild the spatial grid in PyVista
    grid = pv.StructuredGrid(X, Y, Z)
    
    # --- THE FIX: Unscramble the matrices ---
    # Force NumPy to flatten the arrays using Fortran-order ('F')
    # so they perfectly match PyVista's physical coordinates.
    Bx_flat = B_grid[..., 0].flatten(order='F')
    By_flat = B_grid[..., 1].flatten(order='F')
    Bz_flat = B_grid[..., 2].flatten(order='F')
    
    grid["B_field"] = np.column_stack((Bx_flat, By_flat, Bz_flat))
    grid.set_active_vectors("B_field")
    
    # Generate seed points around all magnets
    seed_points = pv.PolyData()
    has_magnet = False
    for obj in scene_objects:
        if isinstance(obj, Magnet):
            has_magnet = True
            # Get the actual world position of the magnet to center the seed sphere
            b = obj.actors[0].GetBounds()
            center = [(b[0]+b[1])/2, (b[2]+b[3])/2, (b[4]+b[5])/2]
            
            # High-density sphere to spawn lots of lines
            sphere = pv.Sphere(radius=1.0, 
                                center=center, 
                                theta_resolution=15, 
                                phi_resolution=15)
            seed_points = seed_points.merge(sphere)
            
    if not has_magnet:
         print("No magnets to seed lines from! Add at least one magnet.")
         return

    # Trace the physics!
    streamlines = grid.streamlines_from_source(
        seed_points,
        vectors="B_field",
        max_length=50.0,  # CRITICAL: Let the lines travel across the new giant box
        initial_step_length=0.1,
        integration_direction="both"
    )
    
    # --- Show the results ---
    p = pv.Plotter()
    
    # Draw the original meshes as semi-transparent ghosts for reference
    for obj in scene_objects:
        # Skip drawing the cannon in the static magnetic field view
        if isinstance(obj, Cannon):
            continue
            
        for actor in obj.actors:
            # Extract the raw mesh and apply the current VTK transform
            mesh = actor.mapper.dataset.copy()
            mesh.transform(obj.transform)
            p.add_mesh(mesh, color=actor.prop.color, opacity=0.3)
            
    # Draw the streamlines as colored tubes based on field strength
    p.add_mesh(streamlines.tube(radius=0.005), cmap="plasma", render_lines_as_tubes=False)
    p.add_axes()
    p.show(title="Magnetic Field Simulation")

def delete_selected(state=None):
    if selected["current"] is not None:
        obj_to_delete = selected["current"]
        
        # --- THE FIX: Shield the Cannon ---
        if isinstance(obj_to_delete, Cannon):
            print("Cannot delete the Cannon!") # Optional console warning
            return 
        
        obj_to_delete.destroy()                
        scene_objects.remove(obj_to_delete)    
        
        selected["current"] = None             
        strength_slider.Off()                  
        mu_slider.Off()
        plotter.render()

# -----------------------------
# Build Contextual UI
# -----------------------------
strength_slider = plotter.add_slider_widget(callback=on_strength_change, rng=[0.1, 5.0], value=1.0, title="Magnet Strength (T)", pointa=(0.65, 0.9), pointb=(0.95, 0.9), style="modern")
strength_slider.Off() 

mu_slider = plotter.add_slider_widget(callback=on_mu_change, rng=[1.0, 10000.0], value=1000.0, title="Permeability (mu_r)", pointa=(0.65, 0.9), pointb=(0.95, 0.9), style="modern")
mu_slider.Off()

# --- NEW DELETE BUTTON AND KEYBIND ---
plotter.add_checkbox_button_widget(delete_selected, position=(10, 155), size=30, color_on="red", color_off="red")
plotter.add_text("Delete Selected", position=(50, 160), font_size=10)
plotter.add_key_event("Delete", delete_selected) # Bind the physical Delete key

plotter.add_checkbox_button_widget(spawn_magnet, position=(10, 110), size=30, color_on="orange", color_off="orange")
plotter.add_text("Add Magnet", position=(50, 115), font_size=10)

plotter.add_checkbox_button_widget(spawn_material, position=(10, 65), size=30, color_on="silver", color_off="silver")
plotter.add_text("Add Material", position=(50, 70), font_size=10)

plotter.add_checkbox_button_widget(toggle_cannon, position=(10, 20), size=30, color_on="yellow", color_off="darkgreen")
plotter.add_text("Toggle Cannon", position=(50, 25), font_size=10)

plotter.add_checkbox_button_widget(compute_and_show_grid, position=(10, 200), size=30, color_on="cyan", color_off="cyan")
plotter.add_text("Compute Grid", position=(50, 205), font_size=10)

plotter.add_checkbox_button_widget(simulate_and_draw_lines, position=(10, 245), size=30, color_on="green", color_off="green")
plotter.add_text("Simulate Field", position=(50, 250), font_size=10)

# -----------------------------
# Custom Selection Engine
# -----------------------------
def on_left_click(obj, event):
    click_pos = plotter.iren.get_event_position()
    picker = vtk.vtkPropPicker()
    picker.PickProp(click_pos[0], click_pos[1], plotter.renderer)
    picked_actor = picker.GetViewProp()
    
    if picked_actor is None:
        if selected["current"] is not None:
            selected["current"].set_active(False)
            selected["current"] = None
            strength_slider.Off()
            mu_slider.Off()
        return

    for scene_obj in scene_objects:
        if picked_actor in scene_obj.actors:
            if selected["current"] == scene_obj:
                return 
            
            if selected["current"] is not None:
                selected["current"].set_active(False)
            
            selected["current"] = scene_obj
            scene_obj.set_active(True)
            
            strength_slider.Off()
            mu_slider.Off()
            
            if isinstance(scene_obj, Magnet):
                strength_slider.GetRepresentation().SetValue(scene_obj.strength)
                strength_slider.On()
            elif isinstance(scene_obj, Material):
                mu_slider.GetRepresentation().SetValue(scene_obj.mu)
                mu_slider.On()
            return

plotter.iren.add_observer("LeftButtonPressEvent", on_left_click)

scene_objects.append(Cannon(plotter, center=(-3, 0, 0)))

plotter.add_axes()
plotter.show()

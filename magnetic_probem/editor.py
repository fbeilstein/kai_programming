import pyvista as pv
import vtk
from magnet import Magnet
from material import Material
from cannon import Cannon

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

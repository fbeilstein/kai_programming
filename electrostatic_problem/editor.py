import pyvista as pv
import vtk
import numpy as np

from charge import ElectricCharge
from electrostatics_solver import compute_electrostatic_streamlines

# -----------------------------
# Main Setup
# -----------------------------
plotter = pv.Plotter()
plotter.enable_depth_peeling()

scene_objects = []
selected = {"current": None}

# -----------------------------
# The Live-Update Engine
# -----------------------------
def update_field(*args):
    lines = compute_electrostatic_streamlines(scene_objects)
    
    if lines is not None and lines.n_points > 0:
        plotter.add_mesh(lines.tube(radius=0.01), color="gold", name="E_field", reset_camera=False)
    else:
        try:
            plotter.remove_actor("E_field", render=False)
        except Exception:
            pass
            
    plotter.render()

# -----------------------------
# UI Callbacks
# -----------------------------
def on_charge_change(value):
    if selected["current"] and isinstance(selected["current"], ElectricCharge):
        selected["current"].charge = value
        update_field()

def spawn_charge(state):
    offset = len(scene_objects) * 1.5 
    new_charge = ElectricCharge(plotter, center=(offset, offset, 0), charge=1.0)
    
    new_charge.widget.AddObserver("InteractionEvent", lambda w, e: update_field())
    scene_objects.append(new_charge)
    
    if selected["current"] is not None:
        selected["current"].set_active(False)
    
    selected["current"] = new_charge
    new_charge.set_active(True)
    
    charge_slider.GetRepresentation().SetValue(1.0)
    charge_slider.On()
    
    update_field()
    plotter.reset_camera_clipping_range()

def delete_selected(state=None):
    if selected["current"] is not None:
        obj_to_delete = selected["current"]
        
        obj_to_delete.destroy()                
        scene_objects.remove(obj_to_delete)    
        
        selected["current"] = None             
        charge_slider.Off()                  
        
        update_field()

# -----------------------------
# Build Contextual UI
# -----------------------------
charge_slider = plotter.add_slider_widget(
    callback=on_charge_change, 
    rng=[-5.0, 5.0], 
    value=1.0, 
    title="Charge (Coulombs)", 
    pointa=(0.65, 0.9), 
    pointb=(0.95, 0.9), 
    style="modern",
    interaction_event="always"
)
charge_slider.Off() 

plotter.add_checkbox_button_widget(delete_selected, position=(10, 65), size=30, color_on="red", color_off="red")
plotter.add_text("Delete Selected", position=(50, 70), font_size=10)
plotter.add_key_event("Delete", delete_selected) 

plotter.add_checkbox_button_widget(spawn_charge, position=(10, 20), size=30, color_on="blue", color_off="red")
plotter.add_text("Add Point Charge", position=(50, 25), font_size=10)

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
            charge_slider.Off()
        return

    for scene_obj in scene_objects:
        if picked_actor in scene_obj.actors:
            if selected["current"] == scene_obj:
                return 
            
            if selected["current"] is not None:
                selected["current"].set_active(False)
            
            selected["current"] = scene_obj
            scene_obj.set_active(True)
            
            charge_slider.Off()
            
            if isinstance(scene_obj, ElectricCharge):
                charge_slider.GetRepresentation().SetValue(scene_obj.charge)
                charge_slider.On()
            return

plotter.iren.add_observer("LeftButtonPressEvent", on_left_click)

plotter.add_axes()
plotter.show()

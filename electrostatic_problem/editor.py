import pyvista as pv
import vtk
import numpy as np

from charge import ElectricCharge
from cannon import Cannon
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
    """Calculates and redraws the field lines instantly."""
    lines = compute_electrostatic_streamlines(scene_objects)
    
    if lines is not None and lines.n_points > 0:
        plotter.add_mesh(lines.tube(radius=0.01), color="gold", name="E_field", reset_camera=False)
    else:
        # Prevent crashing if all charges are deleted or set to 0
        try:
            plotter.remove_actor("E_field", render=False)
        except Exception:
            pass
            
    plotter.render() # <-- This forces the instant UI update!

# -----------------------------
# UI Callbacks
# -----------------------------
def on_charge_change(value):
    if selected["current"] and isinstance(selected["current"], ElectricCharge):
        selected["current"].charge = value
        update_field()

def toggle_cannon(state):
    for obj in scene_objects:
        if isinstance(obj, Cannon):
            obj.is_firing = state
            if state: 
                obj.actors[1].GetProperty().SetColor(1.0, 1.0, 0.0) 
            else: 
                obj.actors[1].GetProperty().SetColor(0.0, 0.39, 0.0)
            plotter.render()

def spawn_charge(state):
    offset = len(scene_objects) * 1.5 
    new_charge = ElectricCharge(plotter, center=(offset, offset, 0), charge=1.0)
    
    # Bind the physical drag interaction directly to the live update
    new_charge.widget.AddObserver("InteractionEvent", lambda w, e: update_field())
    scene_objects.append(new_charge)
    
    # QoL FIX: Automatically select and enable the widget for the new charge
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
        
        if isinstance(obj_to_delete, Cannon):
            print("Cannot delete the Cannon!") 
            return 
        
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

plotter.add_checkbox_button_widget(delete_selected, position=(10, 110), size=30, color_on="red", color_off="red")
plotter.add_text("Delete Selected", position=(50, 115), font_size=10)
plotter.add_key_event("Delete", delete_selected) 

plotter.add_checkbox_button_widget(spawn_charge, position=(10, 65), size=30, color_on="blue", color_off="red")
plotter.add_text("Add Point Charge", position=(50, 70), font_size=10)

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

# Spawn the permanent cannon
scene_objects.append(Cannon(plotter, center=(-5, 0, 0)))

plotter.add_axes()
plotter.show()

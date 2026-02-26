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

# Initialize our physical objects
scene_objects = [
    Cannon(plotter, center=(-2, 0, 0)),
    Material(plotter, center=(0, 0, 0), mu=5000),
    Magnet(plotter, center=(2, 0, 0))
]

selected = {"current": None}

# -----------------------------
# Custom Selection Engine
# -----------------------------
def on_left_click(obj, event):
    click_pos = plotter.iren.get_event_position()
    picker = vtk.vtkPropPicker()
    picker.PickProp(click_pos[0], click_pos[1], plotter.renderer)
    picked_actor = picker.GetViewProp()
    
    # Clicked empty space
    if picked_actor is None:
        if selected["current"] is not None:
            selected["current"].set_active(False)
            selected["current"] = None
        return

    # Find which logical object owns the clicked actor
    for scene_obj in scene_objects:
        if picked_actor in scene_obj.actors:
            if selected["current"] == scene_obj:
                return 
            
            if selected["current"] is not None:
                selected["current"].set_active(False)
            
            selected["current"] = scene_obj
            scene_obj.set_active(True)
            return

plotter.iren.add_observer("LeftButtonPressEvent", on_left_click)

plotter.add_axes()
plotter.show()

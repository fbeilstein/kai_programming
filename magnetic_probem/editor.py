import pyvista as pv
import vtk
import numpy as np

plotter = pv.Plotter()
plotter.background_color = "white"

# --- Setup Meshes ---
# Red at +0.5, Blue at -0.5. The center of the pair is (0,0,0)
north = pv.Cube(center=(0.5, 0, 0))
south = pv.Cube(center=(-0.5, 0, 0))

n_actor = plotter.add_mesh(north, color="red", opacity=0.7)
s_actor = plotter.add_mesh(south, color="blue", opacity=0.7)
actors = [n_actor, s_actor]

# This is our persistent "Source of Truth" matrix
# It starts as Identity because the magnets are already at their starting spots
current_master_matrix = vtk.vtkMatrix4x4()

# --- Configure the Widget ---
vtk_box = vtk.vtkBoxWidget()
vtk_box.SetInteractor(plotter.iren.interactor)
vtk_box.SetPlaceFactor(1.0)

# Place the widget once on the combined geometry
b1, b2 = north.bounds, south.bounds
init_bounds = (
    min(b1[0], b2[0]), max(b1[1], b2[1]),
    min(b1[2], b2[2]), max(b1[3], b2[3]),
    min(b1[4], b2[4]), max(b1[5], b2[5])
)
vtk_box.PlaceWidget(init_bounds)

# REQUIRED: Explicitly enable rotation
vtk_box.RotationEnabledOn()
vtk_box.TranslationEnabledOn()
vtk_box.ScalingEnabledOn()

# Hide visual junk
vtk_box.OutlineCursorWiresOff()
vtk_box.GetOutlineProperty().SetOpacity(0)
vtk_box.GetFaceProperty().SetOpacity(0)

def box_widget_callback(obj, event):
    # 1. Get the raw matrix from the widget
    raw_tf = vtk.vtkTransform()
    obj.GetTransform(raw_tf)
    raw_mat = raw_tf.GetMatrix()
    
    # 2. Extract Translation (T) and Rotation/Scale (RS)
    # We apply the formula: Final = T * R * T_inv
    # But since we want to rotate around the CENTER of the object's current position:
    # We simply apply the raw_mat directly to the actors because 
    # the widget's internal 'PlaceWidget' already set the pivot to the center of init_bounds.
    
    for a in actors:
        a.SetUserMatrix(raw_mat)
    
    current_master_matrix.DeepCopy(raw_mat)
    print(current_master_matrix)
    plotter.render()

vtk_box.AddObserver("InteractionEvent", box_widget_callback)

def handle_click(obj, event):
    print("fdfdfd")
    click_pos = plotter.iren.interactor.GetEventPosition()
    picker = vtk.vtkPropPicker()
    picker.Pick(click_pos[0], click_pos[1], 0, plotter.renderer)
    picked_actor = picker.GetActor()

    if picked_actor in actors:
        # Sync the gizmo to our master matrix
        temp_tf = vtk.vtkTransform()
        temp_tf.SetMatrix(current_master_matrix)
        print(current_master_matrix)
        vtk_box.SetTransform(temp_tf)
        vtk_box.SetEnabled(1)
    elif picked_actor is None:
        vtk_box.SetEnabled(0)
    
    plotter.render()

plotter.iren.interactor.AddObserver("LeftButtonPressEvent", handle_click, 1.0)
plotter.iren.interactor.AddObserver("RightButtonPressEvent", handle_click, 1.0)

print("Formula applied: Rotation is now forced around the center of mass.")
plotter.show()

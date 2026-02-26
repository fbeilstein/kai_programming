import pyvista as pv
import numpy as np
import vtk

# --- Helper to convert numpy 4x4 to vtkMatrix4x4 ---
def numpy_to_vtk_matrix(mat):
    vtk_mat = vtk.vtkMatrix4x4()
    for i in range(4):
        for j in range(4):
            vtk_mat.SetElement(i, j, float(mat[i, j]))
    return vtk_mat

# --- Create a PyVista plotter ---
plotter = pv.Plotter()
plotter.background_color = "white"

# --- Add two cubes (magnets) ---
north = pv.Cube(center=(0.5, 0, 0))
south = pv.Cube(center=(-0.5, 0, 0))
north_actor = plotter.add_mesh(north, color="red")
south_actor = plotter.add_mesh(south, color="blue")
actors = [north_actor, south_actor]

# --- Compute initial bounding box ---
def get_bounds():
    b = list(north.bounds)
    for i in range(3):
        b[2*i] = min(b[2*i], south.bounds[2*i])
        b[2*i+1] = max(b[2*i+1], south.bounds[2*i+1])
    return tuple(b)

initial_bounds = get_bounds()

def expand_bounds(bounds, factor=1.3):
    xmin, xmax, ymin, ymax, zmin, zmax = bounds
    cx = (xmin + xmax)/2
    cy = (ymin + ymax)/2
    cz = (zmin + zmax)/2
    dx = (xmax - xmin) * factor / 2
    dy = (ymax - ymin) * factor / 2
    dz = (zmax - zmin) * factor / 2
    return (cx - dx, cx + dx, cy - dy, cy + dy, cz - dz, cz + dz)

# --- Callback to apply vtkTransform to actors ---
def apply_transform(vtk_transform):
    mat = vtk_transform.GetMatrix()
    for a in actors:
        a.SetUserMatrix(mat)
    plotter.render()

# --- Create raw VTK box widget ---
# Access the underlying interactor
interactor = plotter.ren_win.GetInteractor()

vtk_box = vtk.vtkBoxWidget()
vtk_box.SetInteractor(interactor)  # now works
vtk_box.SetPlaceFactor(1.3)
vtk_box.PlaceWidget(initial_bounds)
vtk_box.HandlesOn()
vtk_box.On()
vtk_box.EnabledOff()  # start hidden

print(f"{dir(vtk_box)}")

# --- Hide the box outline, keep handles ---
#rep = vtk_box.GetRepresentation()
#rep.GetOutlineProperty().SetOpacity(0.0)  # make the cube invisible
#rep.GetOutlineProperty().SetColor(0,0,0)  # optional, black if it appears
#rep.HandlesOn()  # ensure handles still visible


# --- Callback when the box is interacted with ---
def box_widget_callback(obj, event):
    transform = vtk.vtkTransform()
    obj.GetTransform(transform)
    apply_transform(transform)

vtk_box.AddObserver("InteractionEvent", box_widget_callback)

# --- Picking callback to show/hide box widget ---
def pick_callback(picked):
    if picked in actors:
        vtk_box.PlaceWidget(expand_bounds(get_bounds(), 1.3))
        vtk_box.EnabledOn()
    else:
        vtk_box.EnabledOff()

plotter.enable_mesh_picking(callback=pick_callback, use_actor=True)

plotter.show()

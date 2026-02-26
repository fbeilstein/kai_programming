import pyvista as pv
import vtk

# -----------------------------
# Scene Setup
# -----------------------------
plotter = pv.Plotter()

# FIX: Enable depth peeling. This is crucial for VTK to render translucent widget faces properly.
plotter.enable_depth_peeling()

def create_magnet_actors(center):
    cubes = []
    offsets = [(0, 0, 0), (0.4, 0, 0), (-0.4, 0, 0)]
    for dx, dy, dz in offsets:
        cube = pv.Cube(
            center=(center[0]+dx, center[1]+dy, center[2]+dz), 
            x_length=0.4, y_length=0.4, z_length=0.4
        )
        actor = plotter.add_mesh(
            cube,
            color="orange",
            opacity=0.35,
            show_edges=True,
            edge_color="black",
            line_width=1
        )
        cubes.append(actor)
    return cubes

# -----------------------------
# Stateful Magnet Class
# -----------------------------
class MagnetGroup:
    def __init__(self, actors):
        self.actors = actors
        self.transform = vtk.vtkTransform()
        
        bounds = [None]*6
        for i, actor in enumerate(self.actors):
            b = actor.GetBounds()
            if i == 0:
                bounds = list(b)
            else:
                bounds = [
                    min(bounds[0], b[0]), max(bounds[1], b[1]),
                    min(bounds[2], b[2]), max(bounds[3], b[3]),
                    min(bounds[4], b[4]), max(bounds[5], b[5])
                ]
        
        scale = 1.1
        cx, cy, cz = (bounds[0]+bounds[1])/2, (bounds[2]+bounds[3])/2, (bounds[4]+bounds[5])/2
        dx = (bounds[1]-bounds[0]) * scale
        dy = (bounds[3]-bounds[2]) * scale
        dz = (bounds[5]-bounds[4]) * scale
        
        pad_bounds = [
            cx - (dx if dx else 0.1), cx + (dx if dx else 0.1),
            cy - (dy if dy else 0.1), cy + (dy if dy else 0.1),
            cz - (dz if dz else 0.1), cz + (dz if dz else 0.1)
        ]

        self.widget = vtk.vtkBoxWidget2()
        self.rep = vtk.vtkBoxRepresentation()
        self.rep.PlaceWidget(pad_bounds)
        
        # --- THE OPACITY FIX ---
        # Force the faces to render as physical surfaces
        self.rep.GetFaceProperty().SetRepresentationToSurface()
        self.rep.GetFaceProperty().SetOpacity(0.3)
        self.rep.GetFaceProperty().SetColor(0.5, 0.5, 0.5)
        
        # Force the handles (spheres) to be highly visible as a scaling backup
        self.rep.GetHandleProperty().SetRepresentationToSurface()
        self.rep.GetHandleProperty().SetOpacity(1.0)
        self.rep.GetHandleProperty().SetColor(1.0, 0.0, 0.0) 
        
        self.widget.SetRepresentation(self.rep)
        self.widget.SetInteractor(plotter.iren.interactor)
        self.widget.AddObserver("InteractionEvent", self.on_interact)
        
    def on_interact(self, widget, event):
        self.rep.GetTransform(self.transform)
        for actor in self.actors:
            actor.SetUserTransform(self.transform)
            
    def set_active(self, active):
        self.widget.SetEnabled(active)

magnets = [
    MagnetGroup(create_magnet_actors((0, 0, 0))),
    MagnetGroup(create_magnet_actors((2, 0, 0)))
]

selected_group = {"current": None}

# -----------------------------
# Pure VTK Picker (Kills the Pink Box)
# -----------------------------
def on_left_click(obj, event):
    # Grab the raw 2D pixel coordinates of the click
    click_pos = plotter.iren.get_event_position()
    
    # Cast a ray into the scene to find the actor
    picker = vtk.vtkPropPicker()
    picker.PickProp(click_pos[0], click_pos[1], plotter.renderer)
    picked_actor = picker.GetViewProp()
    
    # Deselect if clicking empty space
    if picked_actor is None:
        if selected_group["current"] is not None:
            selected_group["current"].set_active(False)
            selected_group["current"] = None
        return

    # Select the correct magnet group
    for group in magnets:
        if picked_actor in group.actors:
            if selected_group["current"] == group:
                return 
            
            if selected_group["current"] is not None:
                selected_group["current"].set_active(False)
            
            selected_group["current"] = group
            group.set_active(True)
            return

# Bind the raw VTK click event instead of PyVista's wrapper
plotter.iren.add_observer("LeftButtonPressEvent", on_left_click)

plotter.add_axes()
plotter.show()

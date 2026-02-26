import vtk

class SceneObject:
    # ADDED allow_scaling=True parameter
    def __init__(self, plotter, actors, allow_scaling=True): 
        self.plotter = plotter
        self.actors = actors
        self.transform = vtk.vtkTransform()
        
        # 1. Compute bounds
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
        
        scale = 1.05
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
        
        # --- NEW LOGIC FOR SCALING ---
        if allow_scaling:
            self.rep.GetHandleProperty().SetColor(1, 0, 0)
            self.rep.GetHandleProperty().SetOpacity(1.0)
        else:
            # Disable uniform scaling (right-click drag)
            self.widget.SetScalingEnabled(0)
            # Disable individual axis scaling (grabbing the faces/handles)
            self.widget.SetMoveFacesEnabled(0)
            # Hide the spherical handles completely so they don't even appear
            self.rep.GetHandleProperty().SetOpacity(0.0)
        
        self.widget.SetRepresentation(self.rep)
        self.widget.SetInteractor(plotter.iren.interactor)
        self.widget.AddObserver("InteractionEvent", self.on_interact)
        
    def on_interact(self, widget, event):
        self.rep.GetTransform(self.transform)
        for actor in self.actors:
            actor.SetUserTransform(self.transform)
            
    def set_active(self, active):
        self.widget.SetEnabled(active)

    def get_rasterization_data(self):
        return {"transform": self.transform}
        
    def destroy(self):
        """Safely removes the widget and actors from the VTK scene."""
        self.widget.SetEnabled(0)
        self.widget.RemoveAllObservers()
        for actor in self.actors:
            self.plotter.remove_actor(actor)

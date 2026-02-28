import pyvista as pv
from scene_object import SceneObject

class ElectricCharge(SceneObject):
    def __init__(self, plotter, center=(0, 0, 0), charge=1.0):
        self._charge = charge
        
        # Create a simple, smooth sphere for the point charge
        self.mesh = pv.Sphere(radius=0.5, center=center)
        
        # Initial color based on the sign of the charge
        initial_color = "red" if self._charge > 0 else "blue"
        if self._charge == 0:
            initial_color = "grey"
            
        # Render it onto the PyVista plotter
        actor = plotter.add_mesh(self.mesh, color=initial_color, smooth_shading=True)
        
        # Pass to the parent class with scaling and rotation STRICTLY disabled
        super().__init__(plotter, [actor], allow_scaling=False, allow_rotation=False)

    @property
    def charge(self):
        return self._charge

    @charge.setter
    def charge(self, value):
        self._charge = value
        
        # Dynamically update the visual color if the slider crosses zero
        new_color = "red" if self._charge > 0 else "blue"
        if self._charge == 0:
            new_color = "grey"
            
        # Get the VTK float RGB tuple and apply it to the actor
        rgb = pv.Color(new_color).float_rgb
        self.actors[0].GetProperty().SetColor(rgb)
        self.plotter.render()

    def get_rasterization_data(self):
        # Override to append electrostatics data for the solver
        data = super().get_rasterization_data()
        data["type"] = "charge"
        data["charge"] = self._charge
        return data

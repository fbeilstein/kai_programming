import pyvista as pv
from scene_object import SceneObject

class Material(SceneObject):
    def __init__(self, plotter, center=(0, 0, 0), mu=1000.0):
        cube = pv.Cube(center=center, x_length=0.4, y_length=0.4, z_length=0.4)
        actor = plotter.add_mesh(cube, color="silver", opacity=0.3, show_edges=True)
        
        super().__init__(plotter, [actor])
        
        # Magnetic permeability
        self.mu = mu

    def get_rasterization_data(self):
        return {
            "type": "material",
            "transform": self.transform,
            "mu": self.mu
        }

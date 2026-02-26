import pyvista as pv
from scene_object import SceneObject

class Magnet(SceneObject):
    def __init__(self, plotter, center=(0, 0, 0), strength=1.0):
        # Two opaque cubes side-by-side
        cube_red = pv.Cube(center=(center[0]+0.2, center[1], center[2]), x_length=0.4, y_length=0.4, z_length=0.4)
        cube_blue = pv.Cube(center=(center[0]-0.2, center[1], center[2]), x_length=0.4, y_length=0.4, z_length=0.4)
        
        actor_red = plotter.add_mesh(cube_red, color="red", opacity=0.3, show_edges=True)
        actor_blue = plotter.add_mesh(cube_blue, color="blue", opacity=0.3, show_edges=True)
        
        super().__init__(plotter, [actor_red, actor_blue])
        
        # Physics properties
        self.strength = strength

    def get_rasterization_data(self):
        return {
            "type": "magnet",
            "transform": self.transform,
            "strength": self.strength
        }

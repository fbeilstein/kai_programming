import pyvista as pv
from scene_object import SceneObject

class Cannon(SceneObject):
    def __init__(self, plotter, center=(0, 0, 0)):
        base = pv.Sphere(center=(center[0]-0.2, center[1], center[2]), radius=0.15)
        barrel = pv.Cylinder(center=(center[0]+0.1, center[1], center[2]), direction=(1, 0, 0), radius=0.08, height=0.4)
        
        actor_base = plotter.add_mesh(base, color="black", opacity=0.5)
        actor_barrel = plotter.add_mesh(barrel, color="darkgreen", opacity=0.5)
        
        # FIX: Explicitly disable scaling for this object type
        super().__init__(plotter, [actor_base, actor_barrel], allow_scaling=False)
        
        self.is_firing = False
        self.particle_charge = -1.0
        self.velocity = 500.0

    def get_rasterization_data(self):
        return {
            "type": "cannon",
            "transform": self.transform,
            "is_firing": self.is_firing,
            "charge": self.particle_charge,
            "velocity": self.velocity
        }

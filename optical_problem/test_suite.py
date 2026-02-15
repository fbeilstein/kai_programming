import unittest
import numpy as np
# Import your functions from your implementation file
from implementation_tasks import (
    intersect_line_infinite, intersect_segment,
    intersect_circle, intersect_arc,
    calculate_normal_segment, calculate_normal_arc,
    refract_vector
)

class TestOpticsMath(unittest.TestCase):

    # --- LEVEL 1: Infinite Line Math ---
    def test_l1_infinite_hit(self):
        """Checks if the math finds the intersection point."""
        origin, rd = np.array([0.0, 0.0]), np.array([1.0, 0.0])
        p1, p2 = np.array([5.0, -5.0]), np.array([5.0, 5.0])
        result = intersect_line_infinite(origin, rd, p1, p2)
        
        # We now check for a coordinate point [5, 0], not a scalar
        self.assertIsInstance(result, np.ndarray)
        np.testing.assert_allclose(result, [5.0, 0.0], atol=1e-7)

    def test_l1_parallel(self):
        """Parallel lines should return None."""
        origin, rd = np.array([0.0, 0.0]), np.array([1.0, 0.0])
        p1, p2 = np.array([0.0, 5.0]), np.array([10.0, 5.0])
        result = intersect_line_infinite(origin, rd, p1, p2)
        self.assertIsNone(result)

    # --- LEVEL 2: Segment Logic ---
    def test_l2_segment_bounds(self):
        """Check that a point is returned only if within the [0, 1] segment."""
        # HIT: Segment crosses the ray path
        p_hit = intersect_segment(np.array([0,0]), np.array([1,0]), np.array([5,-1]), np.array([5,1]))
        self.assertIsInstance(p_hit, np.ndarray)
        np.testing.assert_allclose(p_hit, [5.0, 0.0], atol=1e-7)
        
        # MISS: Segment is out of bounds
        p_miss = intersect_segment(np.array([0,0]), np.array([1,0]), np.array([5,2]), np.array([5,4]))
        self.assertIsNone(p_miss)

    # --- LEVEL 3: Circle Math ---
    def test_l3_circle_solutions(self):
        """Should return a list of coordinate points."""
        # Ray starts at x=-10, moves right. Hits circle radius 5 at center [0,0]
        pts = intersect_circle(np.array([-10.0, 0.0]), np.array([1.0, 0.0]), np.array([0.0, 0.0]), 5.0)
        
        # Expect a list/tuple of two points: [-5, 0] and [5, 0]
        self.assertEqual(len(pts), 2)
        x_coords = sorted([p[0] for p in pts])
        self.assertAlmostEqual(x_coords[0], -5.0)
        self.assertAlmostEqual(x_coords[1], 5.0)

    # --- LEVEL 4: Arc Sector Logic ---
    def test_l4_arc_sector(self):
        """Check if angular sector logic correctly returns the valid hit point."""
        center = np.array([0.0, 0.0])
        axis = np.array([-1.0, 0.0]) # Arc faces left
        cos_half = np.cos(np.radians(45))
        
        # Ray from left hits [-5, 0] and [5, 0]. Only [-5, 0] is inside the left-facing arc.
        p = intersect_arc(np.array([-10.0, 0.0]), np.array([1.0, 0.0]), center, 5.0, axis, cos_half)
        self.assertIsInstance(p, np.ndarray)
        np.testing.assert_allclose(p, [-5.0, 0.0], atol=1e-7)

    # --- LEVEL 5 & 6: Normals ---
    def test_l5_segment_normal_flip(self):
        """Verify normal points against the ray direction."""
        rd = np.array([1.0, 0.0])
        p1, p2 = np.array([5.0, -1.0]), np.array([5.0, 1.0])
        n = calculate_normal_segment(rd, p1, p2)
        self.assertLess(np.dot(rd, n), 0)
        self.assertAlmostEqual(n[0], -1.0)

    def test_l6_arc_normal_radial(self):
        """Verify normal is radial and flipped to face ray."""
        center = np.array([0.0, 0.0])
        hit = np.array([-5.0, 0.0])
        rd = np.array([1.0, 0.0])
        n = calculate_normal_arc(hit, rd, center)
        self.assertAlmostEqual(n[0], -1.0)

    # --- LEVEL 7: Refraction ---
    def test_l7_refraction(self):
        """Verifies Snell's Law vector output."""
        ray_dir = np.array([1.0, -1.0]) / np.sqrt(2) # 45 degrees
        normal = np.array([0.0, 1.0])
        new_dir = refract_vector(ray_dir, normal, 1.0, 1.5)
        self.assertIsNotNone(new_dir)
        self.assertAlmostEqual(np.linalg.norm(new_dir), 1.0)
        self.assertLess(abs(new_dir[0]), abs(ray_dir[0])) # Bends toward normal

if __name__ == '__main__':
    unittest.main()

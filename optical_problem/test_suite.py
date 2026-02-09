import unittest
import numpy as np
# Import your functions from your implementation file
from implementation_tasks import (
    intersect_line_infinite, intersect_segment,
    intersect_circle_infinite, intersect_arc,
    calculate_normal_segment, calculate_normal_arc
)

class TestOpticsMath(unittest.TestCase):

    # --- LEVEL 1: Infinite Line Math ---
    def test_l1_infinite_hit(self):
        """Checks if the math finds the intersection of two infinite lines."""
        origin, rd = np.array([0, 0]), np.array([1, 0])
        p1, p2 = np.array([5, -5]), np.array([5, 5])
        t, u = intersect_line_infinite(origin, rd, p1, p2)
        self.assertAlmostEqual(t, 5.0)
        self.assertAlmostEqual(u, 0.5)

    def test_l1_parallel(self):
        """Parallel lines should return infinity to avoid division by zero."""
        origin, rd = np.array([0, 0]), np.array([1, 0])
        p1, p2 = np.array([0, 5]), np.array([10, 5])
        t, u = intersect_line_infinite(origin, rd, p1, p2)
        self.assertEqual(t, float('inf'))

    # --- LEVEL 2: Segment Logic ---
    def test_l2_segment_bounds(self):
        """Check that t is returned only if the hit is within the [0, 1] segment."""
        # Horizontal ray at y=0. Segment goes from x=5, y=-1 to y=1. (HIT)
        t_hit = intersect_segment(np.array([0,0]), np.array([1,0]), np.array([5,-1]), np.array([5,1]))
        self.assertAlmostEqual(t_hit, 5.0)
        # Same ray, but segment is shifted to y=2 to y=4. (MISS)
        t_miss = intersect_segment(np.array([0,0]), np.array([1,0]), np.array([5,2]), np.array([5,4]))
        self.assertEqual(t_miss, float('inf'))

    # --- LEVEL 3: Circle Math ---
    def test_l3_circle_solutions(self):
        """Quadratic should return two points for a ray passing through a circle."""
        ts = intersect_circle_infinite(np.array([-10, 0]), np.array([1, 0]), np.array([0, 0]), 5.0)
        self.assertEqual(len(ts), 2)
        self.assertTrue(5.0 in ts and 15.0 in ts)

    # --- LEVEL 4: Arc Sector Logic ---
    def test_l4_arc_sector(self):
        """Check if the angular sector logic correctly ignores the back of the circle."""
        center = np.array([0, 0])
        axis = np.array([-1, 0]) # Arc faces left
        cos_half = np.cos(np.radians(45))
        # Ray from left hits front of circle at t=5, back at t=15.
        # Only t=5 is inside the 'left-facing' arc.
        t = intersect_arc(np.array([-10, 0]), np.array([1, 0]), center, 5.0, axis, cos_half)
        self.assertAlmostEqual(t, 5.0)

    # --- LEVEL 5: Segment Normal ---
    def test_l5_segment_normal_flip(self):
        """The normal must ALWAYS point against the incoming ray."""
        rd = np.array([1.0, 0.0]) # Ray moving right
        p1, p2 = np.array([5, -1]), np.array([5, 1])
        n = calculate_normal_segment(rd, p1, p2)
        # Normal should be [-1, 0]. Dot product must be negative.
        self.assertLess(np.dot(rd, n), 0)
        self.assertAlmostEqual(n[0], -1.0)

    # --- LEVEL 6: Arc Normal ---
    def test_l6_arc_normal_radial(self):
        """Normal should be the radial vector (Hit - Center) flipped to face the ray."""
        center = np.array([0, 0])
        hit = np.array([-5, 0])
        rd = np.array([1, 0]) # Ray moving right hits left side of circle
        n = calculate_normal_arc(hit, rd, center)
        # Vector (Hit-Center) is [-5, 0]. rd is [1, 0]. Dot is -5 (already faces ray).
        self.assertAlmostEqual(n[0], -1.0)

if __name__ == '__main__':
    unittest.main()

#!/usr/bin/env python
print("In module products __package__, __name__ ==", __package__, __name__)
from abc import ABC, abstractmethod
import math

class Region(ABC):
    @abstractmethod
    def contains(self, point):
        """
        Check if the given point (x, y, z) is inside the region.
        
        Parameters:
        - point: A tuple representing a point in space (x, y, z).
        
        Returns:
        - bool: True if the point is inside the region, False otherwise.
        """
        pass

class Sphere(Region):
    def __init__(self, center, radius):
        """
        Initialize a Sphere instance.
        
        Parameters:
        - center: A tuple representing the center of the sphere (x, y, z).
        - radius: The radius of the sphere.
        """
        self.center = center
        self.radius = radius

    def contains(self, point):
        """
        Check if the given point is inside the sphere.
        """
        dx = point[0] - self.center[0]
        dy = point[1] - self.center[1]
        dz = point[2] - self.center[2]
        distance_squared = dx**2 + dy**2 + dz**2
        return distance_squared <= self.radius**2

class Block(Region):
    def __init__(self, corner):
        """
        Initialize a Block instance with one corner at the origin (0, 0, 0)
        and the opposite corner in the first quadrant.
        
        Parameters:
        - corner: A tuple representing the opposite corner of the block (x, y, z)
                  from the origin, ensuring x, y, z are positive.
        """
        assert all(c > 0 for c in corner), "Corner must be in the first quadrant"
        self.corner = corner

    def contains(self, point):
        """
        Check if the given point is inside the block defined by the origin
        and the corner in the first quadrant.
        """
        return all(0 < point[i] < self.corner[i] for i in range(3))

import os
import sys

TEST_MODE = (v:=os.environ.get('SIM_TEST_MODE','').upper()) \
    and v[0] in ('T', 'Y', '1')

if __name__ == '__main__':
    sys.exit(0) ### avoid running tests in production

if TEST_MODE:
    import pytest
    
    def test_block_contains_point_inside():
        block = Block((10, 10, 10))
        assert block.contains((5, 5, 5)) == True, \
        "The point is inside the block but was not recognized as such."
    
    def test_block_contains_point_outside():
        block = Block((10, 10, 10))
        assert block.contains((11, 11, 11)) == False, \
        "The point is outside the block but was incorrectly recognized as inside."
    
    def test_block_contains_point_on_edge():
        block = Block((10, 10, 10))
        assert block.contains((10, 10, 10)) == False, \
        "The point is on the edge of the block but was recognized as inside."
    
    def test_block_contains_point_at_origin():
        block = Block((10, 10, 10))
        assert block.contains((0, 0, 0)) == False, \
        "The origin point is outside the block but was not recognized as such."
    
    def test_block_corner_in_first_quadrant():
        with pytest.raises(AssertionError):
            Block((-1, -1, -1))  # This should raise an AssertionError because the corner is not in the first quadrant
    
    # Optional: You might want to test with floating point numbers as well
    def test_block_contains_point_with_floats():
        block = Block((10.5, 10.5, 10.5))
        assert block.contains((10.1, 10.1, 10.1)) == True, \
        "The point is inside the block but was not recognized as such with floating point coordinates."

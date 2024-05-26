#!/usr/bin/env python

import random


from abc import ABC, abstractmethod
from typing import List, Tuple

from montecarlosim.container import Region, Block

class Initializer(ABC):
    def __init__(self, container):
        """
        Initialize the Initializer with a given container.

        Parameters:
        - container: An instance of a class that inherits from Region.
        """
        self.container = container

    @abstractmethod
    def initialize(self, N: int) -> List[Tuple[float, float, float]]:
        raise NotImplementedError

def get_initializer(
    initializer_type: str
    , container: Region
    , *args, **kwargs
) -> Initializer:
    """
    Create an instance of an Initializer subclass based on the given type name.

    Parameters:
    - initializer_type: The name of the Initializer subclass to instantiate ('random').
    - container: An instance of a class that inherits from Region.

    Returns:
    - An instance of the specified Initializer subclass.

    Raises:
    - ValueError: If the specified initializer_type is not recognized.
    """
    initializer_classes = {
        'random': RandomInitializer
    }

    if initializer_type not in initializer_classes:
        raise ValueError(f"Unknown initializer type: {initializer_type}")

    return initializer_classes[initializer_type](container, *args, **kwargs)

class RandomInitializer(Initializer):
    def initialize(self, N: int) -> List[Tuple[float, float, float]]:
        """
        Generate a set of N points contained within the container using random sampling.

        Parameters:
        - N: The number of points to generate.

        Returns:
        - List[Tuple[float, float, float]]: A list of points contained within the container.
        """
        points = []
        while len(points) < N:
            point = (
                random.uniform(0, self.container.corner[0]),
                random.uniform(0, self.container.corner[1]),
                random.uniform(0, self.container.corner[2])
            )
            if self.container.contains(point):
                points.append(point)

        return points

#############################################################################
########################## TESTING CODE BELOW  ##############################
#############################################################################


import os
import sys

TEST_MODE = (v:=os.environ.get('SIM_TEST_MODE','').upper()) \
    and v[0] in ('T', 'Y', '1')

if __name__ == '__main__':
    sys.exit(0) ### avoid running tests in production

if TEST_MODE:
    import pytest

    def test_get_initializer_random():
        block = Block((10, 10, 10))
        initializer = get_initializer('random', block)
        assert isinstance(initializer, RandomInitializer), \
        "The created object is not an instance of RandomInitializer"
        
        NUM_OF_POINTS = 100000 ### Generating more freezes the pc

        points = initializer.initialize( NUM_OF_POINTS )
        assert len(points) == NUM_OF_POINTS, "The number of generated points is incorrect"
        assert all(block.contains(point) for point in points), \
        "Not all points are contained within the block"

    ### TODO: Test with a sphere. This will reveal a bug

    def test_get_initializer_invalid_type():
        block = Block((10, 10, 10))
        with pytest.raises(ValueError, match="Unknown initializer type: invalid"):
            get_initializer('invalid', block)
            

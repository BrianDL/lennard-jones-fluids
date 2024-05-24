#!/usr/bin/env python

from abc import ABC, abstractmethod
import random
from typing import List, Tuple

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
        """
        Generate a set of N points contained within the container.

        Parameters:
        - N: The number of points to generate.

        Returns:
        - List[Tuple[float, float, float]]: A list of points contained within the container.
        """
        pass

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

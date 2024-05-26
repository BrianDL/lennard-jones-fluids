#!/usr/bin/env python

import math
import numpy as np

from montecarlosim.container import Region, get_region
from montecarlosim.initializer import Initializer

def lj_potential(r) -> float: ### verify this is correct
    return 4 * (math.pow(r, -12) - math.pow(r, -6))


class Simulation():
    def __init__(
        self 
        , container:Region
        , num_of_particles:int = 0
        , beta = 1 ### how to select this?
        , step_size:float = 0.01
        , dynamic_factor:float = 0.0 
        , stop_condition:str = 'max_steps_10K' ### Change this to something better
        , initializer:Initializer = None ### Change this to something better
    ):

        self.step_size = step_size
        self.beta = beta
        self.stop_condition = stop_condition
        self.dynamic_factor = dynamic_factor
        self.initializer = initializer

        self.container = container

        assert num_of_particles > 0, \
            "num_of_particles must be a positive number"

        self.num_of_particles = num_of_particles
        
        self.steps = []
        self.energies = []
        self.system = None

    def energy(self, system)->float:
        total_energy = 0
        for i in range(self.num_of_particles):
            pi = (system[0][i], system[1][i], system[2][i])
            if not self.container.contains(pi):
                return np.inf

            for j in range(i +1, self.num_of_particles):
                pj = (system[0][j], system[1][j], system[2][j])
                if not self.container.contains(pj):
                    return np.inf

                r = np.sqrt(
                    (pi[0] - pj[0])**2 +
                    (pi[1] - pj[1])**2 +
                    (pi[2] - pj[2])**2
                )

                if r < 1e-6: return np.inf

                total_energy += lj_potential(r)
        
        return total_energy
    
    def initialize(self):
        if self.initializer:
            self.system = \
                self.initializer.initialize(
                    self.num_of_particles
                )
            
            return True
        
        max_x, max_y, max_z = self.container.corner

        xs = max_x * np.random.rand(self.num_of_particles)
        ys = max_y * np.random.rand(self.num_of_particles)
        zs = max_z * np.random.rand(self.num_of_particles)

        for i in range(self.num_of_particles):
            pi = (xs[i], ys[i], zs[i])

            while not self.container.contains(pi):
                pi = (
                    max_x * np.random.rand(),
                    max_y * np.random.rand(),
                    max_z * np.random.rand()
                )
            
            xs[i], ys[i], zs[i] = pi[0], pi[1], pi[2]
        
        system = (xs, ys, zs)
        total_energy = self.energy( system )

        if total_energy == np.inf:
            return self.initialize()

        self.system = system
        self.energies.append( total_energy )
        self.steps.append( system )

        return True


    def is_finished(self) -> bool:
        return any([
            self.stop_condition =='max_steps_10' and len(self.steps) >= 10
            , self.stop_condition =='max_steps_100' and len(self.steps) >= 100
            , self.stop_condition =='max_steps_1K' and len(self.steps) >= 1000
            , self.stop_condition =='max_steps_10K' and len(self.steps) >= 10000
            , self.stop_condition =='max_steps_100K' and len(self.steps) >= 100000
        ])

    def start(self):
        while not self.is_finished(): self.move()

    def move(self):
        
        if not self.system:
            self.initialize()

        ### Calculating next step
        xs = self.system[0].copy()
        ys = self.system[1].copy()
        zs = self.system[2].copy()
        
        ### Randomly move a sample of particles
        num_of_particles_to_move = \
            int(self.num_of_particles * self.dynamic_factor) \
            if self.dynamic_factor else 1

        for j in range(num_of_particles_to_move):
            ### Choose a random particle
            i = np.random.randint(self.num_of_particles)
            xs[i] += np.random.uniform(-self.step_size, self.step_size)
            ys[i] += np.random.uniform(-self.step_size, self.step_size)
            zs[i] += np.random.uniform(-self.step_size, self.step_size)

        new_step = (xs, ys, zs)
        new_step_energy = self.energy(new_step)

        de = new_step_energy - self.energies[-1]

        is_accepted = de < 0 \
            or np.exp(-self.beta * de) * np.random.rand()

        if is_accepted:
            self.system = new_step
            self.steps.append(new_step)
            self.energies.append(new_step_energy)
        
        return is_accepted
        

#########################################################################
###################### TESTING CODE BELOW ###############################
#########################################################################

import os
import sys

TEST_MODE = (v:=os.environ.get('SIM_TEST_MODE','').upper()) \
    and v[0] in ('T', 'Y', '1')

if __name__ == '__main__':
    if TEST_MODE: sys.exit(0)

if TEST_MODE:
    import pytest

    def test_initialize_without_initializer():
        container = get_region('block', (10, 10, 10))
        sim = Simulation(container, 1000)
        sim.initialize()

        ### Verify that the system is initialized
        assert sim.system is not None, \
            "The system is not initialized"

        ### Verify the number of particles
        assert len(sim.system[0]) == sim.num_of_particles, \
            "The number of particles is incorrect"

        ### Verify that all particles are contained within the container
        for i in range(sim.num_of_particles):
            pi = (sim.system[0][i], sim.system[1][i], sim.system[2][i])
            assert container.contains(pi), \
                f"Particle {i} is not contained within the container"
    
    ### test move function
    def test_move():
        container = get_region('block', (10, 10, 10))
        sim = Simulation(container, 100)
        move_result = sim.move()

        if not move_result:
            assert len(sim.steps) == 1, "There should be one step"
            assert len(sim.energies) == 1, "There should be one energy"

        else:
            assert len(sim.steps) == len(sim.energies), \
                "There should be an equal number of steps and energies"
            
            assert len(sim.steps) == 2, \
                "There should be two steps"
            
            xs1, ys1, zs1 = sim.steps[0][0], sim.steps[0][1], sim.steps[0][2]
            xs2, ys2, zs2 = sim.steps[1][0], sim.steps[1][1], sim.steps[1][2]
            steps_are_equal = len(xs1) == len(xs2) \
                and len(ys1) == len(ys2) \
                and len(zs1) == len(zs2) \
                and all([ xs1[i] == xs2[i] for i in range(len(xs1))]) \
                and all([ ys1[i] == ys2[i] for i in range(len(ys1))]) \
                and all([ zs1[i] == zs2[i] for i in range(len(zs1))])

            assert not steps_are_equal, "Steps should be different"
            
    def test_energy_decreases_along_simulation():
        container = get_region('block', (10, 10, 10))
        sim = Simulation(container, 5000
            , stop_condition='max_steps_10')
        sim.start()
    
        # Verify that the simulation has exactly 10 steps
        assert len(sim.steps) == 10, "The simulation should have exactly 10 steps"
    
        # Verify that the energy decreases from the first step to the last step
        initial_energy = sim.energies[0]
        final_energy = sim.energies[-1]
        

        print(sim.energies)
        assert final_energy < initial_energy, \
            "The energy should decrease from the first step to the last step"
            

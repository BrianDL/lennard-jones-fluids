#!/usr/bin/env python3

import time

from montecarlosim.simulation import Simulation
from montecarlosim.container import Region, get_region
from montecarlosim.plotter import Plotter

def main():

    L = 0.01
    N = 500

    print(f"Simulation started at {time.asctime()}...")
    t0 = time.time()

    container = get_region('block', (L, L, L))
    sim = Simulation(container, N
        , step_size=0.1
        , epsilon=1
        , beta=1
        , stop_condition='max_steps_5K'
        )
    sim.start()

    dt = time.time() - t0
    print(f"Simulation took {dt} ms")

    plotter = Plotter(sim)
    plotter.plot_energy()
    plotter.plot_xz_projection(0)
    plotter.plot_xz_projection()

if __name__ == '__main__':
    main()

#!/usr/bin/env python3

import time

from montecarlosim.simulation import Simulation
from montecarlosim.container import Region, get_region
from montecarlosim.plotter import Plotter

def main():

    L = 500
    N = 1200

    print(f"Simulation started at {time.asctime()}...")
    t0 = time.time()

    container = get_region('sphere', radius=L)
    sim = Simulation(container, N
        , step_size=0.2
        , stop_condition='max_steps_6K'
        , init_pressure=0.2
        )
    sim.start()

    dt = time.time() - t0
    print(f"Simulation took {dt} s")

    plotter = Plotter(sim)
    plotter.plot_energy()
    
    plotter.plot_xz_projection(0)
    plotter.plot_xz_projection()
    
    plotter.plot_3d_state()

if __name__ == '__main__':
    main()

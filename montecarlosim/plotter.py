import matplotlib.pyplot as plt

class Plotter:
    def __init__(self, simulation):
        self.simulation = simulation

    def plot_energy(self):
        energies = self.simulation.energies
        steps = range(len(energies))
        
        plt.figure(figsize=(10, 6))
        plt.plot(steps, energies, marker='.', linestyle='-', color='b', markersize=2)
        plt.title('Energy Variation Over Time')
        plt.xlabel('Step')
        plt.ylabel('Energy')
        plt.grid(True)
        plt.savefig('energy_over_time.png')
        plt.close()
        

    def plot_xz_projection(self, step_index=-1):
        """
        Plots the projection of the system on the xz plane for a given snapshot,
        fitting the plot to the whole capacity of the container.
        If step_index is not provided, it plots the last snapshot by default.
        """
        system_snapshot = self.simulation.steps[step_index]
        xs, _, zs = system_snapshot  # Extract xs and zs from the snapshot
        max_x, _, max_z = self.simulation.container.corner  # Get container dimensions

        plt.figure(figsize=(10, 6))
        plt.scatter(xs, zs, color='r')
        plt.xlim(0, max_x)
        plt.ylim(0, max_z)
        plt.title('Projection of the System on the XZ Plane')
        plt.xlabel('X')
        plt.ylabel('Z')
        plt.grid(True)
        plt.savefig(f'xz_projection_step_{step_index}.png')
        plt.close()
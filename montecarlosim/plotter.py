import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # Import the 3D plotting toolkit

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
        system_snapshot = self.simulation.steps[step_index]
        xs, _, zs = system_snapshot
        max_x, _, max_z = self.simulation.container.corner

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

    def plot_3d_state(self, step_index=-1):
        """
        Plots the 3D state of the system for a given snapshot.
        If step_index is not provided, it plots the last snapshot by default.
        """
        system_snapshot = self.simulation.steps[step_index]
        xs, ys, zs = system_snapshot  # Extract xs, ys, and zs from the snapshot

        fig = plt.figure(figsize=(10, 8))
        ax = fig.add_subplot(111, projection='3d')  # Create a 3D subplot
        ax.scatter(xs, ys, zs, color='r', marker='o')  # Plot the points in 3D space

        ax.set_title('3D State of the System')
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        ax.grid(True)
        plt.savefig(f'3d_state_step_{step_index}.png')
        plt.close()
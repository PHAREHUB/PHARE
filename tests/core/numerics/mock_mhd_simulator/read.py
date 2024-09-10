import h5py
import matplotlib.pyplot as plt
import numpy as np

file_path = "hall_harris3.h5"


def read_and_plot_h5(file_path, timestep, quantity):
    """
    Reads 2D data from an HDF5 file and plots it.

    Parameters:
        file_path (str): Path to the HDF5 file.
        timestep (str): Timestep to read, e.g., "0.499800".
        quantity (str): Quantity to read, e.g., "p".
    """
    with h5py.File(file_path, "r") as h5_file:
        # Check if the specified timestep exists
        if timestep not in h5_file:
            raise ValueError(f"Timestep '{timestep}' not found in the file.")

        # Check if the specified quantity exists
        if quantity not in h5_file[timestep]:
            raise ValueError(
                f"Quantity '{quantity}' not found in timestep '{timestep}'."
            )

        # Read the data
        data = h5_file[f"{timestep}/{quantity}"][:]

        # Plot the data
        plt.figure(figsize=(8, 6))
        plt.imshow(data.T, origin="lower", cmap="coolwarm", aspect="auto")
        plt.colorbar(label=quantity)
        plt.title(f"{quantity} at {timestep}")
        plt.xlabel("X-axis")
        plt.ylabel("Y-axis")
        plt.show()


timestep_to_plot = "0.002000"
quantity_to_plot = "vz"

read_and_plot_h5(file_path, timestep_to_plot, quantity_to_plot)

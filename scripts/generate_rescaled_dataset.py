import os
import numpy as np
from jsonargparse import CLI
from typing import Optional


def get_ipi_potential(positions, r0, D, a):
    """
    positions: np.array of shape (n_beads, 3)
    """
    if positions.ndim == 3:
        distances = np.linalg.norm(positions, axis=2)
    elif positions.ndim == 2:
        distances = np.linalg.norm(positions, axis=1)
    adjusted_distances = distances - r0
    potentials = D * (
        np.exp(-2 * a * adjusted_distances) - 2 * np.exp(-a * adjusted_distances)
    )
    return potentials


def generate_rescaled_dataset(
    P: int,
    T_0: int,
    T_1: int,
    source_location: str,
    timestep: float,
    base_stride: int,
    additional_stride: int,
    equilibration_time: float,
    total_time: float,
    bead_stride: int,
    dataset_location: str,
    output_prefix: str,
    input_folder_name: Optional[str] = None,
    traj_file_root: Optional[str] = None,
    mass_0: Optional[float] = 1728.3071,
    mass_1: Optional[float] = 1728.3071,
):

    """
    Generate dataset for Morse potential

    Parameters:
    ----------

    P: int
        Number of beads
    T: int
        Temperature. Should be used to name the simulation folder
    timestep: float
        Trajectory timestep, fs
    base_stride: int
        Stride used to generate the trajectory, steps
    additional_stride: int
        Stride applied to the data during the processing
    equilibration_time: float
        Part of the trajectory corresponding to this time in ps will be removed
    total_time: float
        Total time used for dataset, ps
    bead_stride: int
        Stride applied to  bead index
    dataset_location: str
        Location to put the final dataset
    output_prefix: str
        The output prefix for the output files
    input_folder_name: str
        Name of the input folder within the source location. If none,
        T is used as input_folder_prefix
    traj_file_root: Optional[str]
        Root of the trajectory names. The total name is assumeed to be {root}.pos_{bead_ndx}.
        If None, root = {T}_simulation.
    mass: Optional[float] = 1728.3071
    """
    D = 0.18748511263179304
    a = 1.1562696428501682
    r0 = 1.8323926

    kBT_0 = T_0 * 3.1668105e-06
    kBT_1 = T_1 * 3.1668105e-06

    root_directory = source_location

    if input_folder_name is None:
        input_folder_name = str(T)

    if traj_file_root is None:
        traj_file_root = f"{T}_simulation"
    pos = []
    cwd = os.getcwd()
    os.chdir(f"{root_directory}/{input_folder_name}")

    # Initial frame
    frame_min = int(equilibration_time * 1000 / (timestep * base_stride))
    frame_max = int(total_time * 1000 / (timestep * base_stride))

    for iP in range(P):

        iP = str(iP).zfill(2)
        print("Read pos and for for bead: ", iP)

        os.system(
            f"grep H {traj_file_root}.pos_"
            + iP
            + ".xyz"
            + " | "
            + "awk '{print $2, $3, $4}' | awk '!seen[$0]++' >"
            + f"{traj_file_root}.pos_"
            + iP
            + ".data"
        )
        pos.append(
            np.loadtxt(f"{traj_file_root}.pos_" + iP + ".data")[
                frame_min:frame_max:additional_stride
            ]
        )
        os.system(f"rm {traj_file_root}.pos_" + iP + ".data")

    os.chdir(cwd)
    pos = np.array(pos)
    print(pos.shape)

    # Get potential energy for original coordinates
    potential = get_ipi_potential(pos, r0, D, a)
    print(f"Potential shape: {potential.shape}")
    mean_potential = np.sum(potential, axis=0) / P
    print(f"Mean potential shape: {mean_potential.shape}")

    # Get centroid positions
    centroid_pos = np.mean(pos, axis=0)

    # Get rescaled coordinates
    pos_rescaled = centroid_pos + np.sqrt(mass_0 * T_0 / (mass_1 * T_1)) * (
        pos - centroid_pos
    )

    potential_rescaled = get_ipi_potential(pos_rescaled, r0, D, a)
    assert potential_rescaled.shape[0] == P
    mean_potential_rescaled = np.sum(potential_rescaled, axis=0) / P

    # Get Yamamoto weights
    weights = np.exp(-mean_potential_rescaled / kBT_1 + mean_potential / kBT_0)
    weights = weights / np.sum(weights)
    weights = weights.flatten()

    # Get selector
    selection = np.random.choice(
        list(range(len(weights))), size=len(weights), p=weights
    )

    # Calculating spring forces
    hbar = 1

    prefactor = -mass_1 * P * (kBT_1) ** 2 / (hbar) ** 2
    spring_forces = []

    for iP in range(P):
        iPp1 = (iP + 1) % 64
        iPm1 = iP - 1

        spring = prefactor * (
            2 * pos_rescaled[iP] - pos_rescaled[iPm1] - pos_rescaled[iPp1]
        )
        spring_forces.append(spring)

    total_positions = np.array(
        [pos_rescaled[i][selection] for i in range(0, P, bead_stride)]
    )
    spring_forces = np.array(
        [spring_forces[i][selection] for i in range(0, P, bead_stride)]
    )

    spring_forces_data = spring_forces.reshape(1, -1, 3)[0]
    total_coordinates_data = total_positions.reshape(-1, 3)

    np.save(
        f"{dataset_location}/{output_prefix}_total_coordinates.npy",
        total_coordinates_data,
    )
    np.save(f"{dataset_location}/{output_prefix}_spring_forces.npy", spring_forces_data)
    return


if __name__ == "__main__":
    CLI(generate_rescaled_dataset)

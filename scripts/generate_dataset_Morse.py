import os
import numpy as np
from jsonargparse import CLI


def generate_dataset(
    P: int,
    T: int,
    source_location: str,
    timestep: float,
    base_stride: int,
    additional_stride: int,
    equilibration_time: float,
    total_time: float,
    bead_stride: int,
    dataset_location: str,
    output_prefix: str,
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
    """
    root_directory = source_location

    positions = {}
    cwd = os.getcwd()
    os.chdir(f"{root_directory}/{T}")

    # Initial frame
    frame_min = int(equilibration_time * 1000 / (timestep * base_stride))
    frame_max = int(total_time * 1000 / (timestep * base_stride))

    for iP in range(P):

        iP = str(iP).zfill(2)
        print("Read pos and for for bead: ", iP)

        os.system(
            f"grep H {T}_simulation.pos_"
            + iP
            + ".xyz"
            + " | "
            + "awk '{print $2, $3, $4}' | awk '!seen[$0]++' >"
            + f"{T}_simulation.pos_"
            + iP
            + ".data"
        )
        positions[iP] = np.loadtxt(f"{T}_simulation.pos_" + iP + ".data")[
            frame_min:frame_max:additional_stride
        ]
        os.system(f"rm {T}_simulation.pos_" + iP + ".data")

    os.chdir(cwd)

    # Calculating spring forces
    kBT = T * 3.1668105e-06
    hbar = 1
    m = 1.72830710e03

    prefactor = -m * P * (kBT) ** 2 / (hbar) ** 2
    spring_forces = []

    for iP in range(P):
        iPp1 = str((iP + 1) % P).zfill(2)
        iPm1 = str((iP - 1) % P).zfill(2)
        iP = str((iP) % P).zfill(2)

        spring = prefactor * (2 * positions[iP] - positions[iPm1] - positions[iPp1])
        spring_forces.append(spring)

    total_positions = np.array(
        [positions[str(i).zfill(2)] for i in range(0, P, bead_stride)]
    )
    spring_forces = np.array([spring_forces[i] for i in range(0, P, bead_stride)])

    spring_forces_data = spring_forces.reshape(1, -1, 3)[0]
    total_coordinates_data = total_positions.reshape(-1, 3)

    np.save(
        f"{dataset_location}/{output_prefix}_total_coordinates.npy",
        total_coordinates_data,
    )
    np.save(f"{dataset_location}/{output_prefix}_spring_forces.npy", spring_forces_data)
    return


if __name__ == "__main__":
    CLI(generate_dataset)

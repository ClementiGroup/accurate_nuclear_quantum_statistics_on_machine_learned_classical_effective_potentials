"""
Takes raw i-pi output and converts it into
an hdf5 dataset. 
The script assumes the following structure of the input data:
    - each folder represents one temperature
    - within each temperature folder, there is a separate file
      with positions, classical forces and spring forces for each 
      bead. If needed, the average potential energy and total spring
      energy can be added to the dataset. 

The structure of the dataset is the following:
    - There is one group per folder. Naming conventions
     are '{system}_T={temperature}K"
     - within these first-level groups, there is a single group per bead
       the naming conventions for bead groups 'bead_{i}
    - withing these bead groups, there are datasets for positions, forces, and spring energy

It is assumed, that the i-pi output can be read by ase.
However, ase doesn't recognize D and T isotopes of H. Sould output files contain D or T,
they should be preprocessed such that nonstandard isotopes are replaced with H.

TODO: Add handling of different isotopes of H

"""
import h5py
import numpy as np
import logging
from ase.io import read
from ase.units import Bohr, Hartree, kJ, second, kB
from jsonargparse import CLI
from typing import Dict, List, Optional, Union
from jsonargparse.typing import PositiveInt, Path_fc, Path_drw, Path_fr
from utils import compute_spring_energies
from scipy.constants import hbar as hbar_CI
from ase.data import atomic_masses_legacy

from log_utils import log_execution_info

logger = logging.getLogger(__name__)


def get_number_of_frames(
    P: int,
    order_of_beads: int,
    frame_selection: str,
    run_location: Path_drw,
    data_dict: Dict,
    first_frame: int = 0,
    last_frame: int = -1,
    stride: int = 1,
):
    """
    Get number of frames to use in the dataset.
    If last_frame is not specified, use the minimum number of frames
    over all beads, all datasets
    """
    if last_frame < 0:
        n_frame_list = []
        for bead_ndx in range(P):
            # First, need to figure number of frames in each dataset.
            # This is only required if number of frames is unknown a priori (i.e., `last_frame` < 0)
            # The reason we have to do it is that sometimes different beads have slightly different
            # number of frames, so we need to take the minimum number of frames over all beads.

            # Will have to go over data twice: first,  to get number of frames
            # in each dataset, second, use the minimum number of frames to record to the
            # dataset. When different datasets for the same bead have different number of frames, a warning is
            # given. Also, for the purpose of force optimization, different replicas all should have the same number of frames
            logger.info("Unspecified number of frames")
            logger.info("Determine the common number of frames to use")
            for properties in data_dict.values():
                file_template = properties[0]
                full_path = (
                    f"{run_location}/{file_template}_{bead_ndx:>0{order_of_beads}d}.xyz"
                )
                # read corresponding file with ase and
                frames = read(full_path, frame_selection)
                N_frames = len(frames)
                n_frame_list.append(N_frames)
                # Check that all datasets contain the same number of frames
        n_frame_list = np.array(n_frame_list)
        if not np.all(n_frame_list == n_frame_list[0]):
            logger.warning("Different number of frames in datasets ", n_frame_list)
            logger.warning("Setting the number of frames to the minimum value")
        N_frames = np.min(n_frame_list)

    else:
        N_frames = len(list(range(first_frame, last_frame, stride)))
    return N_frames


def get_atomic_numbers(init_file: Path_fr, atomic_numbers: Optional[List[int]]):
    """
    Extract atomic numbers from the init file or use the provided list
    """
    if init_file is None:
        if atomic_numbers is None:
            raise ValueError("Either init_file or atomic_numbers must be specified")
        else:
            atomic_numbers = np.array(atomic_numbers)
    else:
        if atomic_numbers is not None:
            logger.warn(
                "atomic_numbers are disregarded. They will be deduced from the init file instead"
            )
        else:
            init_frame = read(init_file)
            atomic_numbers = np.array(init_frame.numbers)
    return atomic_numbers


def parse_ipi_out_file(file_path: str, fields="all"):
    """
    Parse the name ipi (*.out) file and get data as a dictionary
    """
    # Initialize an empty list for column names
    data_dict = {}
    column_descriptions = []

    # Open the file again to extract column names
    with open(file_path, "r") as file:
        for line in file:
            if line.startswith("#"):
                column_descriptions.append(line.lstrip("#").strip())
            else:
                break  # Stop as soon as you find all column names

    # process column descriptions: need to get index, name and description
    for description in column_descriptions:
        d_split = description.split(":")
        assert len(d_split) == 2
        field = d_split[0].split(">")[1].lstrip().rstrip()
        print(field)
        # ipi has 1-based columns, would like to have 0-based columns
        column = int(d_split[0].split()[1]) - 1
        print(column)

        if (fields == "all") or (field in fields):
            data_dict[field] = np.loadtxt(
                file_path, comments="#", delimiter=None, usecols=column
            )
    return data_dict


def ipi_to_hdf5_single_run(
    dataset_fn: Path_fc,
    P: PositiveInt,
    run_location: Path_drw,
    system: str,
    group_name: str,
    data_dict: Dict,
    init_file: Optional[Path_fr] = None,
    atomic_numbers: Optional[List[int]] = None,
    first_frame: int = 0,
    last_frame: int = -1,
    stride: int = 1,
    temperature: Union[float, int] = None,
    add_spring_energy: bool = False,
    add_potential_energy: bool = False,
    spring_energy_file: Optional[str] = None,
    recompute_spring_energy: bool = False,
    potential_energy_file: Optional[str] = None,
    potential_energy_file_is_replay: bool = False,
):
    """
    Convert a single ipi run into a group in hdf5 dataset

    Parameters:
    -----------
    dataset_fn : str
                Path to a dataset do be created

    P          : int
                Number of replicas in ipi path integral simulation
    run_location : str
                    Location of the run folder

    """
    logger.info(f"Processing data for {P} replicas in PIMD simulation")
    logger.info("Folders used " + str(run_location))

    frame_selection = f"{first_frame}:{last_frame}:{stride}"
    logger.info("The following frames will be used: " + frame_selection)

    atomic_numbers = get_atomic_numbers(
        init_file=init_file, atomic_numbers=atomic_numbers
    )
    message = f"Atomic numbers used: {str(atomic_numbers)}"
    logger.info(message)

    order_of_beads = len(str(P))  # Number of symbols used to encode bead index
    # create bead keys

    bead_keys = [f"bead_{i}" for i in range(P)]

    # First, need to figure number of frames in each dataset.
    # This is only required if number of frames is unknown a priori (i.e., `last_frame` < 0)
    # The reason we have to do it is that sometimes different beads have slightly different
    # number of frames, so we need to take the minimum number of frames over all beads.

    N_frames = get_number_of_frames(
        P=P,
        order_of_beads=order_of_beads,
        frame_selection=frame_selection,
        run_location=run_location,
        data_dict=data_dict,
        first_frame=first_frame,
        last_frame=last_frame,
        stride=stride,
    )
    logger.info("N_frames used: " + str(N_frames))
    # create the file if doesn't exist, open othewise
    f = h5py.File(dataset_fn, "a")
    # create the group
    temp_group = f.create_group(group_name)

    # Get positions and other data into the dataset
    for bead_ndx, bead in enumerate(bead_keys):
        logger.info(f"Process data for  {bead}")
        bead_group = temp_group.create_group(bead)
        for key, properties in data_dict.items():
            file_template = properties[0]
            data_type = properties[1]
            if len(properties) == 2:
                # not replay
                file_is_replay = False
            else:
                file_is_replay = properties[2]
            if file_is_replay:
                start = first_frame + 1
                end = last_frame + 1
                frame_selection = f"{start}:{end}:{stride}"
            else:
                frame_selection = f"{first_frame}:{last_frame}:{stride}"

            full_path = (
                f"{run_location}/{file_template}_{bead_ndx:>0{order_of_beads}d}.xyz"
            )
            # read corresponding file with ase and
            frames = None
            frames = read(full_path, frame_selection)
            all_pos = []
            for frame in frames:
                all_pos.append(frame.get_positions())
            all_pos = np.array(all_pos)[
                :N_frames
            ]  # N_frames is needed for cases when different frames have different number of records
            if data_type == "coordinates":
                all_pos *= Bohr
            elif data_type == "forces":
                all_pos *= Hartree
                all_pos /= Bohr
            else:
                message = f"The data type {data_type} is not supported. Supported values are 'coordinates' or 'forces'"
                logger.error(message)
                raise ValueError(message)
            all_pos = all_pos.astype(np.float32)
            logger.info(f"Created dataset {key} with shape {all_pos.shape}")
            bead_group.create_dataset(key, data=all_pos)
            f.flush()

        bead_group.attrs.create("N_frames", N_frames)
        bead_group.attrs.create("Z", atomic_numbers)

    # Get additional dataset and group-level attributes
    temp_group.attrs.create("temperature", temperature)
    temp_group.attrs.create("P", P)
    temp_group.attrs.create("system", system)

    # Get potential energy from the file

    if add_potential_energy:
        if potential_energy_file is None:
            potential_energy_file_path = f"{run_location}/simulation.out"
        else:
            potential_energy_file_path = f"{run_location}/{potential_energy_file}"
        potential_energy = parse_ipi_out_file(
            potential_energy_file_path, fields=["potential"]
        )["potential"]
        if potential_energy_file_is_replay:
            pot_start = first_frame + 1
            pot_end = last_frame + 1
        else:
            pot_start = first_name
            pot_end = last_frame
        potential_energy = potential_energy[pot_start:pot_end:stride]
        temp_group.create_dataset("potential_energy", data=potential_energy)
        f.flush()

    if add_spring_energy:
        if recompute_spring_energy:
            h_bar = hbar_CI * kJ / 1000 * second
            beta = 1.0 / (temperature * kB)
            logger.info(f"Recomputing spring energy with beta= {beta} a.u.)")
            full_spring_energy = np.zeros(N_frames)
            masses = atomic_masses_legacy[atomic_numbers]
            masses = masses.reshape((1, -1))
            for bead_ndx in range(len(bead_keys)):
                # Get positions for current bead
                # Get positions for the previous bead
                current_bead = bead_keys[bead_ndx]
                previous_bead = bead_keys[bead_ndx - 1]
                current_pos = temp_group[current_bead]["pos"]
                previous_pos = temp_group[previous_bead]["pos"]
                # Compute spring_energy
                spring_energy = compute_spring_energies(
                    current_pos[:],
                    previous_pos[:],
                    masses,
                    P,
                    h_bar,
                    beta,
                    scaling=1.0e0,
                )
                full_spring_energy += np.sum(spring_energy, axis=1)
            dataset_name = "spring"
            full_spring_energy *= Hartree
            temp_group.create_dataset(dataset_name, data=full_spring_energy)
            logger.info(
                f"Created dataset {dataset_name} with shape {full_spring_energy.shape}"
            )

            f.flush()
        else:
            if spring_energy_file is None:
                spring_energy_file_path = f"{run_location}/simulation.out"
            else:
                spring_energy_file_path = f"{run_location}/{spring_energy_file}"
            logger.info(f"Extracting spring energies from: {spring_energy_file_path}")
            spring_energy = parse_ipi_out_file(
                spring_energy_file_path, fields=["spring"]
            )["spring"]
            spring_energy = spring_energy[first_frame:last_frame:stride] * Hartree
            dataset_name = "spring"
            temp_group.create_dataset(dataset_name, data=spring_energy)
            logger.info(
                f"Created dataset {dataset_name} with shape {spring_energy.shape}"
            )

            f.flush()
    f.close()
    message = f"Group {group_name} was processed successfully"
    logger.info(message)
    return


def ipi_to_hdf5_converter(
    dataset_fn: Path_fc,
    P: PositiveInt,
    root_location: Path_drw,
    system: str,
    temp_folders: List,
    data_dict: Dict,
    init_file: Optional[Path_fr] = None,
    atomic_numbers: Optional[List[int]] = None,
    first_frame: int = 0,
    last_frame: int = -1,
    stride: int = 1,
    temperatures: Union[List[float], List[int]] = None,
    add_spring_energy: bool = False,
    add_potential_energy: bool = False,
    spring_energy_file: Optional[str] = None,
    recompute_spring_energy: bool = False,
    potential_energy_file: Optional[str] = None,
    potential_energy_file_is_replay: bool = False,
):
    """
        Convert output of ipi PIMD simulation into an hdf5 dataset.
        It is assumed, that ipi results are stored in folders with the
        following structure:
        root_location/temp_folders[@]/output_file

        Parameters:
        ----------
        dataset_fn : str
                 Path to a dataset do be created
    P          : int
                 Number of replicas in ipi path integral simulation
    root_location : str
                   Location of the root folder
    temp_folder: List[str]
                 List of strings, representing subfolders in the root_location.
                 Each subfolder should contain the results of the simulation for a single system
    atomic_numbers: List[int]
                   List of atomic numbers in order as they are in the frame
    init_file: Optional[Path_fr]
               A single frame file corresponding to the data. If provided, is used to deduce atomic numbers and will
               take priority over atomic_numbers
    data_dict: Dict
               Dictionary of pairs key:value, where each key represents a key for the dataset,
               and value is  a two- or three- element list. In this list, the first
               element represents a file name template for the ipi output
               that contains corresponding data, and the second represents type of
               the data (coordinates or force). The third element is optional, and is  a boolean value. If set to true, the corresponding
               files are replays, so the first frame of them will be ignored.
               If  Data about type are used to convert
               from ipi units (Bohr/Hartree) to ASE units (Angstrom/eV)
               file. It is assumed that the ipi output files are named as
               <template>_<bead index}>.xyz

    first_frame : int
               Index of the first frame to use, default 0

    last_frame : int
               Index of the last frame to use, default -1

    add_spring_energy: bool, default False
                          If True, add spring energy to the dataset.
                          Requires that the spring energy is present
                          in the ipi output
    add_potential_energy: bool, default False
                            If True, add potential energy to the dataset.
                            Requires that the potential energy is present
                            in the ipi output.

    Example configuration
    ----------------------
    P = 64  # Number of beads
    dataset_fn = 'dataset_d5o2+.h5'
    system = ['d5o2+']
    atomic_numbers = np.array([8, 8, 1, 1, 1, 1 ,1])
    root_location = '/scratch/iryna/PI_CG/zundel_cation/ZundelD_sampling' # Location of the file
    temp_folders = ['0100', '0200', '0300', '0400', '0500', '0600']
    # data_dict contains pairs  template :  data_key, where
    # template represents file name template  <template>_XX.xyz, data_key represents key that should be used
    # for these data in the dataset

    data_dict = {'pos': ['simulation.pos', 'coordinate'],
                 'for': ['simulation.for', 'force', true],
                 'fsp': ['simulation.fsp', 'force']}

    first_frame  =  0
    last_frame = -1
    """
    # ###################################################
    groups = [f"{system}_T={i}K" for i in temp_folders]

    for temp_folder, group, temperature in zip(temp_folders, groups, temperatures):
        logger.info(f"Processing data for {group}")
        run_location = f"{root_location}/{temp_folder}"
        ipi_to_hdf5_single_run(
            dataset_fn=dataset_fn,
            P=P,
            run_location=run_location,
            system=system,
            group_name=group,
            data_dict=data_dict,
            init_file=init_file,
            atomic_numbers=atomic_numbers,
            first_frame=first_frame,
            last_frame=last_frame,
            stride=stride,
            temperature=temperature,
            add_spring_energy=add_spring_energy,
            add_potential_energy=add_potential_energy,
            spring_energy_file=spring_energy_file,
            recompute_spring_energy=recompute_spring_energy,
            potential_energy_file=potential_energy_file,
            potential_energy_file_is_replay=potential_energy_file_is_replay,
        )

    logger.info("All groups were processed successfully")


if __name__ == "__main__":

    # Set up logging configuration
    logging.basicConfig(
        level=logging.DEBUG,
        format="%(asctime)s %(levelname)s %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    logger = logging.getLogger()
    log_execution_info(logger)

    CLI(ipi_to_hdf5_converter)
    logger.info("Execution has ended successfully")

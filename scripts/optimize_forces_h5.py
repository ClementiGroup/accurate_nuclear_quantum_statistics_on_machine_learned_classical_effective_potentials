"""
Script takes h5 dataset with spring and classical forces already in place,
and performs force optimization. Resulting datasets are writen back to the 
original dataset
"""

import h5py
import numpy as np


from aggforce import linearmap as lm
from aggforce import agg as ag
from aggforce import constfinder as cf
from aggforce import featlinearmap as p
from jsonargparse import CLI
from typing import Union, List, Optional
import logging

from log_utils import log_execution_info

logger = logging.getLogger(__name__)


def index_replica(N, N_atoms):
    """
    Get indexes that correspond to replica N with N atoms
    """
    inds = [[i] for i in range(N * N_atoms, (N + 1) * N_atoms)]
    return inds


def optimize_forces(
    P: int,
    filename: str,
    t0: Union[float, int],
    temperatures: List[Union[float, int]],
    groups: List[str],
    split_by_molecule: bool = False,
    split: Optional[int] = None,
    dry_run: Optional[bool] = False,
    optimized_group: Optional[str] = "total_forces",
    eps_abs: float = 1e-8,
    l2_regularizations: List[float] = [5e-2],
    optimization_force_stride: int = 1,
    transform_in_chunks: bool = False,
    chunk_size: Optional[int] = None,
):
    """
    Optimize forces from the hdf5 dataset
    P : int
      Number of beads in the path integral
    filename : str
       Path to the dataset
    t0 : float
        Base temperature used for prior rescaling
    groups: list
         List of groups from the dataset to optimize.
    temperatures: list
        List of temperatures for groups listed in `groups'
    split_by_molecule: bool, default False
        If True, identical molecules will share the same map.
        How exactly the molecules are defined is given by split
    split: int, optional
        Is used when split_by_molecule = True and defines number of
        atoms in a molecule. It is assumed that the system consists of identical
        molecules, listed sequentially and with atoms within each molecule
        listed in the same order.
    dry_run: bool, optional, default=False
        if True, optimized forces are printed out, but no datasets are created
    optimized_group: str, default='total_forces'
        Group used to derive optimum force mapping. In most cases, should be  total_forces.
        A different group can be used mostly for debugging
    eps_abs: float, 1e-8
        Convergence paramter of the optimizer
    l2_regularizations: list of float, default [1e-1, 5e-2]
        Regularization values to use during optimization
    optimization_force_stride: 1
        Striding forces before the optimization
    transform_in_chunks: bool, default False
        If true, the transformation matrix will be applied in chunks. Usefull when calculations
        are expensive
    chunk_size: Optional[int]
        Number of frames in a chunk

    """
    f = h5py.File(filename, "a")
    # Assume, that groups of interest are
    bead_keys = [f"bead_{i}" for i in range(P)]

    if optimized_group != "total_forces":
        logger.warning(
            "You optimized forces based on {}, and not total_forces !. Make sure you know what you are doing"
        )
    if split_by_molecule:
        assert split is not None

    for l2_regularization in l2_regularizations:
        for t, group in zip(temperatures, groups):
            logger.info(f"Group: {group}")
            scaling = 1.0 * t / t0
            logger.info(f"Prior scaling: {scaling}")
            total_force_list = []
            for bead in bead_keys:
                total_forces = f[group][bead].get(optimized_group)
                dataset_shape = total_forces.shape
                N_atoms = total_forces.shape[1]
                if split_by_molecule:
                    assert (
                        N_atoms % split == 0
                    ), "Total number of atoms is not devided by number of atoms per molecule"
                    n_splits = N_atoms // split
                    logger.info(f"Splitting is used. Number of splits: {n_splits}")
                    reshaped_forces = np.concatenate(
                        np.split(total_forces, n_splits, axis=1), axis=0
                    )
                    total_force_list.append(reshaped_forces)
                else:
                    total_force_list.append(total_forces)
            full_traj_total_forces = np.concatenate(total_force_list, axis=1)
            logger.info(
                f"The total force trajectory has shape {full_traj_total_forces.shape}"
            )
            del total_force_list
            # Do force optimization
            for resndx, bead in enumerate(bead_keys):
                logger.info(f"Bead: {bead}")
                logger.info(f"split: {split}")
                if split_by_molecule:
                    inds = index_replica(resndx, split)
                    assert len(inds) == split
                else:
                    inds = index_replica(resndx, N_atoms)
                    assert len(inds) == N_atoms
                logger.debug(f"Indexing of the replica: {inds}")
                cmap = lm.LinearMap(inds, n_fg_sites=full_traj_total_forces.shape[1])
                try:
                    optim_results = ag.project_forces(
                        xyz=None,
                        forces=full_traj_total_forces[
                            ::optimization_force_stride, :, :
                        ],
                        config_mapping=cmap,
                        constrained_inds={},
                        l2_regularization=l2_regularization,
                        solver_args=dict(solver="osqp", eps_abs=eps_abs),
                    )
                    new_map = optim_results["map"]
                    matrix = new_map.standard_matrix
                except Exception as error:
                    logger.warning(f"Optimization was not successfull: {error}")
                    continue
                logger.info(f"Transform map has shape: {matrix.shape}")
                logger.debug(f"Transform matrix: {matrix}")

                if not transform_in_chunks:
                    total_forces_projected = np.matmul(matrix, full_traj_total_forces)
                    logger.info(
                        f"Projected total forces have shape: {total_forces_projected.shape}"
                    )
                    # Need to reshape total_forces_projected
                    if split_by_molecule:
                        total_forces_projected = np.concatenate(
                            np.split(total_forces_projected, n_splits, axis=0), axis=1
                        )
                    if dry_run:
                        logger.info("Running with dry run")
                        logger.info(str(total_forces_projected))
                        continue

                    total_forces_projected_scaled_classical_prior = (
                        total_forces_projected - f[group][bead].get("for")[:] * scaling
                    )
                    try:
                        f[group][bead].create_dataset(
                            f"total_forces_projected_scaled_classical_prior_scaling_{scaling:.3e}_l2reg_{l2_regularization:.3e}",
                            data=total_forces_projected_scaled_classical_prior,
                        )
                    except ValueError:
                        print("Dataset already exists")

                    f.flush()
                else:
                    N_frames = full_traj_total_forces.shape[0]
                    # Each chunk should contain an integer number of splits
                    # Create all the datasets in advance

                    assert chunk_size % n_splits == 0
                    if not dry_run:
                        try:
                            total_forces_projected_scaled_classical_prior_ds = f[group][
                                bead
                            ].create_dataset(
                                f"total_forces_projected_scaled_classical_prior_scaling_{scaling:.3e}_l2reg_{l2_regularization:.3e}",
                                shape=dataset_shape,
                                dtype=np.float32,
                            )
                            total_forces_projected_scaled_classical_prior_filled = False
                        except ValueError:
                            print("Dataset already exists")
                            total_forces_projected_scaled_classical_prior_filled = False
                            total_forces_projected_scaled_classical_prior_ds = f[group][
                                bead
                            ][
                                f"total_forces_projected_scaled_classical_prior_scaling_{scaling:.3e}_l2reg_{l2_regularization:.3e}"
                            ]

                    for i in range(0, N_frames, chunk_size):
                        logger.info(f"Processing frames {i}:{i+chunk_size}")
                        current_chunk_size = min(chunk_size, N_frames - i)
                        total_forces_projected_chunk = np.matmul(
                            matrix,
                            full_traj_total_forces[i : i + current_chunk_size, :, :],
                        )
                        # Need to reshape total_forces_projected
                        if split_by_molecule:
                            total_forces_projected_chunk = np.concatenate(
                                np.split(
                                    total_forces_projected_chunk, n_splits, axis=0
                                ),
                                axis=1,
                            )

                            logger.info(
                                f"Processed Projected total forces have shape: {total_forces_projected_chunk.shape}"
                            )
                            reshaped_chunk_size = int(current_chunk_size / n_splits)
                            logger.info(
                                f"The reshaped chunk_size: {chunk_size/n_splits}"
                            )
                            loc = int(i / n_splits)
                            logger.info(f"locator position: {i/n_splits}")
                        else:
                            loc = i
                            reshaped_chunk_size = current_chunk_size
                        if dry_run:
                            logger.info("Running with dry run")
                            logger.info(
                                f"Shape of projected total forces: {total_forces_projected_chunk.shape}"
                            )
                            continue

                        if not total_forces_projected_scaled_classical_prior_filled:
                            total_forces_projected_scaled_classical_prior_chunk = (
                                total_forces_projected_chunk
                                - f[group][bead].get("for")[
                                    loc : loc + reshaped_chunk_size, :, :
                                ]
                                * scaling
                            )
                            total_forces_projected_scaled_classical_prior_ds[
                                loc : loc + reshaped_chunk_size, :, :
                            ] = total_forces_projected_scaled_classical_prior_chunk
                            del total_forces_projected_scaled_classical_prior_chunk

                    f.flush()
    f.close()
    logger.info("The execution has finished successfully")
    return ()


if __name__ == "__main__":
    # Set up logging configuration
    logging.basicConfig(
        level=logging.DEBUG,
        format="%(asctime)s %(levelname)s %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    logger = logging.getLogger(__name__)
    log_execution_info(logger)
    CLI(optimize_forces)
    logger.info("Execution has ended successfully")

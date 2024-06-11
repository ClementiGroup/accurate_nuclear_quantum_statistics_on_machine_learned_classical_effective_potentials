"""
Scripts takes a  a list of model files and extracts internal coordinates for them
"""

from ase.io import read
from ase.geometry import get_distances, get_dihedrals

import numpy as np
from jsonargparse import CLI


def extract_cv_h2o(position_file):
    """
    Extracts data from a position file
    """
    atomic_positions = read(position_file, ":")
    reference_frame = atomic_positions[0]  #
    symbols = reference_frame.get_chemical_symbols
    data = {}
    data["doh1"] = np.concatenate([i.get_distances(0, 1) for i in atomic_positions])
    data["doh2"] = np.concatenate([i.get_distances(0, 2) for i in atomic_positions])
    data["dh1h2"] = np.concatenate([i.get_distances(1, 2) for i in atomic_positions])
    data["angles"] = np.array([i.get_angle(1, 0, 2) for i in atomic_positions])
    return data


def analyze_zundel(atomic_object):
    """
    Function takes a single ase atomic object and

    These quantities are:
     - all the OH distances
     - O-O distance
     - Proton transfere coordinate
     - the torsion angle φ with respect to XA-OA-OB-XB
       where XA and XB are the midpoints of two non-hydrogen
       bonded atoms in the same H3O+fragment
     - out-of-plane distance d(Oa) and d(Ob) measured between the
       corresponding O atom and the plane of three hydrogen atoms
       comprising the H3O+ fragment

    The coordinates of interest are defined as

    Marx, D., Tuckerman, M., Hutter, J. et al.
    The nature of the hydrated excess proton in water.
    Nature 397, 601–604 (1999). https://doi.org/10.1038/17579
    (see figure 2)

    These quantities are:
     - all the OH and O-O distances
     - Proton transfere coordinates

    The function:

     - determines indexes of O and H (assume the order is the same in all the frames)
     - calculate O-O distances and all the O-H distances
     - extract proton transfer coordinates
     - calculates dihedral angle, tat defines rotation
    """
    atomic_positions = atomic_object.positions
    symbols = atomic_object.get_chemical_symbols()

    # find indices of O and H
    o_ndxs = [i for i, s in enumerate(symbols) if s == "O"]
    h_ndxs = [i for i, s in enumerate(symbols) if s == "H"]
    assert len(o_ndxs) == 2
    assert len(h_ndxs) == 5

    #  calculate O-O distance
    d_oo = atomic_object.get_distances(*o_ndxs)

    #  calculate O-H distances
    #  First, will compute distances for all O-H pairs

    # calculate O-H distances
    # distances o_h: List of two lists.
    # The first list corresponds to O1-H distances( to the 0th O atom)
    # The second list  corrsspons to O2-H distances (to the 1st atom)
    # H atoms should be in the same order for both inner lists
    positions_o = atomic_positions[o_ndxs]
    positions_h = atomic_positions[h_ndxs]

    # ase.geometry.get_distances returns distance matrix (element 0)
    # and distances themselves (element 1)
    distances_oh = get_distances(positions_o, positions_h)[1]

    # calculate H-H distances
    distances_hh = get_distances(positions_h)[1]
    # bounded_distances_oh contains distance to the closest O for each H
    bounded_distances_oh = distances_oh.min(axis=0)
    print(bounded_distances_oh)
    # Index of the oxygen bound to a particular hydrogen
    bound2oxygen_ndx = distances_oh.argmin(axis=0)

    # PROTON TRANSFER COORDINATE
    # Distance difference for all the
    delta_all = distances_oh[1] - distances_oh[0]

    # index of the transfered H
    delta_ndx = np.argmin(np.abs(delta_all))

    # Protein transfewr coordinate
    delta = delta_all[delta_ndx]

    # DIHEDRAL ANGLE
    # First, need to determine pairs of H, that are not proton-transfer
    # List of indexes of H atoms, bound to the 0th and 1st atoms
    # (excluding the transfer H atom)
    h_o0_ndx = []
    h_o1_ndx = []
    for ndx, val in enumerate(bound2oxygen_ndx):
        if ndx == delta_ndx:
            continue
        elif val == 0:
            h_o0_ndx.append(ndx)
        elif val == 1:
            h_o1_ndx.append(ndx)
        else:
            raise (ValueError)

    try:
        assert (
            len(h_o0_ndx) == 2
        ), f"Wrong number of H bound to the Oxygen 0 : {len(h_o0_ndx)} instead of 2 {delta_ndx}"
        assert len(h_o1_ndx) == 2, "Wrong number of H bound to the Oxygen 1"
        # Calculate middle points
        XA = np.mean(positions_h[np.array(h_o0_ndx)], axis=0)
        XB = np.mean(positions_h[np.array(h_o1_ndx)], axis=0)

        # Get vectors for dihedral calculation XA - O0 - O1 - XB
        v1 = positions_o[0] - XA
        v2 = positions_o[1] - positions_o[0]
        v3 = XB - positions_o[1]

        dihedral = get_dihedrals([v1], [v2], [v3])
        phi = min(dihedral, 360 - dihedral)

        # OUT-OF-PLAIN DISTANCES
        # Will consider deltaH as the center of coordinates, so
        # distance between plane spanned by 3 H and O can be calculated as
        # abs((v1 x v2) dot v3/ || v1 x v2||, where  v1, v2 are vectors from transfer H to bound H,
        # and v3 is a vector from transfer H to the corresponding O atom
        # Atom Oo
        v1 = positions_h[h_o0_ndx[0]] - positions_h[delta_ndx]
        v2 = positions_h[h_o0_ndx[1]] - positions_h[delta_ndx]
        v3 = positions_o[0] - positions_h[delta_ndx]
        v1crossv2 = np.cross(v1, v2)
        plane_norm = v1crossv2 / np.linalg.norm(v1crossv2)
        dO0 = np.abs(np.dot(plane_norm, v3))

        # Atom O1
        v1 = positions_h[h_o1_ndx[0]] - positions_h[delta_ndx]
        v2 = positions_h[h_o1_ndx[1]] - positions_h[delta_ndx]
        v3 = positions_o[1] - positions_h[delta_ndx]
        v1crossv2 = np.cross(v1, v2)
        plane_norm = v1crossv2 / np.linalg.norm(v1crossv2)
        dO1 = np.abs(np.dot(plane_norm, v3))

    except:
        print("UNRELIABLE ANALYSIS: wrong OH connectivity")
        delta = np.nan
        phi = np.nan
        dO0 = np.nan
        dO1 = np.nan

    return d_oo, distances_oh, distances_hh, delta, phi, dO0, dO1


def extract_cv_zundel_cat(position_file):
    """
    Extracts data from a position file for zundel cation
    """
    atomic_positions = read(position_file, ":")
    d_oo_total = []
    distances_oh_total = []
    distances_hh_total = []
    delta_total = []
    phi_total = []
    dO0_total = []
    dO1_total = []

    for structure in atomic_positions:
        d_oo, distances_oh, distances_hh, delta, phi, dO0, dO1 = analyze_zundel(
            structure
        )

        d_oo_total.append(d_oo)
        distances_oh_total.append(distances_oh)
        distances_hh_total.append(distances_hh)
        delta_total.append(delta)
        phi_total.append(phi)
        dO0_total.append(dO0)
        dO1_total.append(dO1)

    data = {}
    data["doo"] = np.array(d_oo_total).flatten()
    data["distances_oh"] = np.array(distances_oh_total)
    data["distances_hh"] = np.array(distances_hh_total)
    data["delta"] = np.array(delta_total).flatten()
    data["phi"] = np.array(phi_total).flatten()
    data["dO0"] = np.array(dO0_total).flatten()
    data["dO1"] = np.array(dO1_total).flatten()
    return data


systems = {"h2o": extract_cv_h2o, "zundel_cat": extract_cv_zundel_cat}


def run_analysis(
    position_file: str, output_file: str = "data.npz", system: str = "h2o"
):

    cv_function = systems[system]
    data = cv_function(position_file)
    np.savez(output_file, **data)
    return


if __name__ == "__main__":
    CLI(run_analysis)
import argparse
import numpy as np
from HeuslerPlot.search import FindEs, FindBands
from HeuslerPlot.parseVASP import ParseOutcar, ParseOszicar, ParseEigenval

def PlotBands(ks, eigenvals, magmom, E_Fermi, k_labels, R, out_path):
    # Remove duplicate k-point pairs; they correspond to symmetry points
    # (where VASP moves from one k1->k2 path to another).
    # Note the indices of the symmetry points in the new k and eigenval lists.
    sym_indices, ks_cut, eigenvals_cut = _cut_duplicates(ks, eigenvals)

    recip_dists = _recip_dist(sym_indices, ks_cut, R)

    xs = _scaled_k_xs(sym_indices, ks_cut, recip_dists)

def _cut_duplicates(ks, eigenvals):
    '''Return lists sym_indices, ks_cut, and eigenvals_cut, where ks_cut and
    eigenvals_cut have had one member of duplicate pair entries removed and
    sym_indices gives indices in the cut lists which correspond to symmetry
    points (at the ends and where the duplicate pairs were).
    '''
    # First symmetry point is at the initial k-point.
    sym_indices = [0]
    ks_cut, eigenvals_cut = [], []

    for ks_index, k in enumerate(ks):
        # Always keep the first point.
        if ks_index == 0:
            ks_cut.append(k)
            eigenvals_cut.append(eigenvals[0])
            continue
        # If we get here, ks_index > 0.
        # Check if this point is the same as the last one.
        prev_k = ks[ks_index - 1]
        if _vec_equal(k, prev_k):
            # We are at the second point in a duplicate pair.
            # Cut out this point (by not adding it to ks_cut, eigenvals_cut)
            # and note the location of the first point of the pair.
            sym_indices.append(len(ks_cut) - 1)
        else:
            # Not at the second point in a duplicate pair.
            ks_cut.append(k)
            eigenvals_cut.append(eigenvals[ks_index])

    # Last symmetry point is at the last k-point.
    sym_indices.append(len(ks_cut) - 1)

    return sym_indices, ks_cut, eigenvals_cut

def _recip_dist(sym_indices, ks_cut, R):
    '''Return a list of values which give the Cartesian distance between each
    pair of symmetry point.
    '''
    dists = []
    for point_index, k_index in enumerate(sym_indices):
        # Skip first point (making [k, previous k] pairs).
        if point_index == 0:
            continue

        k_Cart = np.dot(ks_cut[k_index], R)
        prev_k_index = sym_indices[point_index-1]
        prev_k_Cart = np.dot(ks_cut[prev_k_index], R)

        k_to_prev_k = np.subtract(prev_k_Cart, k_Cart)
        this_dist = np.linalg.norm(k_to_prev_k)
        dists.append(this_dist)

    return dists

def _scaled_k_xs(sym_indices, ks_cut, recip_dists):
    total_dist = sum(recip_dists)
    xs = []

    current_x = 0.0
    base_x = 0.0
    panel_start_index = 0
    panel_number = 0
    for x_index in range(len(ks_cut)):
        points_in_panel = sym_indices[panel_number+1] - sym_indices[panel_number] + 1
        step = (recip_dists[panel_number] / total_dist) / (points_in_panel - 1)
        current_x = base_x + (x_index - panel_start_index)*step
        xs.append(current_x)

        if x_index in sym_indices and x_index != 0:
            base_x = current_x
            panel_start_index = x_index
            panel_number += 1

    return xs

def swap_channels_if_mag_neg(magmom, eigenvals):
    '''If magmom < 0, eigenvals[k_index][0] gives "down" eigenvalues and
    eigenvals[k_index][1] gives "up" eigenvalues. To have consistant up/down
    identification, need to swap these (since the sign of magmom is arbitrary
    in the absence of an applied field).

    Returns a new eigenval list with this swap made if necessary.
    '''
    if magmom < 0.0:
        corrected = []
        for k_index in range(len(eigenvals)):
            corrected.append((eigenvals[k_index][1], eigenvals[k_index][0]))
        return corrected
    else:
        return eigenvals

# Perform exact equality comparison between vectors.
# --> Assumes no floating-point error is present; this may be the case if u
#     and v are obtained from reading a file (where u and v are reported with
#     some not-too-high precision, less than machine epsilon).
def _vec_equal(u, v):
    if len(u) != len(v):
        raise ValueError("Comparing vectors of unequal length")
    for i in range(len(u)):
        if u[i] != v[i]:
            return False
    return True

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot VASP bands.")
    parser.add_argument('dir_path', type=str, help="Path to directory with data to plot.")
    parser.add_argument('--searchEs', action='store_true',
            help="Search subdirectories of dir_path of the form e##.")
    args = parser.parse_args()

    all_data_paths = None
    if args.searchEs:
        all_data_paths = FindEs(args.dir_path)
    else:
        all_data_paths = FindBands(args.dir_path)

    for system_name, system_data_paths in all_data_paths.items():
        E_Fermi, D = ParseOutcar(system_data_paths['outcar_path'])
        magmom = ParseOszicar(system_data_paths['oszicar_path'])
        ks, eigenvals = ParseEigenval(system_data_paths['eigenval_path'])
        
        # The rows of R (the R[i, :]) are the reciprocal lattice vectors.
        R = 2.0 * np.pi * np.linalg.inv(D)

        # TODO - get this from KPOINTS.
        k_labels = ["$\Gamma$", "$X$", "$W$", "$K$", "$\Gamma$", "$L$", "$W$", "$U$", "$X$"]

        # TODO - put somewhere else? Other structure to name?
        out_path = system_name

        nspin = len(eigenvals[0])
        if nspin == 2:
            eigenvals = swap_channels_if_mag_neg(magmom, eigenvals)

        PlotBands(ks, eigenvals, magmom, E_Fermi, k_labels, R, out_path)

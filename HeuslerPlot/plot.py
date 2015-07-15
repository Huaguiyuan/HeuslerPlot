import argparse
import numpy as np
from HeuslerPlot.search import FindEs, FindBands

def PlotBands(ks, eigenvals, E_Fermi, k_labels, R, out_path):
    # Remove duplicate k-point pairs; they correspond to symmetry points
    # (where VASP moves from one k1->k2 path to another).
    # Note the indices of the symmetry points in the new k and eigenval lists.
    sym_indices, ks_cut, eigenvals_cut = _cut_duplicates(ks, eigenvals)

def _cut_duplicates(ks, eigenvals):
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

        PlotBands(ks, eigenvals, E_Fermi, k_labels, R, out_path)

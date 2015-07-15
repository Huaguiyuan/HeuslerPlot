import argparse
import numpy as np
import matplotlib.pyplot as plt
from HeuslerPlot.search import FindEs, FindBands
from HeuslerPlot.parseVASP import ParseOutcar, ParseOszicar, ParseEigenval

def PlotBands(ks, eigenvals, magmom, E_Fermi, k_labels, R, out_path, logo_text=None):
    # Remove duplicate k-point pairs; they correspond to symmetry points
    # (where VASP moves from one k1->k2 path to another).
    # Note the indices of the symmetry points in the new k and eigenval lists.
    sym_indices, ks_cut, eigenvals_cut = _cut_duplicates(ks, eigenvals)

    # Shift values so that E_F = 0.0.
    eigenvals_cut = _shift_Fermi(eigenvals_cut, E_Fermi)

    # Get x-positions of k values, scaled appropriately by the distance
    # between the symmetry points at the ends of panels.
    recip_dists = _recip_dist(sym_indices, ks_cut, R)
    xs = _scaled_k_xs(sym_indices, ks_cut, recip_dists)

    sym_xs = []
    for i in range(len(sym_indices)):
        sym_xs.append(xs[sym_indices[i]])

    # Transform eigenvals[k_index][spin][band_index] into
    # eigenval_ys[spin][band_index][k_index] so we can use it as an argument
    # for plt.plot(xs, ys).
    eigenval_ys = _make_eigenval_series(eigenvals_cut)

    nspin = len(eigenvals[0])
    nbands = len(eigenvals[0][0])
    if nspin == 1:
        fig = plt.figure(figsize=(6, 5))

        plt.ylabel("Energy (eV)")
        plt.xlim(0.0, 1.0)
        plt.ylim(-10.0, 10.0)
        plt.xticks(sym_xs, k_labels)
        plt.yticks(np.arange(-10.0, 10.0, 2.0))
        plt.grid(b=True)

        for b_i in range(nbands):
            plt.plot(xs, eigenval_ys[0][b_i], 'k')

        if logo_text != None:
            pass

        plt.savefig(out_path + '.png', bbox_inches='tight', dpi=500)
    else:
        fig = plt.figure(figsize=(12, 5))

        up_plot = plt.subplot(121)
        plt.title("Up Spin")
        plt.ylabel("Energy (eV)")
        plt.xlim(0.0, 1.0)
        plt.ylim(-10.0, 10.0)
        plt.xticks(sym_xs, k_labels)
        plt.yticks(np.arange(-10.0, 10.0+1e-6, 2.0))
        plt.grid(b=True, linestyle='-', linewidth=0.5, axis='x')
        plt.grid(b=True, axis='y', linewidth=0.5)

        for b_i in range(nbands):
            plt.plot(xs, eigenval_ys[0][b_i], 'k')

        if logo_text != None:
            plt.annotate(logo_text, (0.155, 0.125), xycoords='figure fraction', size=12)

        down_plot = plt.subplot(122)
        plt.title("Down Spin")
        plt.xlim(0.0, 1.0)
        plt.ylim(-10.0, 10.0)
        plt.xticks(sym_xs, k_labels)
        plt.yticks(np.arange(-10.0, 10.0, 2.0))
        plt.grid(b=True, linestyle='-', linewidth=0.5, axis='x')
        plt.grid(b=True, axis='y', linewidth=0.5)

        for b_i in range(nbands):
            plt.plot(xs, eigenval_ys[1][b_i], 'k')

        if logo_text != None:
            plt.annotate(logo_text, (0.6545, 0.125), xycoords='figure fraction', size=12)

        plt.savefig(out_path + '.png', bbox_inches='tight', dpi=500)

def _shift_Fermi(eigenvals, E_Fermi):
    shifted_evs = []
    for k_i in range(len(eigenvals)):
        shifted_evs.append([])
        for s in range(len(eigenvals[0])):
            shifted_evs[k_i].append([])
            for b_i in range(len(eigenvals[0][0])):
                shifted_evs[k_i][s].append(eigenvals[k_i][s][b_i] - E_Fermi)
    return shifted_evs

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

def _make_eigenval_series(eigenvals):
    '''Return a copy of eigenvals rearranged such that it is addressed
    as eigenval_series[spin][band_index][k_index].
    The eigenval_series[spin][band_index] values can then be used as ys in
    plt.plot(xs, ys).
    '''
    nk = len(eigenvals)
    nspin = len(eigenvals[0])
    nbands = len(eigenvals[0][0])

    eigenval_series = []
    for s in range(nspin):
        eigenval_series.append([])
        for b_i in range(nbands):
            eigenval_series[s].append([])

    for k_i in range(nk):
        for s in range(nspin):
            for b_i in range(nbands):
                # Make sure that data is sorted (avoid jumps).
                # Should already be sorted in VASP output, but can
                # double-check here.
                eigenvals_k_s_sorted = sorted(eigenvals[k_i][s])
                # Put in new order.
                eigenval_series[s][b_i].append(eigenvals_k_s_sorted[b_i])

    return eigenval_series

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
    parser.add_argument('--logo', action='store_true',
            help="Add Heusler site URL to plot")
    parser.add_argument('--structure_type', default='L21',
            help="Type of structure contained in subdirectories of dir_path")
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

        # TODO - configurable destination (instead of pwd)?
        out_path = "{}_{}".format(system_name, args.structure_type)

        nspin = len(eigenvals[0])
        if nspin == 2:
            eigenvals = swap_channels_if_mag_neg(magmom, eigenvals)

        logo_text = None
        if args.logo:
            logo_text = "www.heusleralloys.mint.ua.edu"

        PlotBands(ks, eigenvals, magmom, E_Fermi, k_labels, R, out_path, logo_text)

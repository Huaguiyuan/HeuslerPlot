import argparse
import numpy as np
import matplotlib.pyplot as plt
from HeuslerPlot.search import FindEs, FindBands
from HeuslerPlot.parseVASP import ParseOutcar, ParseOszicar, ParseEigenval

def PlotBands(ks, eigenvals, E_Fermi, k_labels, R, out_path, logo_text=None):
    '''Plot the bands given by ks and eigenvals. Plotted eigenvals will
    be shifted so that E = 0 corresponds to E_Fermi. Symmetry points are
    labeled by k_labels, and panels are sized such that their width
    corresponds to the length of the corresponding path in Cartesian
    coordinates. Optionally, a text logo can be displayed over the plot.

    ks = a list of k-points with the form returned by parseVASP.ParseEigenval.
        The elements have the form (ka, kb, kc) giving the value of the
        corresponding k-point in reciprocal lattice coordinates.
    eigenvals = a list of eigenvalues with the form returned by
        parseVASP.ParseEigenval. Eigenvalues are stored with the structure
        eigenvals[k_index][spin_index][band_index]. If the calculation
        is non-spin-polarized or has noncollinear spins, spin_index is
        always 0. If the calculation is spin-polarized with collinear
        spins, spin_index may be 0 or 1 (corresponding to up- and down-states).
    E_Fermi = the Fermi energy, used to shift plotted eigenvalues such that
        E = 0 corresponds to E_Fermi.
    k_labels = a list of strings labelling the symmetry points used in the plot.
        Symmetry points are assumed to lie at the beginning and end of the plot
        and at all points in between at which the same k-point occurs twice
        in a row.
    R = a 3x3 numpy array with rows giving the reciprocal lattice vectors
        (R[i, :] corresponds to the i'th reciprocal lattice vector).
    out_path = path at which the plot file will be written.
    logo_text = if this is not None, its value is written on the plot.
    '''
    # TODO - generalize from 'doubled k-points = symmetry point' approach
    # to allow for discontinuous panels.
    # Add an optional argument giving the number of k-points per panel in
    # eigenvals?

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
        plt.yticks(np.arange(-10.0, 10.0+1e-6, 2.0))
        plt.grid(b=True)

        for b_i in range(nbands):
            plt.plot(xs, eigenval_ys[0][b_i], 'k')

        if logo_text != None:
            # TODO - logo position for nspin = 1.
            pass

        plt.savefig(out_path + '.png', bbox_inches='tight', dpi=500)
        plt.close('all')
    else:
        # TODO - make figure size an argument?
        fig = plt.figure(figsize=(16, 5))

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
        plt.yticks(np.arange(-10.0, 10.0+1e-6, 2.0))
        plt.grid(b=True, linestyle='-', linewidth=0.5, axis='x')
        plt.grid(b=True, axis='y', linewidth=0.5)

        for b_i in range(nbands):
            plt.plot(xs, eigenval_ys[1][b_i], 'k')

        if logo_text != None:
            plt.annotate(logo_text, (0.6545, 0.125), xycoords='figure fraction', size=12)

        plt.savefig(out_path + '.png', bbox_inches='tight', dpi=500)
        plt.close('all')

def _shift_Fermi(eigenvals, E_Fermi):
    '''Return a copy of eigenvals where E_Fermi has been subtracted from
    each value. In the returned set of eigenvalues, E = 0 corresponds to
    E_Fermi.
    '''
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
    eps = 1e-9 # consider k-points equal if components differ by less than eps

    for ks_index, k in enumerate(ks):
        # Always keep the first point.
        if ks_index == 0:
            ks_cut.append(k)
            eigenvals_cut.append(eigenvals[0])
            continue
        # If we get here, ks_index > 0.
        # Check if this point is the same as the last one.
        prev_k = ks[ks_index - 1]
        if _vec_approx_equal(k, prev_k, eps):
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
    pair of consecutive symmetry points.
    '''
    dists = []
    for point_index, k_index in enumerate(sym_indices):
        # Skip first point (need to make [k, previous k] pairs).
        if point_index == 0:
            continue
        # Get k and previous k (prev_k).
        k = ks_cut[k_index]
        prev_k_index = sym_indices[point_index-1]
        prev_k = ks_cut[prev_k_index]
        # Convert k and prev_k from reciprocal lattice coordinates to
        # Cartesian coordinates.
        # Note that k and prev_k here are row vectors (i.e. 1x3 'dual vectors',
        # the transpose of 3x1 column 'vectors').
        k_Cart = np.dot(k, R)
        prev_k_Cart = np.dot(prev_k, R)
        # Get vector from k to prev_k and its length.
        k_to_prev_k = np.subtract(prev_k_Cart, k_Cart)
        this_dist = np.linalg.norm(k_to_prev_k)
        dists.append(this_dist)

    return dists

def _scaled_k_xs(sym_indices, ks_cut, recip_dists):
    '''Return a list of x values at which each k value in ks_cut will be
    plotted. The sizes of the panels between symmetry points are scaled
    such that these sizes correspond to the Cartesian distance between the
    k-points at the ends of the panel (relative to the sum of these distances
    over all panels).
    '''
    total_dist = sum(recip_dists)
    xs = []

    current_x = 0.0
    base_x = 0.0
    panel_start_index = 0
    panel_number = 0
    for x_index in range(len(ks_cut)):
        # At the last value in a panel, x_index - panel_start_index = points_in_panel - 1.
        # Set step such that (points_in_panel - 1) * step = panel_x_length,
        # i.e. such that the last value in a panel has an x value the appropriate distance
        # away from the first value in that panel.
        points_in_panel = sym_indices[panel_number+1] - sym_indices[panel_number] + 1
        panel_x_length = recip_dists[panel_number] / total_dist
        step = panel_x_length / (points_in_panel - 1)
        current_x = base_x + (x_index - panel_start_index)*step
        xs.append(current_x)

        # When we reach the right side of a panel, we are at the left side of the next panel.
        # In this situation, current_x = the left side of the next panel: make this the new base_x.
        # Also set panel_start_index = x_index so that xs[panel_start_index] = the new base_x,
        # and advance the count of which panel we are on.
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

# Perform approximate equality comparison between vectors.
# Returns True if individual components differ by no more than eps;
# otherwise returns False.
def _vec_approx_equal(u, v, eps):
    if len(u) != len(v):
        raise ValueError("Comparing vectors of unequal length")
    for i in range(len(u)):
        if abs(u[i] - v[i]) > eps:
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
    parser.add_argument('--scf_dir', default=None,
            help="Separate directory for SCF files")
    args = parser.parse_args()

    all_data_paths = None
    scf_data_paths = None
    if args.searchEs:
        all_data_paths = FindEs(args.dir_path)
        if args.scf_dir != None:
            scf_data_paths = FindEs(args.scf_dir)
    else:
        all_data_paths = FindBands(args.dir_path)
        if args.scf_dir != None:
            scf_data_paths = FindEs(args.scf_dir)

    for system_name, system_data_paths in all_data_paths.items():
        E_Fermi, D, magmom = None, None, None
        if args.scf_dir != None:
            try:
                print("Reading OUTCAR file for {}".format(system_name))
                E_Fermi, D = ParseOutcar(scf_data_paths[system_name]['outcar_path'])
                magmom = ParseOszicar(scf_data_paths[system_name]['oszicar_path'])
            except KeyError as e:
                print("Error reading OUTCAR file for {}; skipping.".format(system_name))
                print("Content of the error:")
                print(str(e))
                continue
        else:
            E_Fermi, D = ParseOutcar(system_data_paths['outcar_path'])
            magmom = ParseOszicar(system_data_paths['oszicar_path'])

        try:
            print("Reading EIGENVAL file {}".format(system_data_paths['eigenval_path']))
            ks, eigenvals = ParseEigenval(system_data_paths['eigenval_path'])
        except IndexError as e:
            print("Error reading EIGENVAL file {}; skipping.".format(system_data_paths['eigenval_path']))
            print("Content of the error:")
            print(str(e))
            continue
        
        # The rows of R (the R[i, :]) are the reciprocal lattice vectors.
        R = 2.0 * np.pi * np.linalg.inv(D)

        # TODO - get k_labels from KPOINTS?
        k_labels = None
        if args.structure_type in ["L21", "C1b", "XA"]:
            k_labels = ["$\Gamma$", "$X$", "$W$", "$K$", "$\Gamma$", "$L$", "$U$", "$W$", "$L$", "$K$", "$W$", "$U$", "$X$"]
        elif args.structure_type == "D022":
            k_labels = ["$N$", "$P$", "$X$", "$\Gamma$", "$Z$"]
        else:
            raise ValueError("Structure type {} not supported".format(args.structure_type))

        # TODO - configurable destination (instead of pwd)?
        out_path = "{}_{}".format(system_name, args.structure_type)

        nspin = len(eigenvals[0])
        if nspin == 2:
            eigenvals = swap_channels_if_mag_neg(magmom, eigenvals)

        logo_text = None
        if args.logo:
            logo_text = "www.heusleralloys.mint.ua.edu"

        PlotBands(ks, eigenvals, E_Fermi, k_labels, R, out_path, logo_text)

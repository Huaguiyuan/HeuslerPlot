import numpy as np

def ParseOutcar(outcar_path):
    '''Parse the OUTCAR file at outcar_path and return the Fermi energy and
    lattice vectors contained there. The lattice vectors are given as a numpy
    array D with columns D[:, i] giving the i'th lattice vector.
    '''
    E_Fermi = None
    D = np.zeros((3, 3), dtype=np.float64)

    fp = open(outcar_path, 'r')
    lines = fp.readlines()
    fp.close()

    # Fermi energy line.
    # Looks like:
    #
    #  E-fermi :   8.3134     XC(G=0): -12.9975     alpha+bet :-15.4775
    for line in lines:
        if "E-fermi" in line:
            E_Fermi = float(line.split()[2])
            break

    # Lattice vector lines.
    # Looks like:
    #
    #   direct lattice vectors          reciprocal lattice vectors
    # 0.000000000  2.813095533  2.812809989    -0.177740142  0.177740142  0.177758186
    # 2.813095533  0.000000000  2.812809989     0.177740142 -0.177740142  0.177758186
    # 2.813095533  2.813095533  0.000000000     0.177740142  0.177740142 -0.177758186
    for i, line in enumerate(lines):
        if "direct lattice vectors" in line:
            a_line, b_line, c_line = lines[i+1], lines[i+2], lines[i+3]
            for col, lat_line in enumerate([a_line, b_line, c_line]):
                components = lat_line.split()[0:3]
                for row in range(3):
                    val = components[row]
                    D[row, col] = float(val)

    return E_Fermi, D

def ParseOszicar(oszicar_path):
    '''Parse the OSZICAR file at oszicar_path and return the final value
    of the magnetic moment contained there. If no magnetic moment is
    given, return None.
    '''
    mag = None

    fp = open(oszicar_path, 'r')
    lines = fp.readlines()
    fp.close()

    # Magnetization line.
    # Looks like:
    #
    #   1 F= -.30431484E+02 E0= -.30431484E+02  d E =0.000000E+00  mag=    -4.9999
    for line in lines:
        if "mag=" in line:
            mag = float(line.split()[-1])

    return mag

def ParseEigenval(eigenval_path):
    '''Parse the EIGENVAL file at eigenval_path and return two lists.
    The first list contains the k-points given in the EIGENVAL file, with
    elements of the form (ka, kb, kc). The second list contains the
    corresponding eigenvalues for each spin channel present. If there is
    only one spin channel, the list elements have the form (eigenvals,); if
    both spin channels are present the list elements have the form
    (up_eigenvals, down_eigenvals).
    '''
    ks = []
    all_eigenvals = []

    fp = open(eigenval_path, 'r')
    lines = fp.readlines()
    fp.close()

    # EIGENVAL header.
    # Last value on the first line = number of spins.
    nspin = int(lines[0].strip().split()[-1])
    if nspin not in [1, 2]:
        raise ValueError("Got nspin = {}; expected 1 or 2".format(str(nspin)))
    # Second value on the 6th line = number of k-points.
    # Third value on the 6th line = number of bands.
    nks = int(lines[5].strip().split()[1])
    nbands = int(lines[5].strip().split()[2])

    # Eigenvalues start on the 8th line.
    k_line_index = 7
    while len(ks) < nks:
        # First line of group: ka, kb, kc, weight.
        # Ignore weight.
        k = lines[k_line_index].strip().split()[0:-1]
        ka, kb, kc = list(map(float, k))
        ks.append((ka, kb, kc))
        # Band lines, if nspin = 1: band index, value.
        if nspin == 1:
            k_eigenvals = []
            for band_index in range(nbands):
                band_line = lines[k_line_index + 1 + band_index].strip()
                val = float(band_line.split()[1])
                k_eigenvals.append(val)
            all_eigenvals.append((k_eigenvals,))
        # Band lines, if nspin = 2: band index, up_value, down_value.
        else:
            k_eigenvals_up, k_eigenvals_down = [], []
            for band_index in range(nbands):
                band_line = lines[k_line_index + 1 + band_index].strip()
                val_up = float(band_line.split()[1])
                val_down = float(band_line.split()[2])
                k_eigenvals_up.append(val_up)
                k_eigenvals_down.append(val_down)
            all_eigenvals.append((k_eigenvals_up, k_eigenvals_down))
        # k line + band lines + empty line before next k
        k_line_index += 1 + nbands + 1

    return ks, all_eigenvals

def ParseKpoints(kpoints_path):
    # TODO - get k-point labels from end-of-line comments
    pass

import unittest
import os
import numpy as np
from HeuslerPlot.plot import _cut_duplicates, _recip_dist, _scaled_k_xs
from HeuslerPlot.parseVASP import ParseEigenval, ParseOutcar
from HeuslerPlot.parseVASP_test import _vec_equal

class TestCutDuplicates(unittest.TestCase):
    def test_Co2MnSi(self):
        eigenval_path = os.path.join("test_data", "e29", "Co2MnSi", "BANDS", "EIGENVAL")
        ks, eigenvals = ParseEigenval(eigenval_path)
        sym_indices, ks_cut, eigenvals_cut = _cut_duplicates(ks, eigenvals)
        self.assertEqual(len(sym_indices), 9)
        nk_per_sym = 40
        for i, sym_index in enumerate(sym_indices):
            self.assertEqual(sym_index, i*nk_per_sym - i)
        self.assertEqual(len(ks_cut)-1, (len(sym_indices)-1)*nk_per_sym - (len(sym_indices)-1))
        self.assertEqual(len(eigenvals_cut)-1, (len(sym_indices)-1)*nk_per_sym - (len(sym_indices)-1))

        # TODO - check that ks and eigenvalues in cut lists match
        # corresponding ones in non-cut lists.

class TestScaledKXs(unittest.TestCase):
    def test_Co2MnSi(self):
        eigenval_path = os.path.join("test_data", "e29", "Co2MnSi", "BANDS", "EIGENVAL")
        outcar_path = os.path.join("test_data", "e29", "Co2MnSi", "OUTCAR")
        ks, eigenvals = ParseEigenval(eigenval_path)
        E_Fermi, D = ParseOutcar(outcar_path)
        R = 2.0 * np.pi * np.linalg.inv(D)

        sym_indices, ks_cut, eigenvals_cut = _cut_duplicates(ks, eigenvals)

        recip_dists = _recip_dist(sym_indices, ks_cut, R)
        xs = _scaled_k_xs(sym_indices, ks_cut, recip_dists)

        sym_xs = []
        for i in range(len(sym_indices)):
            sym_xs.append(xs[sym_indices[i]])

        sym_xs_dists = []
        for i, x in enumerate(sym_xs):
            if i != 0:
                #print(x - sym_xs[i-1])
                sym_xs_dists.append(x - sym_xs[i-1])

        # Symmetry points from e29/Co2MnSi/BANDS/KPOINTS.
        # Different from those used by Setyawan (Comp. Mater. Sci. 49, 299 (2010)).
        # Lattice vectors also differ.
        G = (0.0, 0.0, 0.0)
        X = (1/2, 0.0, 0.0)
        W = (1/2, 1/4, 0.0)
        K = (3/8, 3/8, 0.0)
        L = (1/4, 1/4, 1/4)
        U = (1/2, 1/8, 1/8)

        path = [G, X, W, K, G, L, W, U, X]
        dists = []
        for i, k in enumerate(path):
            k_Cart = np.dot(k, R)
            #print(k_Cart)
            if i != 0:
                prev_k_Cart = np.dot(path[i-1], R)
                k_to_prev_k = np.subtract(prev_k_Cart, k_Cart)
                dists.append(np.linalg.norm(k_to_prev_k))

        eps = 1e-9
        for i, d in enumerate(dists):
            #print(d / sum(dists))
            d_norm = d/sum(dists)
            self.assertTrue(abs(d_norm - sym_xs_dists[i]) < eps)

if __name__ == "__main__":
    unittest.main()

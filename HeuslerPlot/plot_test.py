import unittest
import os
import numpy as np
from HeuslerPlot.plot import _cut_duplicates, _recip_dist, _scaled_k_xs
from HeuslerPlot.parseVASP import ParseEigenval, ParseOutcar

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
        debug_print = False

        #eigenval_path = os.path.join("test_data", "e29", "Co2MnSi", "BANDS", "EIGENVAL")
        #outcar_path = os.path.join("test_data", "e29", "Co2MnSi", "OUTCAR")
        eigenval_path = os.path.join("test_data", "e28", "Co2CrSi", "BANDS", "EIGENVAL")
        outcar_path = os.path.join("test_data", "e28", "Co2CrSi", "OUTCAR")
        ks, eigenvals = ParseEigenval(eigenval_path)
        E_Fermi, D = ParseOutcar(outcar_path)
        R = 2.0 * np.pi * np.linalg.inv(D)

        if debug_print:
            for i in range(3):
                print("D[:, i] = {}".format(str(D[:, i])))
            for i in range(3):
                print("R[i, :] = {}".format(str(R[i, :])))

        sym_indices, ks_cut, eigenvals_cut = _cut_duplicates(ks, eigenvals)

        if debug_print:
            print("ks at sym_indices:")
            for sym_i in sym_indices:
                print("Reciprocal = {}".format(str(ks_cut[sym_i])))

        recip_dists = _recip_dist(sym_indices, ks_cut, R)
        xs = _scaled_k_xs(sym_indices, ks_cut, recip_dists)

        sym_xs = []
        for sym_i in sym_indices:
            sym_xs.append(xs[sym_i])

        if debug_print:
            print("Distance between sym_xs:")
        sym_xs_dists = []
        for i, x in enumerate(sym_xs):
            if i != 0:
                if debug_print:
                    print(x - sym_xs[i-1])
                sym_xs_dists.append(x - sym_xs[i-1])

        G = (0.0, 0.0, 0.0)
        X = (1/2, 0.0, 1/2)
        W = (1/2, 1/4, 3/4)
        K = (3/8, 3/8, 3/4)
        L = (1/2, 1/2, 1/2)
        U = (5/8, 1/4, 5/8)

        if debug_print:
            print("Expected ks as sym_indices:")
        path = [G, X, W, K, G, L, W, U, X]
        expected_dists = []
        for i, k in enumerate(path):
            k_Cart = np.dot(k, R)
            if debug_print:
                print("Reciprocal = {}; Cartesian = {}".format(str(k), str(k_Cart)))
            if i != 0:
                prev_k_Cart = np.dot(path[i-1], R)
                k_to_prev_k = np.subtract(prev_k_Cart, k_Cart)
                expected_dists.append(np.linalg.norm(k_to_prev_k))

        if debug_print:
            print("Expected distance between sym_xs:")
        eps = 1e-9
        expected_dists_norm = []
        for i, d in enumerate(expected_dists):
            d_norm = d/sum(expected_dists)
            if debug_print:
                print(d_norm)
            expected_dists_norm.append(d_norm)
            self.assertTrue(abs(d_norm - sym_xs_dists[i]) < eps)

        sym_xs_vals = [0.0]
        expected_xs_vals = [0.0]
        for i, dist in enumerate(sym_xs_dists):
            sym_xs_vals.append(sym_xs_vals[i] + dist)
        for i, expected_dist in enumerate(expected_dists_norm):
            expected_xs_vals.append(expected_xs_vals[i] + expected_dist)
        if debug_print:
            print("Normalized xs values:")
            print(sym_xs_vals)
            print("Expected normalized xs values:")
            print(expected_xs_vals)
        for i, expected_val in enumerate(expected_xs_vals):
            self.assertTrue(abs(expected_val - sym_xs_vals[i]) < eps)

if __name__ == "__main__":
    unittest.main()

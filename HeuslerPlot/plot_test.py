import unittest
from HeuslerPlot.plot import _cut_duplicates
from HeuslerPlot.parseVASP import ParseEigenval
from HeuslerPlot.parseVASP_test import _vec_equal

class TestCutDuplicates(unittest.TestCase):
    def test_Co2MnSi(self):
        ks, eigenvals = ParseEigenval("test_data/e29/Co2MnSi/BANDS/EIGENVAL")
        sym_indices, ks_cut, eigenvals_cut = _cut_duplicates(ks, eigenvals)
        self.assertEqual(len(sym_indices), 9)
        nk_per_sym = 40
        for i, sym_index in enumerate(sym_indices):
            self.assertEqual(sym_index, i*nk_per_sym - i)
        self.assertEqual(len(ks_cut)-1, (len(sym_indices)-1)*nk_per_sym - (len(sym_indices)-1))
        self.assertEqual(len(eigenvals_cut)-1, (len(sym_indices)-1)*nk_per_sym - (len(sym_indices)-1))

        # TODO - check that ks and eigenvalues in cut lists match
        # corresponding ones in non-cut lists.

if __name__ == "__main__":
    unittest.main()

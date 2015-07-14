import unittest
from parseVASP import ParseOutcar, ParseOszicar, ParseEigenval

class TestParseOutcar(unittest.TestCase):
    def test_Co2MnSi(self):
        outcar_path = "test_data/e29/Co2MnSi/OUTCAR"
        E_Fermi, D = ParseOutcar(outcar_path)
        self.assertEqual(E_Fermi, 8.0505)
        _vec_equal(self, D[:, 0], (0.0000000000,  -2.8130955325,   2.8128099892))
        _vec_equal(self, D[:, 1], (0.0000000000,   2.8130955325,   2.8128099892))
        _vec_equal(self, D[:, 2], (-2.8130955325,   0.0000000000,  -2.8128099892))

class TestParseOszicar(unittest.TestCase):
    def test_Co2MnSi(self):
        oszicar_path = "test_data/e29/Co2MnSi/OSZICAR"
        mag = ParseOszicar(oszicar_path)
        self.assertEqual(mag, -4.9999)

class TestParseEigenval(unittest.TestCase):
    def test_Co2MnSi(self):
        eigenval_path = "test_data/e29/Co2MnSi/BANDS/EIGENVAL"
        ks, eigenvals = ParseEigenval(eigenval_path)
        self.assertEqual(len(ks), 320) # nks
        _vec_equal(self, ks[0], (0.0, 0.0, 0.0))
        _vec_equal(self, ks[-1], (0.5, 0.0, 0.0))

        self.assertEqual(len(eigenvals[0]), 2) # nspin
        self.assertEqual(len(eigenvals[0][0]), 80) # nbands

        self.assertEqual(eigenvals[0][0][0], -3.715960) # first k, first eigenvalue
        self.assertEqual(eigenvals[0][1][0], -3.763301)

        self.assertEqual(eigenvals[0][0][-1], 70.829748) # first k, last eigenvalue
        self.assertEqual(eigenvals[0][1][-1], 69.349120)

        self.assertEqual(eigenvals[-1][0][0], -1.436383) # last k, first eigenvalue
        self.assertEqual(eigenvals[-1][1][0], -1.417314)

        self.assertEqual(eigenvals[-1][0][-1], 66.116215) # last k, last eigenvalue
        self.assertEqual(eigenvals[-1][1][-1], 65.895179)


def _vec_equal(testcase, v, u):
    testcase.assertEqual(len(v), len(u))
    for i in range(len(v)):
        testcase.assertEqual(v[i], u[i])

if __name__ == "__main__":
    unittest.main()

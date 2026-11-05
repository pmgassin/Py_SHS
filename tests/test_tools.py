"""compare_outputs.py: identical files agree, a perturbed value is detected."""
import os
import tempfile
import unittest

import common as C
import compare_outputs as CO


class TestCompareOutputs(unittest.TestCase):
    def test_roundtrip(self):
        d = tempfile.mkdtemp(prefix="pyshs_cmp_")
        bf = os.path.join(d, "beta.txt"); C.write_tensor(bf, C.tensor_chiral())
        C.cli("HRS", "polarplot_single", "800\n1.33\n0\n90\n", (bf, "a.txt"), workdir=d)
        text = open(os.path.join(d, "a.txt")).read()
        self.assertEqual(CO.main([os.path.join(d, "a.txt"), os.path.join(d, "a.txt")]), 0)
        C.write_text(os.path.join(d, "b.txt"), text.replace("av=", "av=9", 1))
        self.assertEqual(CO.main([os.path.join(d, "a.txt"), os.path.join(d, "b.txt")]), 1)

    def test_fortran_exponent_and_table(self):
        d = tempfile.mkdtemp(prefix="pyshs_cmp_")
        C.write_text(os.path.join(d, "x.txt"), "av = 1.5D-03\n 0.0 1.0 2.0 3.0\n")
        C.write_text(os.path.join(d, "y.txt"), "av = 1.5E-03\n 0.0 1.0 2.0 3.0\n")
        self.assertEqual(CO.main([os.path.join(d, "x.txt"), os.path.join(d, "y.txt")]), 0)


if __name__ == "__main__":
    unittest.main()

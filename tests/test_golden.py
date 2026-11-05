"""Frozen regression snapshot (see golden_cases.py).  Detects unintended numerical changes."""
import json
import os
import unittest

import numpy as np

import golden_cases as G

PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "golden", "golden.json")


class TestGolden(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        with open(PATH) as f:
            cls.golden = json.load(f)

    def test_cases(self):
        for name in G.CASES:
            with self.subTest(case=name):
                self.assertIn(name, self.golden, "run make_golden.py")
                got = np.array(G.run_case(name))
                ref = np.array(self.golden[name])
                self.assertEqual(got.shape, ref.shape)
                # relative to the largest value (rounding noise ~1e-16 can appear in zero entries)
                tol = 1e-6 * max(np.abs(ref).max(), 1e-300)
                self.assertLess(np.abs(got - ref).max(), tol)


if __name__ == "__main__":
    unittest.main()

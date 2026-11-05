"""Intensity units (point 4) and scattering-angle grid (point 5), through the command line."""
import os
import tempfile
import unittest

import numpy as np

import common as C
from common import PI

# module -> (stdin before the angle(s), needs orientation file)
CASES = {
    "HRS": ("HRS", "800\n1.33\n0\n", False),
    "SHS": ("SHS", "2\n800\n1.33\n0\n", True),
    "sphere": ("sphere_SHS", "1\n800\n1.5\n0\n1.33\n0\n80\n", False),
    "cylinder": ("cylinder_SHS", "800\n1.5\n0\n1.33\n0\n40\n300\n", False),
    "ellipsoid": ("ellipsoid_SHS", "1\n800\n1.5\n0\n1.33\n0\n60\n300\n", False),
    "cube": ("cube_SHS", "1\n800\n1.5\n0\n1.33\n0\n150\n200\n250\n", False),
}
RAW_DEFAULT = ("HRS", "SHS")        # historical default of these modules is `raw`
P = (lambda K: (K * K / (4 * PI)) ** 2)(4 * PI * 1.33 / 800.0)


class _Base(unittest.TestCase):
    def setUp(self):
        self.d = tempfile.mkdtemp(prefix="pyshs_ua_")
        self.bf = os.path.join(self.d, "beta.txt"); C.write_tensor(self.bf, C.tensor_chiral())
        self.of = os.path.join(self.d, "orient.txt"); C.write_text(self.of, C.ORIENT3[:C.ORIENT3.index("\n", 20) + 1])

    def run_mod(self, name, keyword, tail, env=None):
        mod, pre, needs = CASES[name]
        if name == "SHS":
            pre = "2\n800\n1.33\n0\n"
            C.write_text(self.of, C.ORIENT2)
        files = (self.bf, self.of, "out.txt") if needs else (self.bf, "out.txt")
        return C.cli(mod, keyword, pre + tail, files, env, self.d)[0]


class TestUnits(_Base):
    def test_forced_units_scale_by_the_radiative_prefactor(self):
        for name in CASES:
            with self.subTest(module=name):
                base = C.values(self.run_mod(name, "polarplot_single", "90\n"))
                other = "rad" if name in RAW_DEFAULT else "raw"
                forced = C.values(self.run_mod(name, "polarplot_single", "90\n", {"PYSHS_UNITS": other}))
                m = np.abs(base) > 1e-9 * np.abs(base).max()
                expected = P if other == "rad" else 1.0 / P
                self.assertLess(C.rel(forced[m] / base[m], np.full(m.sum(), expected)), 1e-5)

    def test_i_scale_applies_to_every_module(self):
        for name in CASES:
            with self.subTest(module=name):
                base = C.values(self.run_mod(name, "polarplot_single", "90\n"))
                sc = C.values(self.run_mod(name, "polarplot_single", "90\n", {"PYSHS_I_SCALE": "3"}))
                m = np.abs(base) > 1e-9 * np.abs(base).max()
                self.assertLess(C.rel(sc[m] / base[m], np.full(m.sum(), 3.0)), 1e-5)

    def test_convention_is_reported(self):
        out = self.run_mod("HRS", "polarplot_single", "90\n")
        self.assertIn("Units convention: raw", out)
        out = self.run_mod("sphere", "polarplot_single", "90\n")
        self.assertIn("Units convention: rad", out)
        out = self.run_mod("sphere", "polarplot_single", "90\n", {"PYSHS_UNITS": "raw"})
        self.assertIn("Units convention: raw (PYSHS_UNITS)", out)

    def test_hrs_scientific_notation_with_rad_units(self):
        out = self.run_mod("HRS", "polarplot_single", "90\n", {"PYSHS_UNITS": "rad"})
        self.assertGreater(abs(C.values(out, ("ah",))[0]), 0.0)       # would print 0.000000 with %f


class TestAngleGrid(_Base):
    N = 4

    def _rows(self, name, env=None):
        self.run_mod(name, "angle_scattering", "%d\n" % self.N, env)
        return C.read_table(os.path.join(self.d, "out_plot"))

    def test_default_is_n_plus_1_angles_from_0_to_180(self):
        for name in CASES:
            with self.subTest(module=name):
                r = self._rows(name)
                np.testing.assert_allclose(r[:, 0], np.linspace(0, 180, self.N + 1), atol=1e-4)

    def test_rows_equal_single_angle_runs(self):
        for name in CASES:
            with self.subTest(module=name):
                r = self._rows(name)
                for ang, ih0, ih90, iv0, iv90 in r[[1, 2, 3]]:
                    s = C.values(self.run_mod(name, "polarplot_single", "%g\n" % ang), ("ah", "ch", "av", "cv"))
                    got = np.array([ih0, ih90, iv0, iv90])
                    m = np.abs(s) > 1e-9 * np.abs(s).max()
                    self.assertLess(C.rel(got[m], s[m]), 2e-6)

    def test_compatibility_ranges_reproduce_former_outputs(self):
        for name, mode, n_rows, first in (("HRS", "negative", self.N + 1, -180.0),
                                          ("cylinder", "negative", self.N + 1, -180.0),
                                          ("sphere", "full", 2 * self.N + 1, -180.0),
                                          ("cube", "full", 2 * self.N + 1, -180.0)):
            with self.subTest(module=name, mode=mode):
                r = self._rows(name, {"PYSHS_ANGLE_RANGE": mode})
                self.assertEqual(len(r), n_rows)
                self.assertAlmostEqual(r[0, 0], first, places=3)
                d = self._rows(name)                              # default half range
                for row in d:                                     # same values at the same |angle|
                    j = np.argmin(np.abs(np.abs(r[:, 0]) - row[0]))
                    m = np.abs(row[1:]) > 1e-9 * np.abs(d[:, 1:]).max()
                    self.assertLess(C.rel(r[j, 1:][m], row[1:][m]), 1e-6)

    def test_graph_angle_handles_all_ranges(self):
        import subprocess, sys
        for mode in ("half", "negative", "full"):
            self._rows("cylinder", {"PYSHS_ANGLE_RANGE": mode})
            env = dict(os.environ, MPLBACKEND="Agg")
            r = subprocess.run([sys.executable, os.path.join(C.SRC, "graph_angle.py"), "out_plot", "g.png"],
                               cwd=self.d, capture_output=True, text=True, env=env)
            if "No module named 'matplotlib'" in r.stderr:
                self.skipTest("matplotlib not installed")
            self.assertEqual(r.returncode, 0, r.stderr[-400:])
            self.assertTrue(os.path.exists(os.path.join(self.d, "g.png")))


if __name__ == "__main__":
    unittest.main()

"""Minor points of the review: generators, helper scripts, side-file location, warnings,
Gaussian converter, GUI helpers (tested without tkinter)."""
import importlib.util
import os
import re
import subprocess
import sys
import tempfile
import unittest

import numpy as np

import common as C
from common import PI, K_FUND, K_SH, PREF2

import SHS
import ellipsoid_SHS as EL
import script_gauss2PySHS as G2P
import graph_assembly as GA

HAVE_SCIPY = importlib.util.find_spec("scipy") is not None


def run_script(name, args, stdin="", cwd=None):
    env = dict(os.environ, MPLBACKEND="Agg")
    return subprocess.run([sys.executable, os.path.join(C.SRC, name)] + list(args), input=stdin,
                          capture_output=True, text=True, cwd=cwd, env=env)


class TestGenerators(unittest.TestCase):
    def setUp(self):
        self.d = tempfile.mkdtemp(prefix="pyshs_gen_")

    def test_sphere_scripts_write_the_dipole_file_and_do_not_overwrite_it(self):
        for script, stdin in (("script_sphere.py", "50\n40\n"), ("script_sphere_random.py", "50\n40\n0.1\n")):
            out = os.path.join(self.d, script + ".txt")
            r = run_script(script, [out], stdin, self.d)
            self.assertEqual(r.returncode, 0, r.stderr[-400:])
            d = np.loadtxt(out)                      # still a numeric table (not a picture)
            self.assertEqual(d.shape, (40, 6))
            np.testing.assert_allclose(np.linalg.norm(d[:, 3:], axis=1), 50.0, rtol=1e-5)

    def test_cylinder_generator(self):
        out = os.path.join(self.d, "cyl.txt")
        r = run_script("script_cylinder_random.py", [out], "20\n80\n12\n10\n0\n", self.d)
        self.assertEqual(r.returncode, 0, r.stderr[-400:])
        d = np.loadtxt(out)
        self.assertEqual(d.shape, (120, 6))
        np.testing.assert_allclose(np.hypot(d[:, 3], d[:, 4]), 20.0, rtol=1e-5)
        np.testing.assert_allclose(d[:, 1], PI / 2, rtol=1e-6)      # file has 6 decimals

    def test_spheroid_generator_reproduces_the_ellipsoid_model(self):
        """The generator the GUI refers to (it was missing): dipoles with chi = sigma*beta on a
        spheroid must give the absolute result of ellipsoid_SHS."""
        a, c, N = 80.0, 150.0, 4000
        out = os.path.join(self.d, "spheroid.txt")
        r = run_script("script_spheroid_random.py", [out], "%g\n%g\n%d\n0\n" % (a, c, N), self.d)
        self.assertEqual(r.returncode, 0, r.stderr[-400:])
        d = np.loadtxt(out)
        self.assertEqual(d.shape, (N, 6))
        # points are on the surface
        np.testing.assert_allclose((d[:, 3] ** 2 + d[:, 4] ** 2) / a ** 2 + d[:, 5] ** 2 / c ** 2, 1.0, rtol=1e-6)   # file has 6 decimals
        beta = C.tensor_surface()
        u = np.linspace(-1, 1, 200001)
        area = 2 * PI * a * np.trapezoid(np.sqrt(c * c * (1 - u * u) + a * a * u * u), u)
        sigma = N / area
        P = SHS.rotate_tensors(SHS.euler_matrices(d[:, 0], d[:, 1], d[:, 2]), beta).reshape(N, 27)
        th = np.radians([40.0, 90.0, 150.0])
        with C.quiet():
            ds = SHS.compute(P, d[:, 3:6], K_FUND, K_SH, th, 4, False, (30, 15, 30)) / (8 * PI ** 2)
        r_, w_, Pe = EL.build_surface(a, c, beta, 40, 40)
        with C.quiet():
            cV, cH = EL.compute_coefficients(th, K_FUND, K_SH, PREF2, r_, w_, Pe, 30, 30)
        for i in range(3):
            got = np.array([ds[i, 0, 1], ds[i, 1, 0], ds[i, 1, 1], ds[i, 1, 2]]) / sigma ** 2
            ref = np.array([cV[i, 1], cH[i, 0], cH[i, 1], cH[i, 2]]) / PREF2
            self.assertLess(C.rel(got, ref), 1e-2)

    def test_every_generator_of_the_gui_exists(self):
        src = open(os.path.join(C.SRC, "pyshs_gui.py"), encoding="utf-8").read()
        for name in sorted(set(re.findall(r"script_[a-z_0-9]+\.py", src))):
            self.assertTrue(os.path.exists(os.path.join(C.SRC, name)), name + " is missing")


class TestHelperScripts(unittest.TestCase):
    def test_graph_assembly_column_order(self):
        d = tempfile.mkdtemp(prefix="pyshs_ga_")
        p = os.path.join(d, "asm.txt")
        C.write_text(p, "5.5 1.2 0 1 2 3\n0.1 0.2 0 -10 4 7\n")      # phi > 3.2 rad: valid file
        phi, theta, psi, x, y, z = GA.load(p)
        np.testing.assert_allclose(phi, [5.5, 0.1]); np.testing.assert_allclose(x, [1, -10])
        phi, theta, psi, x, y, z = GA.load(p, positions_first=True)
        np.testing.assert_allclose(x, [5.5, 0.1])

    def test_gauss_converter(self):
        v = np.random.default_rng(1).normal(size=18)

        def block(vals, w):
            return ([" Beta(-2w;w,w)   w= %s\n" % w] + [" filler line %d\n" % i for i in range(7)]
                    + [" %d  %.6fD+00\n" % (i + 1, x) for i, x in enumerate(vals)])
        d = tempfile.mkdtemp(prefix="pyshs_g2p_")
        log = os.path.join(d, "job.log")
        C.write_text(log, " text\n" + "".join(block(v, "0.05")) + " more\n" + "".join(block(2 * v, "0.10")))
        out = os.path.join(d, "beta.txt")
        self.assertEqual(G2P.main([log, out]), 0)
        b = np.array([float(l.split()[1]) for l in open(out) if not l.startswith("#")]).reshape(3, 3, 3)
        np.testing.assert_allclose(b, b.transpose(0, 2, 1))                      # symmetric in (j,k)
        np.testing.assert_allclose(b[0, 0, 0], v[0], atol=1e-6)                  # xxx
        np.testing.assert_allclose(b[1, 0, 2], v[10], atol=1e-6)                 # yxz = value 3*3+1
        self.assertEqual(G2P.main([log, out, "--occurrence", "2"]), 0)
        b2 = np.array([float(l.split()[1]) for l in open(out) if not l.startswith("#")]).reshape(3, 3, 3)
        np.testing.assert_allclose(b2, 2 * b, atol=1e-5)
        self.assertEqual(G2P.main([log, out, "--offset", "3"]), 1)                # wrong layout -> error
        self.assertEqual(G2P.main([log, out, "--marker", "Beta(0;0,0)"]), 1)      # absent block -> error


class TestWarningsAndFiles(unittest.TestCase):
    def test_shs_warns_about_surplus_lines_in_the_orientation_file(self):
        d = tempfile.mkdtemp(prefix="pyshs_w_")
        bf = os.path.join(d, "beta.txt"); C.write_tensor(bf, C.tensor_chiral())
        of = os.path.join(d, "orient.txt"); C.write_text(of, C.ORIENT3)         # 3 dipoles
        out, _ = C.cli("SHS", "polarplot_single", "2\n800\n1.33\n0\n90\n", (bf, of, "out.txt"), workdir=d)
        self.assertIn("WARNING: the orientation file contains 1 more line", out)
        out, _ = C.cli("SHS", "polarplot_single", "3\n800\n1.33\n0\n90\n", (bf, of, "out.txt"), workdir=d)
        self.assertNotIn("more line", out)

    def test_side_files_go_next_to_the_output_file_not_in_the_current_directory(self):
        d = tempfile.mkdtemp(prefix="pyshs_side_")
        cwd, res = os.path.join(d, "cwd"), os.path.join(d, "results")
        os.makedirs(cwd); os.makedirs(res)
        bf = os.path.join(d, "beta.txt"); C.write_tensor(bf, C.tensor_chiral())
        of = os.path.join(d, "orient.txt"); C.write_text(of, C.ORIENT2)
        C.cli("HRS", "polarplot_single", "800\n1.33\n0\n90\n", (bf, os.path.join(res, "o1.txt")), workdir=cwd)
        C.cli("SHS", "polarplot_single", "2\n800\n1.33\n0\n90\n", (bf, of, os.path.join(res, "o2.txt")), workdir=cwd)
        C.cli("cylinder_SHS", "polarplot_single", "800\n1.5\n0\n1.33\n0\n50\n300\n90\n",
              (bf, os.path.join(res, "o3.txt")), workdir=cwd)
        self.assertEqual(os.listdir(cwd), [])                                    # nothing left behind
        self.assertTrue(os.path.exists(os.path.join(res, "out_plot")))
        self.assertTrue(os.path.exists(os.path.join(res, "out_arrows_gnuplot")))

    def test_cylinder_report_titles_are_specific(self):
        d = tempfile.mkdtemp(prefix="pyshs_cyl_")
        bf = os.path.join(d, "beta.txt"); C.write_tensor(bf, C.tensor_surface())
        C.cli("cylinder_SHS", "polarplot_integrate", "800\n1.5\n0\n1.33\n0\n50\n300\n60\n120\n3\n", (bf, "o.txt"), workdir=d)
        t = open(os.path.join(d, "o.txt")).read()
        self.assertIn("average over 60 ... 120 deg (3 points)", t)
        self.assertNotIn("CONFIGURATION AT 0", t)

    def test_hrs_has_no_thread_machinery_left_and_is_angle_independent_in_cost(self):
        src = open(os.path.join(C.SRC, "HRS.py"), encoding="utf-8").read()
        self.assertNotIn("ThreadPoolExecutor", src)
        self.assertNotIn("Python threads", src)

    @unittest.skipUnless(HAVE_SCIPY, "scipy not installed")
    def test_mie_warns_about_ignored_tensor_components(self):
        import shs_mie
        self.assertEqual(shs_mie._check_beta_usage({"zzz": 1.0, "zxx": .4, "zyy": .4, "xzx": .3, "xxz": .3,
                                                    "yzy": .3, "yyz": .3}), [])
        msgs = shs_mie._check_beta_usage({"zzz": 1.0, "xyz": .5, "zxx": .4, "zyy": .2})
        text = " ".join(msgs)
        self.assertIn("IGNORED: xyz", text)
        self.assertIn("zxx != zyy", text)


class TestGuiHelpers(unittest.TestCase):
    """The GUI cannot be imported without tkinter: the two helper pieces are extracted from
    the source and exercised on their own."""

    SRC = open(os.path.join(C.SRC, "pyshs_gui.py"), encoding="utf-8").read()

    def test_old_run_directories_are_removed_but_recent_and_foreign_ones_are_kept(self):
        import shutil
        import time
        from pathlib import Path
        m = re.search(r"WORK_KEEP_DAYS = 7.*?\n\n\nUNITS_CHOICES", self.SRC, re.S)
        self.assertIsNotNone(m)
        wd = Path(tempfile.mkdtemp(prefix="pyshs_wd_"))
        ns = {"shutil": shutil, "WORK_DIR": wd}
        exec(m.group(0).rsplit("\n\n\nUNITS_CHOICES", 1)[0], ns)
        old, new, other = wd / "pyshs_old", wd / "pyshs_new", wd / "keep_me"
        for p in (old, new, other):
            p.mkdir()
        t = time.time() - 10 * 86400
        os.utime(old, (t, t)); os.utime(other, (t, t))
        ns["_cleanup_old_runs"]()
        self.assertFalse(old.exists()); self.assertTrue(new.exists()); self.assertTrue(other.exists())

    def test_polar_file_unit_is_decided_once_per_file(self):
        from pathlib import Path
        m = re.search(r"    def _parse_out_plot\(self, path: Path\) -> dict:.*?\n        return result\n", self.SRC, re.S)
        self.assertIsNotNone(m)
        ns = {"np": np, "Path": Path}
        exec("class G:\n" + m.group(0), ns)
        g = ns["G"]()
        d = tempfile.mkdtemp(prefix="pyshs_gp_")
        p = os.path.join(d, "pp.txt")
        C.write_text(p, "\n".join(["gamma IV IH"] + ["%f 1.0 2.0" % (2 * np.pi * i / 100) for i in range(100)]) + "\n")
        gam = np.array(g._parse_out_plot(Path(p))["gamma"])
        self.assertAlmostEqual(gam.max(), 356.4, places=3)             # radians -> degrees
        C.write_text(p, "\n".join(["gamma IV IH"] + ["%f 1.0 2.0" % (3.6 * i) for i in range(100)]) + "\n")
        gam = np.array(g._parse_out_plot(Path(p))["gamma"])
        self.assertAlmostEqual(gam[1], 3.6, places=6)                   # degrees untouched (no row-wise mixing)


if __name__ == "__main__":
    unittest.main()

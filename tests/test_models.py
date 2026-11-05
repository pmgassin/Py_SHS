"""Cross-checks between the physical models (absolute scale, point 4) and between the
modules' handling of beta (point 1).

Equivalence used: a continuum surface susceptibility chi corresponds to a dipole assembly of
surface density sigma [dipoles/nm^2] and molecular beta with chi = sigma * beta, so
  I_SHS(assembly) / sigma^2 = I_continuum / PREF2         (PREF2 = (|K|^2/4pi)^2)
"""
import importlib.util
import os
import tempfile
import unittest

import numpy as np

import common as C
from common import PI, K_FUND, K_SH, PREF2

import SHS
import sphere_SHS as SP
import cylinder_SHS as CY
import ellipsoid_SHS as EL
import cube_SHS as CU

HAVE_SCIPY = importlib.util.find_spec("scipy") is not None
if HAVE_SCIPY:
    import shs_mie as MIE

THETAS = np.radians([40.0, 90.0, 150.0])
THETAS_DEG = (40.0, 90.0, 150.0)


def _shs_assembly(P, pos, grid=None):
    with C.quiet():
        return SHS.compute(P, pos, K_FUND, K_SH, THETAS, 4, False, grid) / (8 * PI ** 2)


class TestAbsoluteScale(unittest.TestCase):
    beta = C.tensor_surface()

    def test_sphere_vs_dipole_shell(self):
        R, N = 100.0, 2500
        k = 2 * PI * C.N_M / C.LAM                 # index-matched sphere: k = K/2
        K = 2 * k
        ph, th, pos = C.fibonacci_sphere(N, R)
        P = SHS.rotate_tensors(SHS.euler_matrices(ph, th, np.zeros(N)), self.beta).reshape(N, 27)
        sigma = N / (4 * PI * R * R)
        with C.quiet():
            d = SHS.compute(P, pos, k, K, THETAS, 4, False) / (8 * PI ** 2)
        for i, t in enumerate(THETAS):
            g = SP.surface_gamma(R, t, self.beta, k, K, 120)
            _, bv, _, ah, bh, ch = SP.coefficients_from_gamma(g, t)
            got = np.array([d[i, 0, 1], d[i, 1, 0], d[i, 1, 1], d[i, 1, 2]]) / sigma ** 2
            self.assertLess(C.rel(got, [bv, ah, bh, ch]), 3e-3)

    @unittest.skipUnless(HAVE_SCIPY, "scipy not installed")
    def test_mie_equals_sphere_at_matched_index(self):
        R = 100.0
        res = MIE.shs_sphere(R, C.LAM, C.N_M, C.N_M, C.N_M, C.N_M,
                             chi_rrr=1.0, chi_rtt=0.4, chi_trt=0.3)
        sv = MIE.shs_vs_gamma(res, theta_scat_deg=90.0, gamma_deg=[0, 90])
        mie = np.array([sv["I_H"][0, 0], sv["I_H"][1, 0]])            # VH (= ah), HH (= ch)
        k = 2 * PI * C.N_M / C.LAM; K = 2 * k
        g = SP.surface_gamma(R, PI / 2, self.beta, k, K, 120)
        _, _, _, ah, _, ch = SP.coefficients_from_gamma(g, PI / 2)
        pref2 = (K * K / (4 * PI)) ** 2
        self.assertLess(C.rel(mie, np.array([ah, ch]) * pref2), 1e-5)

    def test_ellipsoid_sphere_limit(self):
        a = c = 100.0
        r, w, P = EL.build_surface(a, c, self.beta, 40, 40)
        with C.quiet():
            cV, cH = EL.compute_coefficients(THETAS, K_FUND, K_SH, PREF2, r, w, P, 24, 24)
        for i, t in enumerate(THETAS):
            g = SP.surface_gamma(a, t, self.beta, K_FUND, K_SH, 120)
            _, bv, _, ah, bh, ch = SP.coefficients_from_gamma(g, t)
            got = np.array([cV[i, 1], cH[i, 0], cH[i, 1], cH[i, 2]]) / PREF2
            self.assertLess(C.rel(got, [bv, ah, bh, ch]), 3e-3)

    def test_cylinder_vs_dipoles_on_lateral_surface(self):
        R, L, nb, nz = 50.0, 200.0, 48, 64
        b = (np.arange(nb) + 0.5) * 2 * PI / nb
        z = -L / 2 + (np.arange(nz) + 0.5) * L / nz
        B, Z = [x.ravel() for x in np.meshgrid(b, z, indexing="ij")]
        N = B.size
        Tm = SHS.euler_matrices(B, np.full(N, PI / 2), np.zeros(N))
        P = SHS.rotate_tensors(Tm, self.beta).reshape(N, 27)
        pos = np.column_stack([R * np.cos(B), R * np.sin(B), Z])
        sigma = N / (2 * PI * R * L)
        d = _shs_assembly(P, pos, (30, 15, 30))
        with C.quiet():
            acc, wsum = CY.compute(R, L, self.beta, K_FUND, K_SH, THETAS)
        c = acc / wsum
        for i in range(len(THETAS)):
            got = np.array([d[i, 0, 1], d[i, 1, 1], d[i, 1, 2]]) / sigma ** 2
            self.assertLess(C.rel(got, [c[i, 0, 1], c[i, 1, 1], c[i, 1, 2]]), 3e-3)

    def test_cube_vs_dipoles_on_the_six_faces(self):
        a, m = 100.0, 24
        faces = CU.build_faces(a, a, a, self.beta)
        g = (np.arange(m) + 0.5) / m - 0.5
        Pd, pd = [], []
        for f in faces:
            for s in g:
                for t in g:
                    pd.append(f["r0"] + s * f["u"] + t * f["v"]); Pd.append(f["P"].reshape(27))
        Pd, pd = np.array(Pd), np.array(pd)
        sigma = len(Pd) / (6 * a * a)
        d = _shs_assembly(Pd, pd, (30, 15, 30))
        with C.quiet():
            cV, cH = CU.compute_coefficients(THETAS, K_FUND, K_SH, PREF2, faces, 30, 30, 30)
        for i in range(len(THETAS)):
            got = np.array([d[i, 0, 1], d[i, 1, 1], d[i, 1, 2]]) / sigma ** 2
            ref = np.array([cV[i, 1], cH[i, 1], cH[i, 2]]) / PREF2
            self.assertLess(C.rel(got, ref), 3e-3)


class TestBetaSymmetrisation(unittest.TestCase):
    """A tensor that is not symmetric in its last two indices must give exactly the result of
    its symmetrised version, in every module (point 1)."""

    CASES = {
        "HRS": ("HRS", "800\n1.33\n0\n90\n", False),
        "SHS": ("SHS", "2\n800\n1.33\n0\n90\n", True),
        "sphere": ("sphere_SHS", "1\n800\n1.5\n0\n1.33\n0\n150\n90\n", False),
        "cylinder": ("cylinder_SHS", "800\n1.5\n0\n1.33\n0\n50\n300\n90\n", False),
        "ellipsoid": ("ellipsoid_SHS", "1\n800\n1.5\n0\n1.33\n0\n80\n400\n90\n", False),
        "cube": ("cube_SHS", "1\n800\n1.5\n0\n1.33\n0\n200\n200\n200\n90\n", False),
    }

    def test_every_module(self):
        b = C.tensor_asym()
        bs = 0.5 * (b + b.transpose(0, 2, 1))
        with tempfile.TemporaryDirectory() as d:
            fa, fs = os.path.join(d, "asym.txt"), os.path.join(d, "sym.txt")
            C.write_tensor(fa, b); C.write_tensor(fs, bs)
            fo = os.path.join(d, "orient.txt"); C.write_text(fo, C.ORIENT2)
            for name, (mod, stdin, needs_orient) in self.CASES.items():
                with self.subTest(module=name):
                    outs = []
                    for f in (fa, fs):
                        files = (f, fo, "out.txt") if needs_orient else (f, "out.txt")
                        out, _ = C.cli(mod, "polarplot_single", stdin, files, workdir=d)
                        outs.append(out)
                    self.assertIn("symmetrised", outs[0])
                    self.assertNotIn("WARNING: the tensor is not symmetric", outs[1])
                    self.assertLess(C.rel(C.values(outs[0], ("bv", "bh")), C.values(outs[1], ("bv", "bh"))), 1e-10)


if __name__ == "__main__":
    unittest.main()

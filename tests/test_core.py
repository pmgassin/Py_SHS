"""HRS / SHS core: quadrature, analytic values, polarisation formulas, legacy mode,
symmetrisation of beta, Theta -> -Theta symmetry."""
import os
import tempfile
import unittest

import numpy as np

import common as C
from common import PI

import HRS
import SHS


def _fake_single_orientation(T):
    class One:
        pass
    o = One()
    o.n, o.T, o.w = 1, T, np.array([1.0])
    return o


class TestQuadrature(unittest.TestCase):
    def test_euler_matrices_orthonormal(self):
        rng = np.random.default_rng(0)
        T = HRS.euler_matrices(rng.uniform(0, 2 * PI, 50), rng.uniform(0, PI, 50),
                               rng.uniform(0, 2 * PI, 50))
        self.assertLess(np.abs(np.einsum("nij,nkj->nik", T, T) - np.eye(3)).max(), 1e-14)
        self.assertLess(np.abs(np.linalg.det(T) - 1).max(), 1e-14)

    def test_orientation_weights_sum_to_8pi2(self):
        # Gauss-Legendre in cos(theta): the 1/(8 pi^2) normalisation is exact (point 3)
        for mod in (HRS, SHS):
            self.assertAlmostEqual(mod.Orientations().w.sum() / (8 * PI ** 2), 1.0, places=13)

    def test_midpoint_rule_available_for_reproducibility(self):
        os.environ["PYSHS_THETA_RULE"] = "midpoint"
        try:
            s = HRS.Orientations().w.sum() / (8 * PI ** 2)
        finally:
            os.environ.pop("PYSHS_THETA_RULE")
        self.assertAlmostEqual(s, 1.0041242, places=6)      # the former +0.41 % bias


class TestHRSAnalytic(unittest.TestCase):
    def test_beta_zzz_only(self):
        # <beta_XXX^2> = 1/7 , <beta_XYY^2>-type terms 1/35, b = 6/35 (V, gamma-expansion)
        ori = HRS.Orientations()
        b = np.zeros((3, 3, 3)); b[2, 2, 2] = 1
        S = np.zeros((ori.n, 54)); S[:, :27] = b.reshape(27)
        ang = HRS.AngleTable(np.array([PI / 2]), 0, 0)
        r = HRS.finalize(S, ori, ang, 0, False) / (8 * PI ** 2)
        np.testing.assert_allclose(r[0, :3], [1 / 7, 6 / 35, 1 / 35], rtol=1e-10)

    def test_hrs_equals_shs_single_molecule(self):
        with tempfile.TemporaryDirectory() as d:
            bf = os.path.join(d, "beta.txt"); C.write_tensor(bf, C.tensor_chiral())
            of = os.path.join(d, "orient.txt"); C.write_text(of, "0 0 0 0 0 0\n")
            h, _ = C.cli("HRS", "polarplot_single", "800\n1.33\n0\n90\n", (bf, "out.txt"), workdir=d)
            s, _ = C.cli("SHS", "polarplot_single", "1\n800\n1.33\n0\n90\n", (bf, of, "out.txt"), workdir=d)
        np.testing.assert_allclose(C.values(h), C.values(s), rtol=1e-12)


class TestPolarisationFormula(unittest.TestCase):
    """I(gamma) = a c^4 + b c^2 s^2 + c s^4 + b1 c^3 s + b2 c s^3 against a brute-force
    evaluation of the orientational average (independent of the coefficient code)."""

    def test_default_mode_matches_brute_force_for_an_assembly(self):
        rng = np.random.default_rng(3)
        lam, n = 800.0, 1.33
        k = 2 * PI * n / lam; K = 2 * k
        beta = C.tensor_chiral(5)
        N = 6
        ph, th_, ps = rng.uniform(0, 2 * PI, N), np.arccos(rng.uniform(-1, 1, N)), rng.uniform(0, 2 * PI, N)
        pos = rng.uniform(-120, 120, (N, 3))
        P = SHS.rotate_tensors(SHS.euler_matrices(ph, th_, ps), beta).reshape(N, 27)
        THETA = np.radians(70.0)
        with C.quiet():
            r = SHS.compute(P, pos, k, K, np.array([THETA]), 2, False)[0] / (8 * PI ** 2)
        ori = SHS.Orientations(); T, w = ori.T, ori.w
        Lab = np.einsum("nia,njb,nkc,dabc->ndijk", T, T, T, P.reshape(N, 3, 3, 3))
        Ll = np.stack([Lab[:, :, :, 0, 0], Lab[:, :, :, 1, 1], Lab[:, :, :, 0, 1]], axis=-1)
        rl = np.einsum("nia,da->ndi", T, pos)
        ph_ = np.exp(1j * ((2 * k - K * np.cos(THETA)) * rl[:, :, 2] - K * np.sin(THETA) * rl[:, :, 1]))
        A = np.einsum("nd,ndiq->niq", ph_, Ll)
        c, s = np.cos(THETA), np.sin(THETA)
        gam = np.linspace(0, 2 * PI, 13)[:-1]
        cg, sg = np.cos(gam), np.sin(gam)
        E = A[:, :, 0, None] * cg ** 2 + A[:, :, 1, None] * sg ** 2 + 2 * A[:, :, 2, None] * cg * sg
        brute_v = (w[:, None] * np.abs(E[:, 0]) ** 2).sum(0) / (8 * PI ** 2)
        brute_h = (w[:, None] * np.abs(c * E[:, 1] - s * E[:, 2]) ** 2).sum(0) / (8 * PI ** 2)
        for row, brute in ((0, brute_v), (1, brute_h)):
            curve = (r[row, 0] * cg ** 4 + r[row, 1] * cg ** 2 * sg ** 2 + r[row, 2] * sg ** 4
                     + r[row, 3] * cg ** 3 * sg + r[row, 4] * cg * sg ** 3)
            self.assertLess(np.abs(curve - brute).max() / np.abs(brute).max(), 1e-12)

    def test_single_orientation_formula_hrs_and_shs(self):
        """For ONE orientation b1, b2 are non-zero, so every term of the V and H expansions
        is tested against a direct evaluation of |E(gamma)|^2 (sensitive to sign errors in
        the cross terms, which the isotropic average hides)."""
        gam = np.linspace(0, 2 * PI, 13)[:-1]
        cg, sg = np.cos(gam), np.sin(gam)
        for mod in (HRS, SHS):
            T = mod.euler_matrices(np.array([0.7]), np.array([1.1]), np.array([2.3]))
            S = np.zeros((1, 54)); S[:, :27] = C.tensor_chiral(0).reshape(27)
            for thd in (35.0, 110.0):
                ang = mod.AngleTable(np.array([np.radians(thd)]), 0, 0)
                r = mod.finalize(S, _fake_single_orientation(T), ang, 0, False)
                A = mod.rotate_needed(T, S[:, :27].reshape(1, 3, 3, 3))[0]          # [i, q]
                E = A[:, 0, None] * cg ** 2 + A[:, 1, None] * sg ** 2 + 2 * A[:, 2, None] * cg * sg
                c, s = ang.cs[0], ang.sn[0]
                brute = (np.abs(E[0]) ** 2, np.abs(c * E[1] - s * E[2]) ** 2)
                for row in (0, 1):
                    curve = (r[row, 0] * cg ** 4 + r[row, 1] * cg ** 2 * sg ** 2 + r[row, 2] * sg ** 4
                             + r[row, 3] * cg ** 3 * sg + r[row, 4] * cg * sg ** 3)
                    self.assertLess(np.abs(curve - brute[row]).max() / np.abs(brute[row]).max(), 1e-12)

    def test_chiral_cross_terms_vanish_after_isotropic_average(self):
        """b1, b2 are identically zero for an isotropic average (the coefficients of the
        odd terms are the real part of a Hermitian-symmetric contraction, which the chiral
        invariants cannot feed).  Consequence: --legacy cannot change any reported number."""
        lam, n = 800.0, 1.33
        k = 2 * PI * n / lam; K = 2 * k
        worst = 0.0
        for seed in range(4):
            rng = np.random.default_rng(seed)
            N = int(rng.integers(2, 7))
            ph, th_, ps = rng.uniform(0, 2 * PI, N), np.arccos(rng.uniform(-1, 1, N)), rng.uniform(0, 2 * PI, N)
            pos = rng.uniform(-150, 150, (N, 3))
            P = SHS.rotate_tensors(SHS.euler_matrices(ph, th_, ps), C.tensor_chiral(seed)).reshape(N, 27)
            with C.quiet():
                acc = SHS.compute(P, pos, k, K, np.radians([25., 70., 140.]), 2, False) / (8 * PI ** 2)
            worst = max(worst, np.abs(acc[:, :, 3:]).max() / np.abs(acc[:, :, :3]).max())
        self.assertLess(worst, 1e-13)

    def test_cross_terms_are_not_trivially_zero_for_one_orientation(self):
        T = SHS.euler_matrices(np.array([0.7]), np.array([1.1]), np.array([2.3]))
        S = np.zeros((1, 54)); S[:, :27] = C.tensor_chiral(0).reshape(27)
        ang = SHS.AngleTable(np.array([np.radians(70.)]), 0, 0)
        d = SHS.finalize(S, _fake_single_orientation(T), ang, 0, False)
        self.assertGreater(np.abs(d[:, 3:]).min(), 1e-2)


class TestLegacyMode(unittest.TestCase):
    """--legacy reproduces two defects of the old C code (documented in HRS.py / SHS.py):
    (1) bv1, bv2, bh1, bh2 not divided by 8 pi^2 ;
    (2) the sin*cos cross terms of bh1, bh2 paired the wrong components."""

    def _single(self, mod):
        T = mod.euler_matrices(np.array([0.7]), np.array([1.1]), np.array([2.3]))
        S = np.zeros((1, 54)); S[:, :27] = C.tensor_chiral(0).reshape(27)
        ang = mod.AngleTable(np.array([np.radians(70.)]), 0, 0)
        o = _fake_single_orientation(T)
        return (mod.finalize(S, o, ang, 0, False), mod.finalize(S, o, ang, 0, True),
                S[:, :27].reshape(1, 3, 3, 3), T, ang)

    def test_a_b_c_and_v_cross_terms_unchanged(self):
        for mod in (HRS, SHS):
            d, l, *_ = self._single(mod)
            np.testing.assert_allclose(d[:, :3], l[:, :3], rtol=1e-12)
            np.testing.assert_allclose(d[0, 3:], l[0, 3:], rtol=1e-12)

    def test_legacy_h_cross_terms_match_the_documented_mispairing(self):
        for mod in (HRS, SHS):
            d, l, S27, T, ang = self._single(mod)
            A = mod.rotate_needed(T, S27) if mod is HRS else None
            if A is None:      # SHS rotates the (already lab-frame) per-dipole amplitudes
                A = mod.rotate_needed(T, S27)
            c, s = ang.cs[0], ang.sn[0]
            gyx, gzx, gyy, gzy = A[0, 1, 0], A[0, 2, 0], A[0, 1, 1], A[0, 2, 1]
            gyxy, gzxy = A[0, 1, 2], A[0, 2, 2]
            r4 = lambda u, v: 4.0 * (u * np.conj(v)).real
            # correct expansion of 4 Re(Hxx Hxy*)
            correct = (c * c * r4(gyx, gyxy) + s * s * r4(gzx, gzxy)
                       - s * c * (r4(gyx, gzxy) + r4(gzx, gyxy)))
            miswired = (c * c * r4(gyx, gyxy) + s * s * r4(gzx, gzxy)
                        - s * c * (r4(gzx, gzxy) + r4(gyx, gyxy)))
            self.assertAlmostEqual(d[1, 3], correct, places=12)       # default is correct
            self.assertAlmostEqual(l[1, 3], miswired, places=12)      # legacy reproduces the old pairing
            self.assertGreater(abs(correct - miswired), 1e-3)

    def test_legacy_cli_changes_nothing_that_is_reported(self):
        with tempfile.TemporaryDirectory() as d:
            bf = os.path.join(d, "beta.txt"); C.write_tensor(bf, C.tensor_chiral())
            of = os.path.join(d, "orient.txt"); C.write_text(of, C.ORIENT3)
            args = ("3\n800\n1.33\n0\n90\n", (bf, of, "out.txt"))
            a, _ = C.cli("SHS", "polarplot_single", *args, workdir=d)
            b, _ = C.cli("SHS", "polarplot_single", *args, env={"PYSHS_LEGACY": "1"}, workdir=d)
        np.testing.assert_allclose(C.values(a), C.values(b), rtol=1e-12)
        # b1/b2 are pure rounding noise in both modes
        n1 = C.values(a, ("bv1", "bv2", "bh1", "bh2")); n2 = C.values(b, ("bv1", "bv2", "bh1", "bh2"))
        self.assertLess(np.abs(n1).max(), 1e-9)
        self.assertLess(np.abs(n2).max(), 1e-9)


class TestSymmetries(unittest.TestCase):
    def test_theta_minus_theta_symmetry_with_chiral_tensor(self):
        with tempfile.TemporaryDirectory() as d:
            bf = os.path.join(d, "beta.txt"); C.write_tensor(bf, C.tensor_chiral())
            of = os.path.join(d, "orient.txt"); C.write_text(of, C.ORIENT3)
            for t in (35, 110, 160):
                p, _ = C.cli("SHS", "polarplot_single", "3\n800\n1.33\n0\n%d\n" % t, (bf, of, "out.txt"), workdir=d)
                m, _ = C.cli("SHS", "polarplot_single", "3\n800\n1.33\n0\n%d\n" % -t, (bf, of, "out.txt"), workdir=d)
                names = ("av", "cv", "ah", "ch")
                self.assertLess(C.rel(C.values(p, names), C.values(m, names)), 1e-6)


if __name__ == "__main__":
    unittest.main()

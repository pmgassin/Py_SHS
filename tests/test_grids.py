"""Automatic quadrature grids (point 2): the default (auto) grids must reproduce a much
finer reference for large objects, where the former fixed grids failed."""
import os
import unittest

import numpy as np

import common as C
from common import PI, K_FUND, K_SH, PREF2

import SHS
import shs_grid
import sphere_SHS as SP
import cylinder_SHS as CY
import ellipsoid_SHS as EL
import cube_SHS as CU

QMAX = 2 * K_FUND + K_SH


class TestRules(unittest.TestCase):
    def test_grids_are_never_coarser_than_the_historical_ones(self):
        # sphere, dipole assemblies and box: small objects keep the historical sizes exactly
        self.assertEqual(shs_grid.sphere_pas(50, K_FUND, K_SH, 40)[0], 40)
        self.assertEqual(shs_grid.dipole_grid(3.0, K_FUND, K_SH)[:3], (20, 10, 20))
        self.assertEqual(shs_grid.box_grid(20, 20, 20, K_FUND, K_SH)[:3], (24, 24, 24))
        # cylinder / ellipsoid: the rule has a small constant margin, so never coarser
        cyl = shs_grid.cylinder_grid(20, 80, K_FUND, K_SH, (10, 10, 10))[:3]
        self.assertTrue(all(c >= 10 for c in cyl))
        ell = shs_grid.ellipsoid_grid(20, 40, K_FUND, K_SH, (30, 20, 20, 20))[:4]
        self.assertTrue(all(e >= b for e, b in zip(ell, (30, 20, 20, 20))))

    def test_grids_grow_with_size(self):
        self.assertGreater(shs_grid.sphere_pas(2000, K_FUND, K_SH)[0], 40)
        self.assertGreater(shs_grid.cylinder_grid(50, 2000, K_FUND, K_SH)[1], 12)
        self.assertGreater(shs_grid.dipole_grid(600, K_FUND, K_SH)[0], 20)
        self.assertGreater(shs_grid.box_grid(1000, 1000, 1000, K_FUND, K_SH)[0], 24)

    def test_fixed_mode_and_cap_warning(self):
        os.environ["PYSHS_GRID"] = "fixed"
        try:
            self.assertEqual(shs_grid.sphere_pas(4000, K_FUND, K_SH, 40), (40, ""))
        finally:
            os.environ.pop("PYSHS_GRID")
        os.environ["PYSHS_GRID_MAX"] = "60"
        try:
            pas, msg = shs_grid.sphere_pas(4000, K_FUND, K_SH, 40)
        finally:
            os.environ.pop("PYSHS_GRID_MAX")
        self.assertEqual(pas, 60)
        self.assertIn("WARNING", msg)


class TestConvergence(unittest.TestCase):
    beta = C.tensor_surface()

    def test_sphere_large_radius(self):
        for R in (300.0, 1000.0, 2000.0):
            pas, _ = shs_grid.sphere_pas(R, K_FUND, K_SH)
            ref_pas = int(QMAX * R) + 150
            for t in np.radians([90.0, 150.0]):
                a = np.array(SP.coefficients_from_gamma(SP.surface_gamma(R, t, self.beta, K_FUND, K_SH, pas), t))[[1, 3, 4, 5]]
                b = np.array(SP.coefficients_from_gamma(SP.surface_gamma(R, t, self.beta, K_FUND, K_SH, ref_pas), t))[[1, 3, 4, 5]]
                self.assertLess(C.rel(a, b), 1e-5, "R=%g" % R)

    def test_cylinder_long_rod(self):
        R, L = 50.0, 2000.0
        nb, nph, nth, _ = shs_grid.cylinder_grid(R, L, K_FUND, K_SH)[0:3] + (0,)
        th = np.radians([90.0, 150.0])
        with C.quiet():
            a, wa = CY.compute(R, L, self.beta, K_FUND, K_SH, th, nb, nph, nth)
            b, wb = CY.compute(R, L, self.beta, K_FUND, K_SH, th, 64, 90, 90)
        self.assertLess(C.rel((a / wa)[:, :, :3], (b / wb)[:, :, :3]), 1e-4)
        # the former fixed 10x10 grid was off by tens of percent here
        with C.quiet():
            f, wf = CY.compute(R, L, self.beta, K_FUND, K_SH, th, 10, 10, 10)
        self.assertGreater(C.rel((f / wf)[:, :, :3], (b / wb)[:, :, :3]), 0.1)

    def test_ellipsoid_elongated(self):
        a, c = 100.0, 1000.0
        n_al, n_be, n_ph, n_th, _ = shs_grid.ellipsoid_grid(a, c, K_FUND, K_SH)
        th = np.radians([90.0, 150.0])

        def run(na, nb, np_, nt):
            r, w, P = EL.build_surface(a, c, self.beta, na, nb)
            with C.quiet():
                return EL.compute_coefficients(th, K_FUND, K_SH, PREF2, r, w, P, np_, nt)
        v1, h1 = run(n_al, n_be, n_ph, n_th)
        v2, h2 = run(n_al + 40, n_be + 40, n_ph + 30, n_th + 30)
        self.assertLess(C.rel(h1[:, :3], h2[:, :3]), 1e-5)

    def test_cube_large(self):
        a = 1000.0
        faces = CU.build_faces(a, a, a, self.beta)
        n = shs_grid.box_grid(a, a, a, K_FUND, K_SH)[0]
        th = np.radians([90.0, 150.0])
        with C.quiet():
            v1, h1 = CU.compute_coefficients(th, K_FUND, K_SH, 1.0, faces, n, n, n)
            v2, h2 = CU.compute_coefficients(th, K_FUND, K_SH, 1.0, faces, 62, 62, 62)
        self.assertLess(C.rel(h1[:, :3], h2[:, :3]), 2e-3)

    def test_coherent_chain_of_dipoles(self):
        k = 2 * PI * C.N_M / C.LAM; K = 2 * k
        z = np.linspace(-600, 600, 150)
        pos = np.column_stack([0 * z, 0 * z, z])
        N = len(z)
        P = SHS.rotate_tensors(SHS.euler_matrices(0.3 * np.ones(N), 0.6 * np.ones(N), 0.2 * np.ones(N)),
                               self.beta).reshape(N, 27)
        th = np.array([PI / 2])
        n_ph, n_th, n_ps, _ = shs_grid.dipole_grid(600.0, k, K)
        with C.quiet():
            a = SHS.compute(P, pos, k, K, th, 4, False, (n_ph, n_th, n_ps))
            b = SHS.compute(P, pos, k, K, th, 4, False, (80, 40, 80))
            f = SHS.compute(P, pos, k, K, th, 4, False, (20, 10, 20))
        self.assertLess(C.rel(a[:, :, :3], b[:, :, :3]), 1e-4)
        self.assertGreater(C.rel(f[:, :, :3], b[:, :, :3]), 0.1)     # the fixed grid was wrong by ~30 %


if __name__ == "__main__":
    unittest.main()

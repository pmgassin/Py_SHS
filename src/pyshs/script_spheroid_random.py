"""Dipoles on an ellipsoid of revolution (spheroid), for SHS.py.

Usage:  python script_spheroid_random.py <output_file>      (answers are read from stdin)

Questions (same order as in the GUI): equatorial radius a (nm), polar semi-axis c (nm)
(the symmetry axis is z, as in ellipsoid_SHS.py), number of dipoles, disorder sigma (rad).

Output (one line per dipole, the format read by SHS.py):   phi theta psi x y z
  * the points are quasi-uniformly distributed in *area* on the surface (golden-angle spiral
    with a latitude spacing taken from the cumulative area, so the surface density sigma is
    constant -- needed to compare with ellipsoid_SHS.py through chi = sigma * beta);
  * each dipole frame is the one used by ellipsoid_SHS.py: z = outward normal,
    x = meridian tangent (towards increasing polar angle), y = azimuthal tangent.  With
    psi = 0 and (phi, theta) = spherical angles of the normal this is exactly Euler(phi, theta, 0);
  * disorder: Gaussian noise of standard deviation sigma added to phi and theta.
    (A noise on phi is a real tilt of sigma*sin(theta) only, so the angular spread is larger
    near the poles of the normal than near the equator.)
A 3D picture <output_file>.png is written when matplotlib is available.
"""
import sys

import numpy as np


def spheroid_points(a, c, n, rng=None, randomize=True):
    """Return phi, theta, psi, x, y, z for n points of constant area density."""
    rng = rng or np.random.default_rng()
    # cumulative area as a function of u = cos(alpha):  dS = 2 pi a N(u) du
    u = np.linspace(-1.0, 1.0, 40001)
    nf = np.sqrt(c * c * (1.0 - u * u) + a * a * u * u)
    cum = np.concatenate([[0.0], np.cumsum(0.5 * (nf[1:] + nf[:-1]) * np.diff(u))])
    cum /= cum[-1]
    i = np.arange(n)
    shift = rng.random() if randomize else 0.0
    ui = np.interp((i + 0.5) / n, cum, u)                    # quasi-uniform in area
    beta = ((i + shift * n) % n) * np.pi * (3.0 - np.sqrt(5.0))
    sa = np.sqrt(1.0 - ui * ui)
    x, y, z = a * sa * np.cos(beta), a * sa * np.sin(beta), c * ui
    nfac = np.sqrt(c * c * sa * sa + a * a * ui * ui)
    # outward normal ~ (c sa cos(beta), c sa sin(beta), a u)
    theta = np.arccos(a * ui / nfac)
    phi = np.mod(beta, 2.0 * np.pi)
    return phi, theta, np.zeros(n), x, y, z


def main():
    if len(sys.argv) < 2:
        print("usage: python script_spheroid_random.py <output_file>")
        return 1
    out = sys.argv[1]
    a = float(input("equatorial radius a of the spheroid in nm? "))
    c = float(input("polar semi-axis c (symmetry axis z) in nm? "))
    n = int(input("number of dipoles onto the spheroid? "))
    sigma = float(input("disorder index, sigma of the Gaussian distribution (rad)? "))
    rng = np.random.default_rng()
    phi, theta, psi, x, y, z = spheroid_points(a, c, n, rng)
    phi_o = phi + rng.normal(size=n) * sigma
    theta_o = theta + rng.normal(size=n) * sigma
    with open(out, "w") as f:
        for k in range(n):
            f.write("%12.6f %12.6f %12.6f %12.6f %12.6f %12.6f\n"
                    % (phi_o[k], theta_o[k], psi[k], x[k], y[k], z[k]))
    try:
        import matplotlib
        if not matplotlib.get_backend():
            matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig = plt.figure()
        ax = fig.add_subplot(111, projection="3d")
        s = max(a, c) / 4.0
        ax.quiver(x, y, z, s * np.cos(phi_o) * np.sin(theta_o), s * np.sin(phi_o) * np.sin(theta_o),
                  s * np.cos(theta_o), color="b")
        plt.savefig(out if out.lower().endswith((".png", ".pdf", ".jpg")) else out + ".png")
    except Exception as exc:                     # the plot is optional
        print("(no 3D picture: %s)" % exc)
    return 0


if __name__ == "__main__":
    sys.exit(main())

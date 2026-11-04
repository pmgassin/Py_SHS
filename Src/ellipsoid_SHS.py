#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ellipsoid_SHS.py  --  Py_SHS : Coherent SHS by an ellipsoid of revolution (Python / NumPy)

Python version of ellipsoid_SHS.c.  Same physics, same input protocol,
same output files.

Physical model
--------------
Spheroid (ellipsoid of revolution) with equatorial radius a and half-length c
along the symmetry axis (c = L/2).

  x = a sin(alpha) cos(beta)
  y = a sin(alpha) sin(beta)       alpha in [0, pi], beta in [0, 2 pi]
  z = c cos(alpha)

Local surface frame of the ki(2) tensor (same convention as the cylinder):
  local x = meridian tangent
  local y = azimuthal tangent
  local z = outward surface normal

Surface integration is coherent (amplitudes); orientation average is
incoherent (intensities), normalised by the sum of orientation weights,
then multiplied by the radiative prefactor (|K_2ω|²/4π)².

Keywords : polarplot_single | polarplot_integrate | angle_scattering
Usage    : python ellipsoid_SHS.py keyword ki2_input_file output_file
stdin    : NLWKB flag, wavelength, refractive index(es), a, length, then angle(s)

Environment
-----------
PYSHS_I_SCALE=x   optional intensity scale factor (default 1).
"""
import os
import sys
import math

import numpy as np

import shs_grid
import shs_units
import shs_angles
import shs_io

PI = math.pi

# ---- numerical sampling (same defaults as the C code) ----
# Minimum / historical grid sizes; enlarged automatically for large particles
# (see shs_grid.py; PYSHS_GRID=fixed keeps these values).
N_ALPHA = 30       # surface points along alpha (minimum)
N_BETA = 20        # surface points along beta (minimum)
N_PHI = 20         # orientations: phi in [0, 2pi) (minimum)
N_THETA = 20       # orientations: theta in [0, pi] (minimum)
GRID_USED = [N_ALPHA, N_BETA, N_PHI, N_THETA]   # grid actually used (set in main)


# --------------------------------------------------------------------------
# Helpers
# --------------------------------------------------------------------------
class TokenReader:
    """Reads whitespace-separated values from stdin like C scanf."""

    def __init__(self):
        self.buf = []

    def next(self):
        while not self.buf:
            line = sys.stdin.readline()
            if line == "":
                raise EOFError("unexpected end of the standard input")
            self.buf = line.split()[::-1]
        return self.buf.pop()

    def real(self):
        return float(self.next())

    def integer(self):
        return int(self.next())


def say(text="", end="\n"):
    print(text, end=end, flush=True)


def euler_zyz(phi, theta, psi):
    """Single ZYZ Euler matrix (3, 3). lab = T @ particle."""
    cp, sp = math.cos(phi), math.sin(phi)
    ct, st = math.cos(theta), math.sin(theta)
    cs, ss = math.cos(psi), math.sin(psi)
    T = np.empty((3, 3))
    T[0, 0] = -sp * ss + ct * cp * cs
    T[0, 1] = -sp * cs - ct * cp * ss
    T[0, 2] = st * cp
    T[1, 0] = cp * ss + ct * sp * cs
    T[1, 1] = cp * cs - ct * sp * ss
    T[1, 2] = st * sp
    T[2, 0] = -st * cs
    T[2, 1] = st * ss
    T[2, 2] = ct
    return T


def rotate_tensor(M, beta):
    """out[i,j,k] = sum_abc M[i,a] M[j,b] M[k,c] beta[a,b,c]."""
    t1 = np.einsum("ia,abc->ibc", M, beta)
    t2 = np.einsum("jb,ibc->ijc", M, t1)
    return np.einsum("kc,ijc->ijk", M, t2)


def read_ki(name):
    """1 header line + 27 lines 'name value'."""
    ki = np.zeros((3, 3, 3))
    try:
        with open(name, "rt") as f:
            f.readline()
            for a in range(3):
                for b in range(3):
                    for c in range(3):
                        ki[a, b, c] = float(f.readline().split()[1])
    except (OSError, IndexError, ValueError):
        return None
    # Intensity formulas assume beta_ijk = beta_ikj: symmetrise in the last two indices
    # so that every module treats any input tensor the same way.
    asym = np.abs(ki - ki.transpose(0, 2, 1)).max()
    if asym > 1e-12 * max(np.abs(ki).max(), 1e-300):
        say("WARNING: the tensor is not symmetric in its last two indices "
            "(max |beta_ijk - beta_ikj| = %.3e): it has been symmetrised, "
            "beta_ijk = beta_ikj = (beta_ijk + beta_ikj)/2." % asym)
        ki = 0.5 * (ki + ki.transpose(0, 2, 1))
    return ki


# --------------------------------------------------------------------------
# Surface grid
# --------------------------------------------------------------------------
def build_surface(a, c, ki, n_alpha=None, n_beta=None):
    """Build surface points: positions, weights, ki(2) in the particle frame.

    Gauss-Legendre in u = cos(alpha) (the area element a*N(u) du dbeta is smooth in u)
    and periodic trapezoid in beta.  PYSHS_THETA_RULE=midpoint restores the former
    mid-point rule in alpha.

    Returns arrays:
      r   (NPTS, 3)   positions
      w   (NPTS,)     surface element weights
      P   (NPTS, 3, 3, 3)  ki(2) in particle frame
    """
    n_alpha = N_ALPHA if n_alpha is None else n_alpha
    n_beta = N_BETA if n_beta is None else n_beta
    hbe = 2.0 * PI / n_beta

    if os.environ.get("PYSHS_THETA_RULE", "gauss").strip().lower() == "midpoint":
        al = (np.arange(n_alpha) + 0.5) * (PI / n_alpha)
        ca, sa = np.cos(al), np.sin(al)
        wal = np.full(n_alpha, PI / n_alpha) * sa            # weight in alpha: sin(al) dal
    else:
        ca, wal = np.polynomial.legendre.leggauss(n_alpha)   # u = cos(alpha), weight du
        sa = np.sqrt(1.0 - ca * ca)
    be = (np.arange(n_beta) + 0.5) * hbe
    sb, cb = np.sin(be), np.cos(be)

    # arrays of shape (n_alpha, n_beta)
    SAg = sa[:, None] * np.ones((1, n_beta))
    CAg = ca[:, None] * np.ones((1, n_beta))
    CBg = np.ones((n_alpha, 1)) * cb[None, :]
    SBg = np.ones((n_alpha, 1)) * sb[None, :]
    Nfac = np.sqrt(c * c * SAg ** 2 + a * a * CAg ** 2)

    r = np.stack([a * SAg * CBg, a * SAg * SBg, c * CAg], axis=-1).reshape(-1, 3)
    # area element: a sin(al) N(al) dal dbeta = a N(u) du dbeta   (u = cos al)
    w = (a * Nfac * wal[:, None] * hbe).reshape(-1)
    t = np.stack([a * CAg * CBg / Nfac, a * CAg * SBg / Nfac, -c * SAg / Nfac], axis=-1)
    ph = np.stack([-SBg, CBg, np.zeros_like(SBg)], axis=-1)
    n = np.stack([c * SAg * CBg / Nfac, c * SAg * SBg / Nfac, a * CAg / Nfac], axis=-1)
    S = np.stack([t, ph, n], axis=-1).reshape(-1, 3, 3)     # columns = local axes
    P = np.einsum("nia,njb,nkc,abc->nijk", S, S, S, ki)      # local -> particle frame
    return r, w, P


# --------------------------------------------------------------------------
# Core calculation
# --------------------------------------------------------------------------
def theta_rule(n):
    """Quadrature nodes/weights for  int_0^pi f(theta) sin(theta) dtheta.

    Default: Gauss-Legendre in cos(theta) (exact for polynomials of degree <= 2n-1,
    weights sum to exactly 2).  PYSHS_THETA_RULE=midpoint restores the former midpoint
    rule (weights h*sin(theta), sum slightly different from 2).
    """
    if os.environ.get("PYSHS_THETA_RULE", "gauss").strip().lower() == "midpoint":
        h = PI / n
        th = (np.arange(n) + 0.5) * h
        return th, h * np.sin(th)
    x, wl = np.polynomial.legendre.leggauss(n)
    return np.arccos(x), wl



def compute_coefficients(Theta, petitk, grandk, scale, r, w, P, n_phi=None, n_theta=None):
    """Orientation average of SHS intensity for a list of scattering angles.

    Orientations: n_phi trapezoid points in phi x n_theta Gauss-Legendre points in
    cos(theta) (the particle is axisymmetric, no psi).
    Returns coefV[nang, 5], coefH[nang, 5] already normalised and scaled.
    """
    n_phi = N_PHI if n_phi is None else n_phi
    n_theta = N_THETA if n_theta is None else n_theta
    nang = len(Theta)
    Theta = np.asarray(Theta, dtype=np.float64)
    dphi = 2.0 * PI / n_phi
    th_nodes, th_w = theta_rule(n_theta)

    coefV = np.zeros((nang, 5))
    coefH = np.zeros((nang, 5))
    wsum = 0.0

    # Precompute angle-dependent q
    cs = np.cos(Theta)
    sn = np.sin(Theta)
    qz = 2.0 * petitk - grandk * cs
    qy = -grandk * sn

    def add_pol(X, weight, coef):          # X: (nang, 3)  xx, yy, xy
        coef[:, 0] += weight * np.abs(X[:, 0]) ** 2
        coef[:, 1] += weight * (4.0 * np.abs(X[:, 2]) ** 2
                                + 2.0 * np.real(X[:, 0] * np.conj(X[:, 1])))
        coef[:, 2] += weight * np.abs(X[:, 1]) ** 2
        coef[:, 3] += weight * 4.0 * np.real(X[:, 0] * np.conj(X[:, 2]))
        coef[:, 4] += weight * 4.0 * np.real(X[:, 1] * np.conj(X[:, 2]))

    for ip in range(n_phi):
        say("calculation in progress : %5.1f pourcent " % (100.0 * ip / n_phi))
        phi = (ip + 0.5) * dphi
        for it in range(n_theta):
            th = th_nodes[it]
            wori = dphi * th_w[it]
            wsum += wori

            T = euler_zyz(phi, th, 0.0)

            # Lab-frame tensors (xx, yy, xy components) and positions, all points at once.
            # M[p,a] = sum_bc P[p,a,b,c] T[x,b] T[y,c]  ->  full[p,i] = sum_a T[i,a] M[p,a]
            X, Y = T[0], T[1]
            Mxx = np.einsum("pabc,b,c->pa", P, X, X)
            Myy = np.einsum("pabc,b,c->pa", P, Y, Y)
            Mxy = 0.5 * (np.einsum("pabc,b,c->pa", P, X, Y)
                         + np.einsum("pabc,b,c->pa", P, Y, X))
            L = np.empty((len(w), 3, 3))      # [p, i, q]  q=xx,yy,xy
            L[:, :, 0] = w[:, None] * (Mxx @ T.T)
            L[:, :, 1] = w[:, None] * (Myy @ T.T)
            L[:, :, 2] = w[:, None] * (Mxy @ T.T)
            yl = r @ T[1]
            zl = r @ T[2]

            ph = np.exp(1j * (qz[:, None] * zl[None, :] + qy[:, None] * yl[None, :]))
            A = np.einsum("ap,piq->aiq", ph, L)                 # (nang, 3, 3)

            add_pol(A[:, 0, :], wori, coefV)                    # V: lab x
            Hs = cs[:, None] * A[:, 1, :] - sn[:, None] * A[:, 2, :]
            add_pol(Hs, wori, coefH)

    # Orientation average + radiative prefactor
    coefV *= scale / wsum
    coefH *= scale / wsum
    return coefV, coefH


# --------------------------------------------------------------------------
# Outputs
# --------------------------------------------------------------------------
def write_header(fout, use_nlwkb, lam, n_p, n_m, a, length, petitk, grandk,
                 rad_pref, rad_pref2, ki):
    fout.write("********************************************************\n")
    fout.write("****           Py_SHS: An open source software      ****\n")
    fout.write("****         about Second Harmonic Scattering       ****\n")
    fout.write("**** Dr Pierre-Marie GASSIN and Dr Gaelle GASSIN    ****\n")
    fout.write("****                  ICGM - ENSCM                  ****\n")
    fout.write("****                   june 2021                    ****\n")
    fout.write("********************************************************\n")
    fout.write("                                                        \n")
    fout.write("           Don't forget to cite this work:              \n")
    fout.write("    J. Chem. Inf. Model. 2020, 60, 12, pp 5912–5917     \n")
    fout.write("                                                        \n")
    fout.write("********************************************************\n")
    fout.write("********************************************************\n")
    fout.write("**       SIMULATION RESULTS  for ellipsoid programm ****\n")
    fout.write("********************************************************\n")
    fout.write("  Mode                : %s\n" % ("NLWKB" if use_nlwkb else "Classic RGD"))
    fout.write("  wavelength = %f nm\n" % lam)
    if use_nlwkb:
        fout.write("  n_p = %.5f + i*%.5f   n_m = %.5f + i*%.5f\n"
                   % (n_p.real, n_p.imag, n_m.real, n_m.imag))
    else:
        fout.write("  n = %.5f + i*%.5f\n" % (n_p.real, n_p.imag))
    fout.write("  equatorial radius a = %f nm   length 2c = %f nm\n" % (a, length))
    fout.write("  petitk  (k_omega)   = %.6f + i*%.6f\n" % (petitk.real, petitk.imag))
    fout.write("  grandk  (k_2omega)  = %.6f + i*%.6f\n" % (grandk.real, grandk.imag))
    fout.write("  Radiative pref.     : |k|^2/(4pi) = %.6e  (I *= pref^2 = %.6e)\n"
               % (rad_pref, rad_pref2))
    fout.write("  Quadrature          : surface %d x %d, orientations %d x %d\n"
               % (GRID_USED[0], GRID_USED[1], GRID_USED[2], GRID_USED[3]))
    fout.write("*********          SURFACE KI(2) TENSOR      ***********\n")
    fout.write("\n")
    labels = "xyz"
    for a1 in range(3):
        for a2 in range(3):
            for a3 in range(3):
                fout.write("    %c%c%c=%f" % (labels[a1], labels[a2], labels[a3],
                                              ki[a1, a2, a3]))
            fout.write("\n")
        fout.write("\n")
    fout.write("************************************************************************\n")


def polar_output(fout, title, V, H):
    av, bv, cv, bv1, bv2 = V
    ah, bh, ch, bh1, bh2 = H
    dv = 3 * av + bv + 3 * cv
    dh = 3 * ah + bh + 3 * ch
    with np.errstate(all="ignore"):
        I4v = (av - bv + cv) / dv if dv != 0 else 0.0
        I2v = 4 * (av - cv) / dv if dv != 0 else 0.0
        I4h = (ah - bh + ch) / dh if dh != 0 else 0.0
        I2h = 4 * (ah - ch) / dh if dh != 0 else 0.0

    say("")
    say("**********calculation complete ************")
    say("av = %e" % av)
    say("bv = %e" % bv)
    say("bv1 = %e" % bv1)
    say("bv2 = %e" % bv2)
    say("cv = %e" % cv)
    say("ah = %e" % ah)
    say("bh = %e" % bh)
    say("bh1 = %e" % bh1)
    say("bh2 = %e" % bh2)
    say("ch = %e" % ch)
    say("********************************************")
    say("I2v = %f " % I2v)
    say("I4v = %f " % I4v)
    say("I2h = %f " % I2h)
    say("I4h = %f " % I4h)

    fout.write("********************************************\n")
    fout.write("                                            \n")
    fout.write("* %s *\n" % title)
    fout.write("                                            \n")
    fout.write("********************************************\n")
    fout.write("*************************\n")
    fout.write("I2v=%f  I4v=%f \n" % (I2v, I4v))
    fout.write("I2h=%f  I4h=%f \n" % (I2h, I4h))
    fout.write("*************************\n")
    fout.write("av=%e   bv=%e   bv1=%e  bv2=%e  cv= %e \n"
               % (av, bv, bv1, bv2, cv))
    fout.write("ah=%e   bh=%e   bh1=%e  bh2=%e  ch= %e \n"
               % (ah, bh, bh1, bh2, ch))
    fout.close()

    try:
        fp = open(shs_io.path("out_plot"), "w", encoding="utf-8")
    except OSError:
        say("impossible d'ouvrir le fichier pour ecrire gnuplot!")
        return
    with fp:
        fp.write("gamma IV(gamma) IH(gamma)\n")
        for iii in range(100):
            g = 2 * PI * iii / 100
            nc, ns = math.cos(g), math.sin(g)
            Iv = (av * nc**4 + bv * nc**2 * ns**2 + cv * ns**4
                  + bv1 * nc**3 * ns + bv2 * nc * ns**3)
            Ih = (ah * nc**4 + bh * nc**2 * ns**2 + ch * ns**4
                  + bh1 * nc**3 * ns + bh2 * nc * ns**3)
            fp.write("%.8e %.8e %.8e\n" % (g, Iv, Ih))
    say("-----------------------------------------------------------------------------")
    say("calculation finished, to plot the results, open the outplot file with gnuplot")
    say("-----------------------------------------------------------------------------")


# --------------------------------------------------------------------------
def main(argv):
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except AttributeError:
        pass

    if len(argv) < 4:
        say("Computation failed: ")
        say("3 arguments expected : the keyword of the calculation, "
            "the Ki(2) inputfile and the outputfile. ")
        say("USE ONLY polarplot_single OR polarplot_integrate OR angle_scattering !!")
        say("For more details see the documentation ")
        return 1

    keyword = argv[1]
    mode_single = keyword == "polarplot_single"
    mode_integ = keyword == "polarplot_integrate"
    mode_angle = keyword == "angle_scattering"
    if not (mode_single or mode_integ or mode_angle):
        say("--------------------- --------------------------------------------------")
        say("             ERROR in the name of the keyword                           ")
        say("USE ONLY polarplot_single OR polarplot_integrate OR angle_scattering !!  ")
        say(" ! no calculation has been done ! Retry with the correct keyword        ")
        say("------------------------------------------------------------------------")
        return 1

    if len(argv) > 4:
        say("Too many arguments supplied (extra arguments ignored).")

    say("********************************************************")
    say("****           Py_SHS: An open source software      ****")
    say("****         about Second Harmonic Scattering       ****")
    say("**** Dr Pierre-Marie GASSIN and Dr Gaelle GASSIN    ****")
    say("****                  ICGM - ENSCM                  ****")
    say("****                   june 2021                    ****")
    say("********************************************************")
    say("      Ellipsoid SHS module with optional NLWKB          ")
    say("           Don't forget to cite this work:              ")
    say("    J. Chem. Inf. Model. 2020, 60, 12, pp 5912–5917     ")
    say("                                                        ")
    say("********************************************************")
    say("********       Enter the input parameters     **********")
    say("********************************************************")

    rd = TokenReader()
    try:
        say("Do you want to use NLWKB approximation? (1 = yes, 0 = no [classic RGD]): ",
            end="")
        use_nlwkb = rd.integer()
        say("")

        say("Enter the wavelength in nm: ", end="")
        lam = rd.real()

        if use_nlwkb:
            say("")
            say("=== NLWKB mode activated ===")
            say("Enter REAL part of particle refractive index n_p: ", end="")
            nr = rd.real()
            say("Enter IMAGINARY part of particle refractive index n_p: ", end="")
            ni = rd.real()
            n_p = complex(nr, ni)
            say("Enter REAL part of medium refractive index n_m: ", end="")
            nr = rd.real()
            say("Enter IMAGINARY part of medium refractive index n_m: ", end="")
            ni = rd.real()
            n_m = complex(nr, ni)
            say("")
            say("n_p (particle) = %.5f + i*%.5f" % (n_p.real, n_p.imag))
            say("n_m (medium)   = %.5f + i*%.5f" % (n_m.real, n_m.imag))
        else:
            say("")
            say("=== Classic RGD mode (single index) ===")
            say("Enter REAL part of the refractive index n: ", end="")
            nr = rd.real()
            say("Enter IMAGINARY part of the refractive index n: ", end="")
            ni = rd.real()
            n_p = complex(nr, ni)
            n_m = n_p
            say("")
            say("n = %.5f + i*%.5f" % (n_p.real, n_p.imag))

        say("Enter the equatorial radius of the ellipsoid in nm: ", end="")
        a = rd.real()
        say("")
        say("Enter the length of the ellipsoid along its symmetry axis in nm (=2c)? ",
            end="")
        length = rd.real()
        say("      ******************************     ")
    except (EOFError, ValueError):
        return 1

    if a <= 0 or length <= 0 or lam <= 0:
        sys.stderr.write("wavelength, radius and length must be strictly positive.\n")
        return 1
    c = length / 2.0

    petitk = 2.0 * PI * n_p / lam
    grandk = 4.0 * PI * n_m / lam
    if use_nlwkb:
        say("NLWKB vectors:")
    else:
        say("Classic RGD vectors:")
    say("petitk (k_omega) = %.6f + i*%.6f" % (petitk.real, petitk.imag))
    say("grandk (k_2omega)= %.6f + i*%.6f" % (grandk.real, grandk.imag))
    say("")

    unit = shs_units.settings(grandk, "rad")
    rad_pref = unit.pref
    rad_pref2 = unit.multiplier
    for _l in unit.lines():
        say(_l)
    say("")
    say("The input ki(2) file is %s,  and the output file is %s "
        % (argv[2], argv[3]))

    ki = read_ki(argv[2])
    if ki is None:
        say("Problem to open or to read the input ki(2) file! \n")
        return 1

    shs_io.set_output(argv[3])

    try:
        fout = open(argv[3], "w", encoding="utf-8")
    except OSError:
        say("problem to open the output file to write results!")
        return 1

    n_al, n_be, n_ph, n_th, grid_msg = shs_grid.ellipsoid_grid(
        a, c, petitk, grandk, (N_ALPHA, N_BETA, N_PHI, N_THETA))
    GRID_USED[:] = [n_al, n_be, n_ph, n_th]
    say("Grid: surface %d x %d, orientations %d x %d" % (n_al, n_be, n_ph, n_th))
    if grid_msg:
        say(grid_msg)

    write_header(fout, use_nlwkb, lam, n_p, n_m, a, length, petitk, grandk,
                 rad_pref, rad_pref2, ki)
    shs_units.write_header_lines(fout, unit)

    r, w, P = build_surface(a, c, ki, n_al, n_be)

    # ---- modes ----
    if mode_single:
        try:
            say("Enter the scattered angle in degree (0 = transmission): ", end="")
            gdeg = rd.real()
        except (EOFError, ValueError):
            return 1
        Theta = np.array([gdeg * 2 * PI / 360])
        say("********** Begin calculation for ellipsoid polar plot at %.1f ° ********"
            % gdeg)
        cV, cH = compute_coefficients(Theta, petitk, grandk, rad_pref2, r, w, P, n_ph, n_th)
        title = "Ellipsoid POLAR PLOT - CONFIGURATION AT %.1f deg" % gdeg
        polar_output(fout, title, cV[0], cH[0])

    elif mode_integ:
        try:
            say("Enter the beginning scattered angle in degree (0 = transmission): ",
                end="")
            gdeg1 = rd.real()
            say("Enter the final scattered angle in degree (0 = transmission): ",
                end="")
            gdeg2 = rd.real()
            say("Number of integrale points: ", end="")
            point = rd.integer()
            if point < 1:
                point = 1
        except (EOFError, ValueError):
            return 1

        th1 = gdeg1 * 2 * PI / 360
        th2 = gdeg2 * 2 * PI / 360
        if point > 1:
            Theta = th1 + np.arange(point) * (th2 - th1) / (point - 1)
        else:
            Theta = np.array([th1])

        cV, cH = compute_coefficients(Theta, petitk, grandk, rad_pref2, r, w, P, n_ph, n_th)
        V = cV.mean(axis=0)
        H = cH.mean(axis=0)
        polar_output(fout, "Ellipsoid POLAR PLOT - integrated over scattering angle",
                     V, H)

    else:
        # angle_scattering: Theta grid from shs_angles.py (default 0 … 180 deg)
        try:
            say("Enter the number of points between 0 and Pi for the angular "
                "SHS analysis: ", end="")
            pas2 = rd.integer()
            if pas2 < 1:
                pas2 = 1
        except (EOFError, ValueError):
            return 1

        nang = shs_angles.count(pas2)
        say("********** Begin Ellipsoid SHS Angular Distribution calculation**************")

        fout.write("***************************************************************\n")
        fout.write("                                            \n")
        fout.write("****          ellipsoid SHS angular distribution           ****\n")
        fout.write("                                            \n")
        fout.write(" ****************OUTPUT****************************************\n")
        fout.write(" grandtheta  IPSS-IH(0°)   IPPP-IH(90°)   ISSS-IV(0°)   ISPP-IV(90°) \n")

        Theta = shs_angles.scattering_angles(pas2)
        cV, cH = compute_coefficients(Theta, petitk, grandk, rad_pref2, r, w, P, n_ph, n_th)

        with open(shs_io.path("out_plot"), "w", encoding="utf-8") as fgr:
            fgr.write(" grandtheta   IH(0°)  IH(90°)  IV(0°)  IV(90°) \n")
            for mmm in range(nang):
                deg = Theta[mmm] * 360 / (2 * PI)
                IV_0, IV_90 = cV[mmm, 0], cV[mmm, 2]
                IH_0, IH_90 = cH[mmm, 0], cH[mmm, 2]
                say("angle value: %f " % deg)
                say("IH_0=%.6e    IH_90=%.6e   IV_0=%.6e    IV_90=%.6e "
                    % (IH_0, IH_90, IV_0, IV_90))
                row = " %f    %.6e   %.6e   %.6e   %.6e \n" % (
                    deg, IH_0, IH_90, IV_0, IV_90)
                fout.write(row)
                fgr.write(row)

        say("-----------------------------------------------------------------------------")
        say("calculation finished, to plot the results, open the outplot file with gnuplot")
        say("-----------------------------------------------------------------------------")
        fout.close()

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))

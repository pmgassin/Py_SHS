#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
cube_SHS.py  --  Py_SHS : Coherent SHS by a rectangular box / cube (Python / NumPy)

Python version of cube_SHS.c.  Same physics, same input protocol, same outputs.

Physical model
--------------
Rectangular box a x b x c (cube if a=b=c).  Six flat faces with normals
+-x, +-y, +-z.  Face integrals are done analytically:

  int_face exp(i Phi) dS = Area * exp(i Phi_centre) * sinc(Xu/2) * sinc(Xv/2)

Orientation average over the full rotation group (phi, theta, psi):
  phi, psi : midpoint rule on [0, 2pi)
  cos(theta) : Gauss-Legendre

Keywords : polarplot_single | polarplot_integrate | angle_scattering
Usage    : python cube_SHS.py keyword ki2_input_file output_file
stdin    : NLWKB flag, wavelength, refractive index(es), a, b, c, then angle(s)

Environment
-----------
PYSHS_I_SCALE=x   optional intensity scale factor (default 1).
"""
import os
import sys
import math

import numpy as np

from pyshs import shs_grid
from pyshs import shs_units
from pyshs import shs_angles
from pyshs import shs_io

PI = math.pi

# ---- orientation sampling: minimum / historical sizes (same defaults as the C code) ----
# enlarged automatically for large boxes (see shs_grid.py; PYSHS_GRID=fixed keeps these)
N_PHI = 24
N_THETA = 24
N_PSI = 24
GRID_USED = [N_PHI, N_THETA, N_PSI]   # grid actually used (set in main)
NFACES = 6


# --------------------------------------------------------------------------
# Helpers
# --------------------------------------------------------------------------
class TokenReader:
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
    t1 = np.einsum("ia,abc->ibc", M, beta)
    t2 = np.einsum("jb,ibc->ijc", M, t1)
    return np.einsum("kc,ijc->ijk", M, t2)


def cross3(a, b):
    return np.array([a[1] * b[2] - a[2] * b[1],
                     a[2] * b[0] - a[0] * b[2],
                     a[0] * b[1] - a[1] * b[0]])


def read_ki(name):
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


def gauss_legendre(n):
    """Gauss-Legendre nodes x in (-1,1) and weights w (sum = 2)."""
    x = np.empty(n)
    w = np.empty(n)
    m = (n + 1) // 2
    for i in range(m):
        z = math.cos(PI * (i + 0.75) / (n + 0.5))
        while True:
            p1, p2 = 1.0, 0.0
            for j in range(n):
                p3 = p2
                p2 = p1
                p1 = ((2.0 * j + 1.0) * z * p2 - j * p3) / (j + 1.0)
            pp = n * (z * p1 - p2) / (z * z - 1.0)
            z1 = z
            z = z1 - p1 / pp
            if abs(z - z1) <= 1e-14:
                break
        x[i] = -z
        x[n - 1 - i] = z
        ww = 2.0 / ((1.0 - z * z) * pp * pp)
        w[i] = ww
        w[n - 1 - i] = ww
    return x, w


def csinc(x):
    """sin(x)/x for complex x."""
    if abs(x) < 1e-4:
        x2 = x * x
        return 1.0 - x2 / 6.0 + x2 * x2 / 120.0
    return np.sin(x) / x


def csinc_array(x):
    """sin(x)/x for complex arrays (series near 0)."""
    x = np.asarray(x, dtype=np.complex128)
    small = np.abs(x) < 1e-4
    xs = np.where(small, 1.0, x)
    x2 = x * x
    return np.where(small, 1.0 - x2 / 6.0 + x2 * x2 / 120.0, np.sin(xs) / xs)


# --------------------------------------------------------------------------
# Faces of the box
# --------------------------------------------------------------------------
def build_faces(a, b, c, ki):
    """Build the 6 faces: geometry, local frame, ki(2) in particle frame.

    Returns list of dicts with keys: r0, u, v, area, P.
    """
    ha, hb, hc = a / 2, b / 2, c / 2
    defs = [
        # n, r0, u, v
        (np.array([1., 0, 0]), np.array([ha, 0, 0]), np.array([0., b, 0]), np.array([0., 0, c])),
        (np.array([-1., 0, 0]), np.array([-ha, 0, 0]), np.array([0., b, 0]), np.array([0., 0, c])),
        (np.array([0, 1., 0]), np.array([0, hb, 0]), np.array([a, 0., 0]), np.array([0., 0, c])),
        (np.array([0, -1., 0]), np.array([0, -hb, 0]), np.array([a, 0., 0]), np.array([0., 0, c])),
        (np.array([0, 0, 1.]), np.array([0, 0, hc]), np.array([a, 0., 0]), np.array([0., b, 0])),
        (np.array([0, 0, -1.]), np.array([0, 0, -hc]), np.array([a, 0., 0]), np.array([0., b, 0])),
    ]
    ez = np.array([0., 0., 1.])
    faces = []

    for n, r0, u, v in defs:
        uxv = cross3(u, v)
        area = float(np.linalg.norm(uxv))

        if abs(n[2]) > 0.5:          # caps
            xl = np.array([n[2], 0., 0.])
            yl = np.array([0., 1., 0.])
        else:                        # lateral (cylinder convention)
            xl = np.array([0., 0., -1.])
            yl = cross3(ez, n)

        S = np.column_stack((xl, yl, n))
        # safety: right-handed
        chk = cross3(xl, yl)
        if np.abs(chk - n).sum() > 1e-12:
            raise RuntimeError("local frame is not right-handed")

        P = rotate_tensor(S, ki)
        faces.append({"r0": r0, "u": u, "v": v, "area": area, "P": P})

    return faces


# --------------------------------------------------------------------------
# Core calculation
# --------------------------------------------------------------------------
def compute_coefficients(Theta, petitk, grandk, scale, faces,
                         n_phi=None, n_theta=None, n_psi=None, max_elems=2_000_000):
    """Orientation average for a list of scattering angles (vectorised over orientations).

    Returns coefV[nang,5], coefH[nang,5] already normalised and scaled.
    """
    n_phi = N_PHI if n_phi is None else n_phi
    n_theta = N_THETA if n_theta is None else n_theta
    n_psi = N_PSI if n_psi is None else n_psi
    nang = len(Theta)
    Theta = np.asarray(Theta, dtype=np.float64)
    ugl, wgl = gauss_legendre(n_theta)
    dphi = 2.0 * PI / n_phi
    dpsi = 2.0 * PI / n_psi

    # all orientations
    i, j, k = np.meshgrid(np.arange(n_phi), np.arange(n_theta), np.arange(n_psi),
                          indexing="ij")
    phi = (i.ravel() + 0.5) * dphi
    th = np.arccos(ugl[j.ravel()])
    psi = (k.ravel() + 0.5) * dpsi
    wori_all = wgl[j.ravel()] * dphi * dpsi
    norient = phi.size
    cp, sp, ct, st, cs_, ss_ = (np.cos(phi), np.sin(phi), np.cos(th), np.sin(th),
                                np.cos(psi), np.sin(psi))
    Tall = np.empty((norient, 3, 3))          # same convention as euler_zyz
    Tall[:, 0, 0] = cp * ct * cs_ - sp * ss_
    Tall[:, 0, 1] = -cp * ct * ss_ - sp * cs_
    Tall[:, 0, 2] = cp * st
    Tall[:, 1, 0] = sp * ct * cs_ + cp * ss_
    Tall[:, 1, 1] = -sp * ct * ss_ + cp * cs_
    Tall[:, 1, 2] = sp * st
    Tall[:, 2, 0] = -st * cs_
    Tall[:, 2, 1] = st * ss_
    Tall[:, 2, 2] = ct

    Pf = np.array([f["P"] for f in faces])          # (F,3,3,3)
    r0 = np.array([f["r0"] for f in faces])          # (F,3)
    uu = np.array([f["u"] for f in faces])
    vv = np.array([f["v"] for f in faces])
    area = np.array([f["area"] for f in faces])      # (F,)

    coefV = np.zeros((nang, 5))
    coefH = np.zeros((nang, 5))
    wsum = float(wori_all.sum())

    cs = np.cos(Theta)
    sn = np.sin(Theta)
    qz_arr = 2.0 * petitk - grandk * cs               # (nang,)
    qy_arr = -grandk * sn

    def add_pol(X, weight, coef):        # X (nang, B, 3) ; weight (B,)
        coef[:, 0] += np.einsum("b,ab->a", weight, np.abs(X[:, :, 0]) ** 2)
        coef[:, 1] += np.einsum("b,ab->a", weight,
                                4.0 * np.abs(X[:, :, 2]) ** 2
                                + 2.0 * np.real(X[:, :, 0] * np.conj(X[:, :, 1])))
        coef[:, 2] += np.einsum("b,ab->a", weight, np.abs(X[:, :, 1]) ** 2)
        coef[:, 3] += np.einsum("b,ab->a", weight, 4.0 * np.real(X[:, :, 0] * np.conj(X[:, :, 2])))
        coef[:, 4] += np.einsum("b,ab->a", weight, 4.0 * np.real(X[:, :, 1] * np.conj(X[:, :, 2])))

    nf = len(faces)
    batch = max(1, min(norient, max_elems // max(1, nang * nf)))
    done = 0
    for b0 in range(0, norient, batch):
        b1 = min(b0 + batch, norient)
        say("calculation in progress : %5.1f pourcent " % (100.0 * b0 / norient))
        T = Tall[b0:b1]
        w = wori_all[b0:b1]
        X, Y = T[:, 0, :], T[:, 1, :]
        # lab-frame tensors, (xx, yy, xy), per face: L[b,f,i,q]
        Mxx = np.einsum("fabc,pb,pc->pfa", Pf, X, X)
        Myy = np.einsum("fabc,pb,pc->pfa", Pf, Y, Y)
        Mxy = 0.5 * (np.einsum("fabc,pb,pc->pfa", Pf, X, Y)
                     + np.einsum("fabc,pb,pc->pfa", Pf, Y, X))
        L = np.stack([np.einsum("pia,pfa->pfi", T, M) for M in (Mxx, Myy, Mxy)], axis=-1)
        # lab-frame y,z of centre and edge vectors: (B, F)
        yc = np.einsum("pa,fa->pf", T[:, 1, :], r0)
        zc = np.einsum("pa,fa->pf", T[:, 2, :], r0)
        yu = np.einsum("pa,fa->pf", T[:, 1, :], uu)
        zu = np.einsum("pa,fa->pf", T[:, 2, :], uu)
        yv = np.einsum("pa,fa->pf", T[:, 1, :], vv)
        zv = np.einsum("pa,fa->pf", T[:, 2, :], vv)

        qz = qz_arr[:, None, None]
        qy = qy_arr[:, None, None]
        Xu = qy * yu[None] + qz * zu[None]
        Xv = qy * yv[None] + qz * zv[None]
        Ff = (area[None, None, :] * np.exp(1j * (qz * zc[None] + qy * yc[None]))
              * csinc_array(0.5 * Xu) * csinc_array(0.5 * Xv))          # (nang, B, F)
        A = np.einsum("apf,pfiq->apiq", Ff, L)                          # (nang, B, 3, 3)

        add_pol(A[:, :, 0, :], w, coefV)
        Hs = cs[:, None, None] * A[:, :, 1, :] - sn[:, None, None] * A[:, :, 2, :]
        add_pol(Hs, w, coefH)

    coefV *= scale / wsum
    coefH *= scale / wsum
    return coefV, coefH


# --------------------------------------------------------------------------
# Outputs
# --------------------------------------------------------------------------
def write_header(fout, use_nlwkb, lam, n_p, n_m, a, b, c, petitk, grandk,
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
    fout.write("**       SIMULATION RESULTS  for box (cube) programm ****\n")
    fout.write("********************************************************\n")
    fout.write("  Mode                : %s\n" % ("NLWKB" if use_nlwkb else "Classic RGD"))
    fout.write("  wavelength = %f nm\n" % lam)
    if use_nlwkb:
        fout.write("  n_p = %.5f + i*%.5f   n_m = %.5f + i*%.5f\n"
                   % (n_p.real, n_p.imag, n_m.real, n_m.imag))
    else:
        fout.write("  n = %.5f + i*%.5f\n" % (n_p.real, n_p.imag))
    fout.write("  box edges  a = %f nm   b = %f nm   c = %f nm\n" % (a, b, c))
    fout.write("  petitk  (k_omega)   = %.6f + i*%.6f\n" % (petitk.real, petitk.imag))
    fout.write("  grandk  (k_2omega)  = %.6f + i*%.6f\n" % (grandk.real, grandk.imag))
    fout.write("  Radiative pref.     : |k|^2/(4pi) = %.6e  (I *= pref^2 = %.6e)\n"
               % (rad_pref, rad_pref2))
    fout.write("  Quadrature          : exact face integrals, orientations "
               "%d (phi) x %d (cos theta, Gauss-Legendre) x %d (psi)\n"
               % tuple(GRID_USED))
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
    say("      Box (cube) SHS module with optional NLWKB         ")
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

        say("Enter the edge a of the box along x in nm: ", end="")
        a = rd.real()
        say("Enter the edge b of the box along y in nm: ", end="")
        b = rd.real()
        say("Enter the edge c of the box along z (axis) in nm (a=b=c : cube): ",
            end="")
        c = rd.real()
        say("      ******************************     ")
    except (EOFError, ValueError):
        return 1

    if a <= 0 or b <= 0 or c <= 0 or lam <= 0:
        sys.stderr.write("wavelength and the three edges must be strictly positive.\n")
        return 1

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

    n_ph, n_th, n_ps, grid_msg = shs_grid.box_grid(
        a, b, c, petitk, grandk, (N_PHI, N_THETA, N_PSI))
    GRID_USED[:] = [n_ph, n_th, n_ps]
    say("Orientation grid: %d x %d x %d" % (n_ph, n_th, n_ps))
    if grid_msg:
        say(grid_msg)

    write_header(fout, use_nlwkb, lam, n_p, n_m, a, b, c, petitk, grandk,
                 rad_pref, rad_pref2, ki)
    shs_units.write_header_lines(fout, unit)

    faces = build_faces(a, b, c, ki)

    # ---- modes ----
    if mode_single:
        try:
            say("Enter the scattered angle in degree (0 = transmission): ", end="")
            gdeg = rd.real()
        except (EOFError, ValueError):
            return 1
        Theta = np.array([gdeg * 2 * PI / 360])
        say("********** Begin calculation for box polar plot at %.1f ° ********"
            % gdeg)
        cV, cH = compute_coefficients(Theta, petitk, grandk, rad_pref2, faces, *GRID_USED)
        title = "Box POLAR PLOT - CONFIGURATION AT %.1f deg" % gdeg
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

        cV, cH = compute_coefficients(Theta, petitk, grandk, rad_pref2, faces, *GRID_USED)
        V = cV.mean(axis=0)
        H = cH.mean(axis=0)
        polar_output(fout, "Box POLAR PLOT - integrated over scattering angle", V, H)

    else:
        try:
            say("Enter the number of points between 0 and Pi for the angular "
                "SHS analysis: ", end="")
            pas2 = rd.integer()
            if pas2 < 1:
                pas2 = 1
        except (EOFError, ValueError):
            return 1

        nang = shs_angles.count(pas2)
        say("********** Begin Box SHS Angular Distribution calculation**************")

        fout.write("***************************************************************\n")
        fout.write("                                            \n")
        fout.write("****             box SHS angular distribution              ****\n")
        fout.write("                                            \n")
        fout.write(" ****************OUTPUT****************************************\n")
        fout.write(" grandtheta  IPSS-IH(0°)   IPPP-IH(90°)   ISSS-IV(0°)   ISPP-IV(90°) \n")

        Theta = shs_angles.scattering_angles(pas2)
        cV, cH = compute_coefficients(Theta, petitk, grandk, rad_pref2, faces, *GRID_USED)

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

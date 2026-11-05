#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
sphere_SHS.py  --  Py_SHS : Coherent SHS by a colloidal sphere (Python / NumPy)

Python version of sphere_SHS.c.  Same physics, same input protocol,
same output files and same numerical results.

Physical model (manual section 2.4)
-----------------------------------
Continuum surface response of a sphere of radius R under the NLWKB
(or classical Rayleigh-Gans-Debye) approximation.

  Gamma_IJK(Theta) = ∫∫_surface  chi_IJK(phi_c, theta_c)
                     * exp(i Phi)  R^2 sin(theta_c) d phi_c d theta_c

with Phi = (2 k_omega - K_2omega cos Theta) R cos theta_c
         + K_2omega sin Theta * R sin theta_c * sin phi_c
  (sign convention of the C code: e^{+i(...)} with the Y-term negative
   in the C source, reproduced here).

chi in the particle frame is obtained from the surface-frame chi^(2)
(beta input) by a ZYZ rotation with psi = 0 at each surface point.

Because a sphere is isotropic, no orientational average over the particle
is required.

Keywords : polarplot_single | polarplot_integrate | angle_scattering
Usage    : python sphere_SHS.py keyword beta_file output_file
stdin    : NLWKB flag, wavelength, refractive index(es), radius, then angle(s)

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

# Surface quadrature: minimum / historical grid (C code: pas = 40 -> nx=80, ny=40).
# The grid is enlarged automatically for large spheres (see shs_grid.py).
PAS_SURFACE = 40


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


def euler_matrices(phi, theta, psi):
    """ZYZ Euler matrices (identical to the C Txx..Tzz).
    phi, theta, psi : arrays  ->  T : (..., 3, 3)"""
    sp, cp = np.sin(phi), np.cos(phi)
    st, ct = np.sin(theta), np.cos(theta)
    ss, cs = np.sin(psi), np.cos(psi)
    T = np.empty(np.shape(phi) + (3, 3))
    T[..., 0, 0] = -sp * ss + ct * cp * cs
    T[..., 0, 1] = -sp * cs - ct * cp * ss
    T[..., 0, 2] = st * cp
    T[..., 1, 0] = cp * ss + ct * sp * cs
    T[..., 1, 1] = cp * cs - ct * ss * sp
    T[..., 1, 2] = st * sp
    T[..., 2, 0] = -st * cs
    T[..., 2, 1] = st * ss
    T[..., 2, 2] = ct
    return T


def rotate_tensor_pointwise(T, beta):
    """out[..., i, j, k] = sum_abc T[..., i, a] T[..., j, b] T[..., k, c] beta[a,b,c]
    Vectorised over the leading dimensions of T."""
    t1 = np.einsum("...ia,abc->...ibc", T, beta)
    t2 = np.einsum("...jb,...ibc->...ijc", T, t1)
    return np.einsum("...kc,...ijc->...ijk", T, t2)


def read_beta(name):
    """1 header line + 27 lines 'name value'."""
    beta = np.zeros((3, 3, 3))
    try:
        with open(name, "rt") as f:
            f.readline()
            for i in range(3):
                for j in range(3):
                    for k in range(3):
                        beta[i, j, k] = float(f.readline().split()[1])
    except (OSError, IndexError, ValueError):
        return None
    # The intensity formulas use beta_{i,jk} only for (j,k) = (x,x), (y,y), (x,y) and
    # assume beta_{ijk} = beta_{ikj} (intrinsic permutation symmetry).  Symmetrise in the
    # last two indices so that every module treats any input tensor the same way.
    asym = np.abs(beta - beta.transpose(0, 2, 1)).max()
    if asym > 1e-12 * max(np.abs(beta).max(), 1e-300):
        say("WARNING: the tensor is not symmetric in its last two indices "
            "(max |beta_ijk - beta_ikj| = %.3e): it has been symmetrised, "
            "beta_ijk = beta_ikj = (beta_ijk + beta_ikj)/2." % asym)
        beta = 0.5 * (beta + beta.transpose(0, 2, 1))
    return beta


# --------------------------------------------------------------------------
# Surface integral → effective Gamma components
# --------------------------------------------------------------------------
# The nine components needed for the polarisation coefficients
# (I, j, k) in particle / lab frame after surface rotation with psi=0:
_COMPONENTS = {
    "xxx": (0, 0, 0),
    "yxx": (1, 0, 0),
    "zxx": (2, 0, 0),
    "xyy": (0, 1, 1),
    "yyy": (1, 1, 1),
    "zyy": (2, 1, 1),
    "xxy": (0, 0, 1),
    "yxy": (1, 0, 1),
    "zxy": (2, 0, 1),
}


def surface_gamma(R, grandtheta, beta, petitk, grandk, pas=PAS_SURFACE):
    """Mid-point quadrature of Gamma_IJK over the sphere surface.

    Returns a dict of complex scalars keyed by component name
    (xxx, yxx, zxx, xyy, yyy, zyy, xxy, yxy, zxy).
    """
    nx = 2 * pas
    ny = pas
    hx = (2.0 * PI) / nx

    # phi_c: periodic mid-point (trapezoid) rule; theta_c: Gauss-Legendre in cos(theta_c)
    # (PYSHS_THETA_RULE=midpoint restores the former mid-point rule in theta_c)
    phi_c = hx / 2.0 + np.arange(nx) * hx          # (nx,)
    if os.environ.get("PYSHS_THETA_RULE", "gauss").strip().lower() == "midpoint":
        hy = PI / ny
        theta_c = hy / 2.0 + np.arange(ny) * hy    # (ny,)
        w_theta = hy * np.sin(theta_c)
    else:
        xg, w_theta = np.polynomial.legendre.leggauss(ny)
        theta_c = np.arccos(xg)
    Phi, Theta = np.meshgrid(phi_c, theta_c, indexing="ij")  # (nx, ny)

    # Rotation matrices at each surface point (psi = 0)
    T = euler_matrices(Phi, Theta, np.zeros_like(Phi))      # (nx, ny, 3, 3)

    # Surface chi in the particle frame
    chi = rotate_tensor_pointwise(T, beta)                   # (nx, ny, 3, 3, 3)

    # Optical phase (identical to the C code convention)
    # Phi_opt = (2*petitk - grandk*cos(Θ)) * R * cos(θ_c)
    #         - (grandk*sin(Θ)) * R * sin(θ_c) * sin(φ_c)
    phase = ((2.0 * petitk - grandk * np.cos(grandtheta)) * R * np.cos(Theta)
             - (grandk * np.sin(grandtheta)) * R * np.sin(Theta) * np.sin(Phi))
    weight = (R * R * w_theta[None, :] * np.exp(1j * phase)
              * hx)                                          # (nx, ny)

    out = {}
    for name, (I, J, K) in _COMPONENTS.items():
        out[name] = np.sum(chi[:, :, I, J, K] * weight)
    return out


def coefficients_from_gamma(g, grandtheta):
    """Compute (av, bv, cv, ah, bh, ch) from the nine Gamma components."""
    c = np.cos(grandtheta)
    s = np.sin(grandtheta)

    gxxx, gxyy, gxxy = g["xxx"], g["xyy"], g["xxy"]
    gyxx, gyyy, gyxy = g["yxx"], g["yyy"], g["yxy"]
    gzxx, gzyy, gzxy = g["zxx"], g["zyy"], g["zxy"]

    av = np.real(gxxx * np.conj(gxxx))
    cv = np.real(gxyy * np.conj(gxyy))
    bv = np.real(4.0 * gxxy * np.conj(gxxy)
                 + gxxx * np.conj(gxyy) + gxyy * np.conj(gxxx))

    ah = np.real(gyxx * np.conj(gyxx) * c * c
                 + s * s * gzxx * np.conj(gzxx)
                 - c * s * (gyxx * np.conj(gzxx) + gzxx * np.conj(gyxx)))

    ch = np.real(gyyy * np.conj(gyyy) * c * c
                 + s * s * gzyy * np.conj(gzyy)
                 - c * s * (gyyy * np.conj(gzyy) + gzyy * np.conj(gyyy)))

    bh = np.real(
        c * c * (4.0 * gyxy * np.conj(gyxy)
                 + gyxx * np.conj(gyyy) + gyyy * np.conj(gyxx))
        + s * s * (4.0 * gzxy * np.conj(gzxy)
                   + gzxx * np.conj(gzyy) + gzyy * np.conj(gzxx))
        - c * s * (4.0 * gyxy * np.conj(gzxy) + 4.0 * gzxy * np.conj(gyxy)
                   + gyxx * np.conj(gzyy) + gzyy * np.conj(gyxx)
                   + gzxx * np.conj(gyyy) + gyyy * np.conj(gzxx))
    )
    return av, bv, cv, ah, bh, ch


def intensity_VH_extremes(g, grandtheta):
    """IV(0), IV(90), IH(0), IH(90) from Gamma (for angle_scattering)."""
    c = np.cos(grandtheta)
    s = np.sin(grandtheta)

    gxxx, gxyy = g["xxx"], g["xyy"]
    gyxx, gyyy = g["yxx"], g["yyy"]
    gzxx, gzyy = g["zxx"], g["zyy"]

    IV_0 = np.real(gxxx * np.conj(gxxx))
    IV_90 = np.real(gxyy * np.conj(gxyy))
    IH_0 = np.real(gyxx * np.conj(gyxx) * c * c
                   + s * s * gzxx * np.conj(gzxx)
                   - c * s * (gyxx * np.conj(gzxx) + gzxx * np.conj(gyxx)))
    IH_90 = np.real(gyyy * np.conj(gyyy) * c * c
                    + s * s * gzyy * np.conj(gzyy)
                    - c * s * (gyyy * np.conj(gzyy) + gzyy * np.conj(gyyy)))
    return IV_0, IV_90, IH_0, IH_90


# --------------------------------------------------------------------------
# Outputs
# --------------------------------------------------------------------------
def write_header(fout, beta, use_nlwkb, lam, n_p, n_m, R, petitk, grandk,
                 rad_pref, rad_pref2, pas):
    fout.write("********************************************************\n")
    fout.write("****           Py_SHS: Sphere SHS module            ****\n")
    fout.write("****         with optional NLWKB approximation      ****\n")
    fout.write("********************************************************\n")
    fout.write("Mode                : %s\n" % ("NLWKB" if use_nlwkb else "Classic RGD"))
    fout.write("Wavelength          : %.2f nm\n" % lam)
    if use_nlwkb:
        fout.write("n_particle (n_p)    : %.5f + i*%.5f\n" % (n_p.real, n_p.imag))
        fout.write("n_medium   (n_m)    : %.5f + i*%.5f\n" % (n_m.real, n_m.imag))
    else:
        fout.write("Refractive index n  : %.5f + i*%.5f\n" % (n_p.real, n_p.imag))
    fout.write("Radius              : %.2f nm\n" % R)
    fout.write("petitk              : %.6f + i*%.6f\n" % (petitk.real, petitk.imag))
    fout.write("grandk              : %.6f + i*%.6f\n" % (grandk.real, grandk.imag))
    fout.write("Radiative pref.     : |k|^2/(4pi) = %.6e  (I *= pref^2 = %.6e)\n"
               % (rad_pref, rad_pref2))
    fout.write("Surface quadrature  : pas = %d (nx=%d, ny=%d)\n"
               % (pas, 2 * pas, pas))
    fout.write("********************************************************\n\n")
    fout.write("*********          SURFACE KI(2) TENSOR      ***********\n\n")
    b = beta
    fout.write("    xxx=%f   xxy=%f   xxz= %f\n" % (b[0, 0, 0], b[0, 0, 1], b[0, 0, 2]))
    fout.write("    xyx=%f   xyy=%f   xyz= %f\n" % (b[0, 1, 0], b[0, 1, 1], b[0, 1, 2]))
    fout.write("    xzx=%f   xzy=%f   xzz= %f\n" % (b[0, 2, 0], b[0, 2, 1], b[0, 2, 2]))
    fout.write("\n")
    fout.write("    yxx=%f   yxy=%f   yxz= %f\n" % (b[1, 0, 0], b[1, 0, 1], b[1, 0, 2]))
    fout.write("    yyx=%f   yyy=%f   yyz= %f\n" % (b[1, 1, 0], b[1, 1, 1], b[1, 1, 2]))
    fout.write("    yzx=%f   yzy=%f   yzz= %f\n" % (b[1, 2, 0], b[1, 2, 1], b[1, 2, 2]))
    fout.write("\n")
    fout.write("    zxx=%f   zxy=%f   zxz= %f\n" % (b[2, 0, 0], b[2, 0, 1], b[2, 0, 2]))
    fout.write("    zyx=%f   zyy=%f   zyz= %f\n" % (b[2, 1, 0], b[2, 1, 1], b[2, 1, 2]))
    fout.write("    zzx=%f   zzy=%f   zzz= %f\n" % (b[2, 2, 0], b[2, 2, 1], b[2, 2, 2]))
    fout.write("\n")
    fout.write("*******************************************************\n")


def report_polar(fout, deg, av, bv, cv, ah, bh, ch, label="Sphere POLAR PLOT"):
    with np.errstate(all="ignore"):
        I4v = (av - bv + cv) / (3 * av + bv + 3 * cv)
        I2v = 4.0 * (av - cv) / (3 * av + bv + 3 * cv)
        I4h = (ah - bh + ch) / (3 * ah + bh + 3 * ch)
        I2h = 4.0 * (ah - ch) / (3 * ah + bh + 3 * ch)

    say("")
    say("********** calculation complete ************")
    say("av = %.6e" % av)
    say("bv = %.6e" % bv)
    say("cv = %.6e" % cv)
    say("ah = %.6e" % ah)
    say("bh = %.6e" % bh)
    say("ch = %.6e" % ch)
    say("I2v = %f" % I2v)
    say("I4v = %f" % I4v)
    say("I2h = %f" % I2h)
    say("I4h = %f" % I4h)

    fout.write("*********************************************\n")
    fout.write("* %s - CONFIGURATION AT %.1f °\n" % (label, deg))
    fout.write("*********************************************\n")
    fout.write("I2v=%f  I4v=%f\nI2h=%f  I4h=%f\n" % (I2v, I4v, I2h, I4h))
    fout.write("av=%f   bv=%f   cv=%f\nah=%f   bh=%f   ch=%f\n"
               % (av, bv, cv, ah, bh, ch))

    with open(shs_io.path("out_plot"), "w", encoding="utf-8") as fp:
        fp.write("gamma IV(gamma) IH(gamma)\n")
        for iii in range(100):
            g = 2.0 * PI * iii / 100.0
            nc, ns = np.cos(g), np.sin(g)
            Iv = av * nc**4 + bv * nc**2 * ns**2 + cv * ns**4
            Ih = ah * nc**4 + bh * nc**2 * ns**2 + ch * ns**4
            fp.write("%.8e %.8e %.8e\n" % (g, Iv, Ih))
    say("-----------------------------------------------------------------------------")
    say("calculation finished, open the out_plot file with gnuplot")
    say("-----------------------------------------------------------------------------")


# --------------------------------------------------------------------------
def main(argv):
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except AttributeError:
        pass

    say("********************************************************")
    say("****           Py_SHS: Sphere SHS module            ****")
    say("****         with optional NLWKB approximation      ****")
    say("**** Dr Pierre-Marie GASSIN and Dr Gaelle GASSIN    ****")
    say("****                  ICGM - ENSCM                  ****")
    say("********************************************************")
    say("")

    if len(argv) < 4:
        say("Usage: %s <keyword> <ki2_inputfile> <outputfile>" % argv[0])
        say("keyword = polarplot_single | polarplot_integrate | angle_scattering")
        return 1

    keyword = argv[1]
    mode_single = keyword == "polarplot_single"
    mode_integ = keyword == "polarplot_integrate"
    mode_angle = keyword == "angle_scattering"
    if not (mode_single or mode_integ or mode_angle):
        say("Unknown keyword. Use: polarplot_single | polarplot_integrate | angle_scattering")
        return 1

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
            nr_p = rd.real()
            say("Enter IMAGINARY part of particle refractive index n_p: ", end="")
            ni_p = rd.real()
            n_p = complex(nr_p, ni_p)

            say("Enter REAL part of medium refractive index n_m: ", end="")
            nr_m = rd.real()
            say("Enter IMAGINARY part of medium refractive index n_m: ", end="")
            ni_m = rd.real()
            n_m = complex(nr_m, ni_m)

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
            n = complex(nr, ni)
            n_p = n
            n_m = n
            say("")
            say("n = %.5f + i*%.5f" % (n.real, n.imag))

        say("Enter the radius of the sphere in nm: ", end="")
        R = rd.real()
        say("")
    except (EOFError, ValueError):
        return 1

    # Wave vectors
    if use_nlwkb:
        petitk = 2.0 * PI * n_p / lam          # k_ω inside particle
        grandk = 4.0 * PI * n_m / lam          # K_2ω in medium
        say("NLWKB vectors:")
    else:
        petitk = 2.0 * PI * n_p / lam
        grandk = 4.0 * PI * n_m / lam
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

    beta = read_beta(argv[2])
    if beta is None:
        say("Problem to open the input ki(2) file!")
        return 1

    shs_io.set_output(argv[3])

    try:
        fout = open(argv[3], "w", encoding="utf-8")
    except OSError:
        say("Problem to open the output file!")
        return 1

    pas, grid_msg = shs_grid.sphere_pas(R, petitk, grandk, PAS_SURFACE)
    say("Surface grid: nx = %d, ny = %d" % (2 * pas, pas))
    if grid_msg:
        say(grid_msg)
    write_header(fout, beta, use_nlwkb, lam, n_p, n_m, R, petitk, grandk,
                 rad_pref, rad_pref2, pas)
    shs_units.write_header_lines(fout, unit)

    # ------------------------------------------------------------------
    # MODE 1 : polarplot_single
    # ------------------------------------------------------------------
    if mode_single:
        try:
            say("Enter the scattered angle in degree (0 = transmission): ", end="")
            grandthetadegree = rd.real()
        except (EOFError, ValueError):
            return 1
        grandtheta = grandthetadegree * 2.0 * PI / 360.0
        say("********** Begin calculation for sphere polar plot at %.1f ° ********"
            % grandthetadegree)

        g = surface_gamma(R, grandtheta, beta, petitk, grandk, pas)
        av, bv, cv, ah, bh, ch = coefficients_from_gamma(g, grandtheta)
        av *= rad_pref2
        bv *= rad_pref2
        cv *= rad_pref2
        ah *= rad_pref2
        bh *= rad_pref2
        ch *= rad_pref2
        report_polar(fout, grandthetadegree, av, bv, cv, ah, bh, ch)

    # ------------------------------------------------------------------
    # MODE 2 : polarplot_integrate
    # ------------------------------------------------------------------
    elif mode_integ:
        try:
            say("Enter the beginning scattered angle in degree (0 = transmission): ",
                end="")
            grandthetadegree_1 = rd.real()
            say("Enter the final scattered angle in degree (0 = transmission): ",
                end="")
            grandthetadegree_2 = rd.real()
            say("Number of integrale points: ", end="")
            point = rd.integer()
            if point < 1:
                return 1
        except (EOFError, ValueError):
            return 1

        grandtheta_1 = grandthetadegree_1 * 2.0 * PI / 360.0
        grandtheta_2 = grandthetadegree_2 * 2.0 * PI / 360.0

        av = bv = cv = ah = bh = ch = 0.0
        for ic in range(point):
            if point == 1:
                grandtheta = grandtheta_1
            else:
                grandtheta = grandtheta_1 + ic * (grandtheta_2 - grandtheta_1) / (point - 1.0)
            grandthetadegree = grandtheta * 360.0 / (2.0 * PI)
            say("********** Begin calculation for sphere polar plot at %.1f ° ********"
                % grandthetadegree)

            g = surface_gamma(R, grandtheta, beta, petitk, grandk, pas)
            a, b, c, aa, bb, cc = coefficients_from_gamma(g, grandtheta)
            av += a
            bv += b
            cv += c
            ah += aa
            bh += bb
            ch += cc

        av /= point
        bv /= point
        cv /= point
        ah /= point
        bh /= point
        ch /= point
        av *= rad_pref2
        bv *= rad_pref2
        cv *= rad_pref2
        ah *= rad_pref2
        bh *= rad_pref2
        ch *= rad_pref2

        fout.write("***********************************************************************\n")
        fout.write("* Sphere POLAR PLOT - Average between %.1f ° and %.1f °\n"
                   % (grandthetadegree_1, grandthetadegree_2))
        report_polar(fout, grandthetadegree_2, av, bv, cv, ah, bh, ch,
                     label="Sphere POLAR PLOT (integrated)")

    # ------------------------------------------------------------------
    # MODE 3 : angle_scattering
    # ------------------------------------------------------------------
    else:
        try:
            say("Enter the number of points between 0 and Pi for the angular SHS analysis: ",
                end="")
            pas2 = rd.integer()
            if pas2 < 1:
                return 1
        except (EOFError, ValueError):
            return 1

        th_grid = shs_angles.scattering_angles(pas2)
        mmmf = len(th_grid)
        say("********** Begin Sphere SHS Angular Distribution calculation **************")

        fout.write("***************************************************************\n")
        fout.write("*****            Sphere SHS angular distribution          *****\n")
        fout.write(" ****************OUTPUT****************************************\n")
        fout.write(" grandtheta  IPSS-IH(0°)   IPPP-IH(90°)   ISSS-IV(0°)   ISPP-IV(90°) \n")

        with open(shs_io.path("out_plot"), "w", encoding="utf-8") as fgr:
            fgr.write(" grandtheta   IH(0°)  IH(90°)  IV(0°)  IV(90°) \n")
            for mmm in range(mmmf):
                grandtheta = th_grid[mmm]
                grandthetadegree = grandtheta * 360.0 / (2.0 * PI)
                say("angle value: %.1f °" % grandthetadegree)

                g = surface_gamma(R, grandtheta, beta, petitk, grandk, pas)
                IV_0, IV_90, IH_0, IH_90 = intensity_VH_extremes(g, grandtheta)
                IH_0 *= rad_pref2
                IH_90 *= rad_pref2
                IV_0 *= rad_pref2
                IV_90 *= rad_pref2

                say("IH_0=%.6e    IH_90=%.6e   IV_0=%.6e    IV_90=%.6e"
                    % (IH_0, IH_90, IV_0, IV_90))
                row = " %f    %.6e   %.6e   %.6e   %.6e\n" % (
                    grandthetadegree, IH_0, IH_90, IV_0, IV_90)
                fout.write(row)
                fgr.write(row)

        say("-----------------------------------------------------------------------------")
        say("calculation finished, open the out_plot file with gnuplot")
        say("-----------------------------------------------------------------------------")

    fout.close()
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))

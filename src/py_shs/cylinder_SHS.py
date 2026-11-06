#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
cylinder_SHS.py  --  Py_SHS : Coherent SHS by a colloidal cylinder (Python / NumPy)

Python version of cylinder_SHS_rapide.c.  Same physics, same input protocol,
same output files and same numerical organisation.

Physical model (manual section 2.5)
-----------------------------------
Continuum surface response of a finite cylinder (radius R, length L) under
the NLWKB approximation.  The surface chi^(2) is integrated over the lateral
surface with the optical phase, then averaged over the orientation of the
cylinder axis (phi, theta) with weight sin(theta).  Only the LATERAL surface is
modelled: the two end caps are not included, so the model is meant for objects
long compared with their radius (for short cylinders, L ~ R, the caps matter).

Optimisation (same as the C rewrite)
------------------------------------
  for each orientation (phi, theta)        <- Euler matrix T once
      rotate surface tensor per azimuth    <- once per (orient, azimuth)
      for each scattering angle
          sum over surface points -> A     <- 1 exp per point
          accumulate intensities

Because the surface tensor depends only on the azimuth (not on z), the
rotation to the laboratory frame is done once per (orientation, azimuth).

Keywords : polarplot_single | polarplot_integrate | angle_scattering
Usage    : python cylinder_SHS.py keyword ki_input_file output_file
stdin    : wavelength, n_p, n_m, R, length, then the angle(s) of the mode

Environment
-----------
PYSHS_I_SCALE=x   optional intensity scale factor (default 1).
"""
import os
import sys
import math

import numpy as np

from py_shs import shs_grid
from py_shs import shs_units
from py_shs import shs_angles
from py_shs import shs_io

# ---- discretisation (identical to the C code) ----
# Minimum / historical grid sizes. The grids are enlarged automatically for large objects
# (see shs_grid.py; PYSHS_GRID=fixed keeps these values).  The axial integral is analytic,
# so the former N_Z (surface points along the axis) no longer exists.
N_BETA = 10          # surface points around the axis (minimum)
N_PHI = 10           # orientation grid, phi in [0, 2pi) (minimum)
N_THETA = 10         # orientation grid, theta in [0, pi] (minimum)
PI = math.pi

# Polar angle of the surface frame: the local z axis (surface normal) is radial, i.e.
# perpendicular to the cylinder axis.  The former C code used the truncated value 1.570796
# (2.7e-7 rad away from pi/2, a relative effect of ~1e-7 on the intensities); set
# PYSHS_THETA_MESO_LEGACY=1 to get it back.
THETA_MESO = 1.570796 if os.environ.get("PYSHS_THETA_MESO_LEGACY") == "1" else math.pi / 2


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


def euler_matrix(phi, theta, psi):
    """Single ZYZ Euler matrix (3, 3)."""
    sp, cp = math.sin(phi), math.cos(phi)
    st, ct = math.sin(theta), math.cos(theta)
    ss, cs = math.sin(psi), math.cos(psi)
    T = np.empty((3, 3))
    T[0, 0] = -sp * ss + ct * cp * cs
    T[0, 1] = -sp * cs - ct * cp * ss
    T[0, 2] = st * cp
    T[1, 0] = cp * ss + ct * sp * cs
    T[1, 1] = cp * cs - ct * ss * sp
    T[1, 2] = st * sp
    T[2, 0] = -st * cs
    T[2, 1] = st * ss
    T[2, 2] = ct
    return T


def rotate_tensor(T, beta):
    """out[i,j,k] = sum_abc T[i,a] T[j,b] T[k,c] beta[a,b,c]."""
    t1 = np.einsum("ia,abc->ibc", T, beta)
    t2 = np.einsum("jb,ibc->ijc", T, t1)
    return np.einsum("kc,ijc->ijk", T, t2)


def read_ki(name):
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



def _csinc(x):
    """sin(x)/x for complex arrays (series near 0)."""
    x = np.asarray(x, dtype=np.complex128)
    small = np.abs(x) < 1e-4
    xs = np.where(small, 1.0, x)
    out = np.sin(xs) / xs
    x2 = x * x
    return np.where(small, 1.0 - x2 / 6.0 + x2 * x2 / 120.0, out)


def compute(R, length, kisurf, petitk, grandk, theta_s,
            n_beta=None, n_phi=None, n_theta=None):
    """Full cylinder SHS calculation.

    The lateral surface is integrated analytically along the axis (sinc factor, exact for
    any length) and with a periodic trapezoid rule around it (n_beta points).  The
    orientation average uses n_phi trapezoid points in phi and n_theta Gauss-Legendre
    points in cos(theta).

    Returns (acc, wsum) where
      acc[nang, 2, 5]  = coefficients (a,b,c,b1,b2) for V and H
      wsum             = sum of orientation weights
    """
    n_beta = N_BETA if n_beta is None else n_beta
    n_phi = N_PHI if n_phi is None else n_phi
    n_theta = N_THETA if n_theta is None else n_theta
    nang = len(theta_s)
    theta_s = np.asarray(theta_s, dtype=np.float64)

    hb = 2.0 * PI / n_beta
    hphi = 2.0 * PI / n_phi
    th_nodes, th_w = theta_rule(n_theta)

    # Angle-dependent constants
    cs = np.cos(theta_s)
    sn = np.sin(theta_s)
    qz = 2.0 * petitk - grandk * cs          # complex, (nang,)
    qy = -grandk * sn

    # Surface tensor in the cylinder frame, one per azimuth
    # (old code: kimesoijk(thetac, 1.570796, 0, ...) for every point)
    Pb = np.empty((n_beta, 3, 3, 3))
    xb = np.empty(n_beta)
    yb = np.empty(n_beta)
    for ib in range(n_beta):
        b = 0.5 * hb + ib * hb
        Tm = euler_matrix(b, THETA_MESO, 0.0)
        Pb[ib] = rotate_tensor(Tm, kisurf)
        xb[ib] = R * math.cos(b)
        yb[ib] = R * math.sin(b)

    wsurf = R * hb * length            # weight of one azimuth strip (axial integral
                                       # supplies the sinc factor)
    acc = np.zeros((nang, 2, 5))   # [angle, {V,H}, {a,b,c,b1,b2}]
    wsum = 0.0

    for ip in range(n_phi):
        say("calculation in progress : %f pourcent " % (ip * 100.0 / n_phi))
        phi = 0.5 * hphi + ip * hphi
        for it in range(n_theta):
            th = th_nodes[it]
            w = hphi * th_w[it]
            wsum += w

            T = euler_matrix(phi, th, 0.0)

            # lab-frame tensors, all azimuths at once: full[n,i,j,k]
            full = np.einsum("ia,jb,kc,nabc->nijk", T, T, T, Pb)
            Lb = np.empty((n_beta, 3, 3))       # [ib, i, q]  q=xx,yy,xy
            Lb[:, :, 0] = wsurf * full[:, :, 0, 0]
            Lb[:, :, 1] = wsurf * full[:, :, 1, 1]
            Lb[:, :, 2] = wsurf * full[:, :, 0, 1]
            byl = T[1, 0] * xb + T[1, 1] * yb   # (n_beta,)
            bzl = T[2, 0] * xb + T[2, 1] * yb
            ty = T[1, 2]
            tz = T[2, 2]

            # phase of each azimuth strip centre, (nang, n_beta)
            ph = np.exp(1j * (qz[:, None] * bzl[None, :] + qy[:, None] * byl[None, :]))
            # axial integral: int_{-L/2}^{L/2} exp(i kappa z) dz / L = sinc(kappa L / 2)
            kappa = qz * tz + qy * ty
            fz = _csinc(0.5 * kappa * length)                       # (nang,)
            A = np.einsum("ab,biq->aiq", ph, Lb) * fz[:, None, None]  # (nang,3,3)

            def pol_sums(xx, yy, xy):
                return np.stack([
                    w * np.abs(xx) ** 2,
                    w * (4.0 * np.abs(xy) ** 2 + 2.0 * np.real(xx * np.conj(yy))),
                    w * np.abs(yy) ** 2,
                    w * 4.0 * np.real(xx * np.conj(xy)),
                    w * 4.0 * np.real(yy * np.conj(xy)),
                ], axis=-1)                                          # (nang, 5)

            # V : lab x ; H : cos(Theta)*y - sin(Theta)*z
            acc[:, 0] += pol_sums(A[:, 0, 0], A[:, 0, 1], A[:, 0, 2])
            Hxx = cs * A[:, 1, 0] - sn * A[:, 2, 0]
            Hyy = cs * A[:, 1, 1] - sn * A[:, 2, 1]
            Hxy = cs * A[:, 1, 2] - sn * A[:, 2, 2]
            acc[:, 1] += pol_sums(Hxx, Hyy, Hxy)

    return acc, wsum


# --------------------------------------------------------------------------
# Outputs
# --------------------------------------------------------------------------
def write_header(f, b):
    f.write("********************************************************\n")
    f.write("****           Py_SHS: An open source software      ****\n")
    f.write("****         about Second Harmonic Scattering       ****\n")
    f.write("**** Dr Pierre-Marie GASSIN and Dr Gaelle GASSIN    ****\n")
    f.write("****                  ICGM - ENSCM                  ****\n")
    f.write("****                   june 2021                    ****\n")
    f.write("********************************************************\n")
    f.write("                                                        \n")
    f.write("           Don't forget to cite this work:              \n")
    f.write("    J. Chem. Inf. Model. 2020, 60, 12, pp 5912–5917     \n")
    f.write("                                                        \n")
    f.write("********************************************************\n")
    f.write("********************************************************\n")
    f.write("**       SIMULATION RESULTS  for cylinder programm  ****\n")
    f.write("********************************************************\n")
    f.write("*********          SURFACE KI(2) TENSOR      ***********\n")
    f.write("\n")
    f.write("    xxx=%f   xxy=%f   xxz= %f            \n"
            % (b[0, 0, 0], b[0, 0, 1], b[0, 0, 2]))
    f.write("    xyx=%f   xyy=%f   xyz= %f  \n"
            % (b[0, 1, 0], b[0, 1, 1], b[0, 1, 2]))
    f.write("    xzx=%f   xzy=%f   xzz= %f  \n"
            % (b[0, 2, 0], b[0, 2, 1], b[0, 2, 2]))
    f.write("\n")
    f.write("    yxx=%f   yxy=%f   yxz= %f  \n"
            % (b[1, 0, 0], b[1, 0, 1], b[1, 0, 2]))
    f.write("    yyx=%f   yyy=%f   yyz= %f  \n"
            % (b[1, 1, 0], b[1, 1, 1], b[1, 1, 2]))
    f.write("    yzx=%f   yzy=%f   yzz= %f  \n"
            % (b[1, 2, 0], b[1, 2, 1], b[1, 2, 2]))
    f.write("\n")
    f.write("    zxx=%f   zxy=%f   zxz= %f  \n"
            % (b[2, 0, 0], b[2, 0, 1], b[2, 0, 2]))
    f.write("    zyx=%f   zyy=%f   zyz= %f  \n"
            % (b[2, 1, 0], b[2, 1, 1], b[2, 1, 2]))
    f.write("    zzx=%f   zzy=%f   zzz= %f  \n"
            % (b[2, 2, 0], b[2, 2, 1], b[2, 2, 2]))
    f.write("\n")
    f.write("************************************************************************\n")


def report_polar(fout, title, V, H):
    """V, H are arrays (5,) = (a, b, c, b1, b2)."""
    av, bv, cv, bv1, bv2 = V
    ah, bh, ch, bh1, bh2 = H
    with np.errstate(all="ignore"):
        I4v = (av - bv + cv) / (3 * av + bv + 3 * cv)
        I2v = 4.0 * (av - cv) / (3 * av + bv + 3 * cv)
        I4h = (ah - bh + ch) / (3 * ah + bh + 3 * ch)
        I2h = 4.0 * (ah - ch) / (3 * ah + bh + 3 * ch)

    say("")
    say("**********calculation complete ************")
    say("av = %.6e" % av)
    say("bv = %.6e" % bv)
    say("bv1 = %.6e" % bv1)
    say("bv2 = %.6e" % bv2)
    say("cv = %.6e" % cv)
    say("ah = %.6e" % ah)
    say("bh = %.6e" % bh)
    say("bh1 = %.6e" % bh1)
    say("bh2 = %.6e" % bh2)
    say("ch = %.6e" % ch)
    say("********************************************")
    say("I2v = %f " % I2v)
    say("I4v = %f " % I4v)
    say("I2h = %f " % I2h)
    say("I4h = %f " % I4h)

    fout.write("********************************************\n")
    fout.write("                                            \n")
    fout.write("%s\n" % title)
    fout.write("                                            \n")
    fout.write("********************************************\n")
    fout.write("*************************\n")
    fout.write("I2v=%f  I4v=%f \n" % (I2v, I4v))
    fout.write("I2h=%f  I4h=%f \n" % (I2h, I4h))
    fout.write("*************************\n")
    fout.write("av=%.6e   bv=%.6e   bv1=%.6e  bv2=%.6e  cv= %.6e \n"
               % (av, bv, bv1, bv2, cv))
    fout.write("ah=%.6e   bh=%.6e   bh1=%.6e  bh2=%.6e  ch= %.6e \n"
               % (ah, bh, bh1, bh2, ch))

    try:
        fp = open(shs_io.path("out_plot"), "w", encoding="utf-8")
    except OSError:
        say("impossible d'ouvrir le fichier pour ecrire gnuplot!", end="")
        return
    with fp:
        fp.write("gamma IV(gamma) IH(gamma)\n")
        for ig in range(100):
            g = 2.0 * PI * ig / 100
            c, s = math.cos(g), math.sin(g)
            Iv = (av * c**4 + bv * c**2 * s**2 + cv * s**4
                  + bv1 * c**3 * s + bv2 * c * s**3)
            Ih = (ah * c**4 + bh * c**2 * s**2 + ch * s**4
                  + bh1 * c**3 * s + bh2 * c * s**3)
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

    say("********************************************************")
    say("****           Py_SHS: An open source software      ****")
    say("****         about Second Harmonic Scattering       ****")
    say("**** Dr Pierre-Marie GASSIN and Dr Gaelle GASSIN    ****")
    say("****                  ICGM - ENSCM                  ****")
    say("****                   june 2021                    ****")
    say("********************************************************")
    say("                cylinder SHS programm                     ")
    say("           Don't forget to cite this work:              ")
    say("    J. Chem. Inf. Model. 2020, 60, 12, pp 5912–5917     ")
    say("                                                        ")
    say("********************************************************")
    say("********       Enter the input parameters     **********")
    say("********************************************************")

    if len(argv) != 4:
        say("Computation failed: ")
        say("3 arguments expected : the keyword of the calculation, "
            "the Ki(2) inputfile and the outputfile. ")
        say("For more details see the documentation ")
        return 1

    keyword = argv[1]
    mode_single = keyword == "polarplot_single"
    mode_integ = keyword == "polarplot_integrate"
    mode_angle = keyword == "angle_scattering"
    if not (mode_single or mode_integ or mode_angle):
        say("--------------------- --------------------------------------------------")
        say("             ERROR in the name of the keyword                           ")
        say("USE ONLY polarplot_single OR polarplot_integrate OR angle_scattering  !!")
        say(" ! no calculation has been done ! Retry with the correct keyword        ")
        say("------------------------------------------------------------------------")
        return 1

    say("The input ki(2) file is %s,  and the output file is %s "
        % (argv[2], argv[3]))

    rd = TokenReader()
    try:
        say("Enter the wavelength in nm: ", end="")
        lam = rd.real()
        say("Enter REAL part of n_p (particle): ", end="")
        nr = rd.real()
        say("Enter IMAGINARY part of n_p (particle): ", end="")
        ni = rd.real()
        n_p = complex(nr, ni)
        say("Enter REAL part of n_m (medium): ", end="")
        nr = rd.real()
        say("Enter IMAGINARY part of n_m (medium): ", end="")
        ni = rd.real()
        n_m = complex(nr, ni)
        say("Enter the radius of the cylinder in nm: ", end="")
        R = rd.real()
        say("Enter the length of the cylinder in nm? ", end="")
        length = rd.real()
        say("      ******************************     ")
    except (EOFError, ValueError):
        return 1

    petitk = 2.0 * PI * n_p / lam
    grandk = 4.0 * PI * n_m / lam
    say("NLWKB vectors:")
    say("petitk (k_omega) = %.6f + i*%.6f" % (petitk.real, petitk.imag))
    say("grandk (k_2omega)= %.6f + i*%.6f" % (grandk.real, grandk.imag))
    say("")

    unit = shs_units.settings(grandk, "rad")
    rad_pref = unit.pref
    rad_pref2 = unit.multiplier
    for _l in unit.lines():
        say(_l)
    say("")

    beta = read_ki(argv[2])
    if beta is None:
        say("Problem to read the input ki(2) file! \n")
        return 1

    shs_io.set_output(argv[3])

    try:
        fout = open(argv[3], "w", encoding="utf-8")
    except OSError:
        say("problem to open the output file to write results!")
        return 1
    write_header(fout, beta)
    shs_units.write_header_lines(fout, unit)

    n_b, n_ph, n_th, grid_msg = shs_grid.cylinder_grid(
        R, length, petitk, grandk, (N_BETA, N_PHI, N_THETA))
    GRID = (n_b, n_ph, n_th)
    say("Grid: azimuth %d, orientation %d (phi) x %d (theta); axial integral analytic"
        % (n_b, n_ph, n_th))
    if grid_msg:
        say(grid_msg)

    # ---- modes ----
    if mode_single or mode_integ:
        try:
            if mode_single:
                say("Enter the scattered angle in degree (0°=transmission) ? ",
                    end="")
                t1 = rd.real()
                t2 = t1
                nang = 1
                say("      ******************************     ")
                say("********** Begin calculation polar plot at  %f ° ********"
                    % t1)
            else:
                say("Enter the beginning scattered angle in degree (0°=transmission) ? ",
                    end="")
                t1 = rd.real()
                say("      ******************************     ")
                say("Enter the final scattered angle in degree (0°=transmission) ? ",
                    end="")
                t2 = rd.real()
                say("      ******************************     ")
                say(" Number of integrale points ? ", end="")
                nang = rd.integer()
                say("      ******************************     ")
                if nang < 1:
                    return 1
        except (EOFError, ValueError):
            return 1

        if mode_single or nang == 1:
            th = np.array([t1 * 2.0 * PI / 360.0])
        else:
            th = (t1 + np.arange(nang) * (t2 - t1) / (nang - 1)) * 2.0 * PI / 360.0

        acc, wsum = compute(R, length, beta, petitk, grandk, th, *GRID)

        # Average over orientations and over scattering angles
        tot = acc.sum(axis=0)          # (2, 5)
        f = rad_pref2 / (wsum * nang)
        tot *= f

        title = ("* Cylinder POLAR PLOT - angle %g deg *" % t1
                 if mode_single
                 else "* Cylinder POLAR PLOT - average over %g ... %g deg (%d points) *"
                 % (t1, t2, nang))
        report_polar(fout, title, tot[0], tot[1])
        fout.close()

    else:
        # angle_scattering
        try:
            say("Enter the number of points between 0 and Pi for the angular "
                "SHS analysis ?", end="")
            pas2 = rd.integer()
            if pas2 < 1:
                return 1
        except (EOFError, ValueError):
            return 1

        nang = shs_angles.count(pas2)
        say("********** Begin Cylinder SHS Angular Distribution calculation**************")

        fout.write("***************************************************************\n")
        fout.write("                                            \n")
        fout.write("****            cylinder SHS angular distribution          ****\n")
        fout.write("                                            \n")
        fout.write(" ****************OUTPUT****************************************\n")
        fout.write(" grandtheta  IPSS-IH(0°)   IPPP-IH(90°)   ISSS-IV(0°)   ISPP-IV(90°) \n")

        th = shs_angles.scattering_angles(pas2)
        acc, wsum = compute(R, length, beta, petitk, grandk, th, *GRID)

        with open(shs_io.path("out_plot"), "w", encoding="utf-8") as fgr:
            fgr.write(" grandtheta   IH(0°)  IH(90°)  IV(0°)  IV(90°) \n")
            f = rad_pref2 / wsum
            for ia in range(nang):
                deg = th[ia] * 360.0 / (2.0 * PI)
                IH_0 = acc[ia, 1, 0] * f
                IH_90 = acc[ia, 1, 2] * f
                IV_0 = acc[ia, 0, 0] * f
                IV_90 = acc[ia, 0, 2] * f
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

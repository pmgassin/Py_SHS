#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
HRS.py  --  Py_SHS : Hyper-Rayleigh Scattering of a single molecule (Python / NumPy)

Derived from SHS.py for the special case N = 1, position (0,0,0),
molecular orientation (phi, theta, psi) = (0,0,0).

No orientation file is required: the molecular hyperpolarizability beta
is used directly.  The isotropic orientational average is performed over
the laboratory-frame Euler angles of the molecule (same grid as SHS).

Because r = 0 the phase factor is identically 1, so the coherent sum
reduces to the single-molecule tensor.  Numerically this recovers exactly
the result of SHS.py run with N = 1 and a one-line orientation file
"0 0 0 0 0 0".

Keywords : polarplot_single | polarplot_integrate | angle_scattering
Usage    : python HRS.py keyword beta_file output_file [--legacy]
stdin    : wavelength (nm), Re(n), Im(n), then the angle(s) of the mode
           (the wavelength and the index are read to keep the same input protocol as the
           other modules; they do not change the default `raw` intensities -- they only
           enter the radiative prefactor (|K|^2/4pi)^2 when PYSHS_UNITS=rad is requested)

Options / environment
---------------------
--legacy  (or PYSHS_LEGACY=1) reproduces two defects of the old C code: the
          normalisation of bv1, bv2, bh1, bh2 (not divided by 8 pi^2) and the
          pairing of the sin*cos cross terms of bh1, bh2.  After the isotropic
          orientational average these four coefficients are identically zero (the
          odd terms are the real part of a Hermitian-symmetric contraction that the
          chiral invariants cannot feed), so --legacy only changes rounding noise in
          them and none of av, bv, cv, ah, bh, ch, I2, I4.  Kept for completeness;
          see tests/test_core.py (TestLegacyMode) for the exact definitions.
"""
import os
import sys
import math

import numpy as np

from pyshs import shs_units
from pyshs import shs_angles
from pyshs import shs_io

# ---- orientation grid (same values as SHS.py / the C code) ----
PAS = 10
N_PHI = 2 * PAS          # phi   in [0, 2 pi)
N_THETA = PAS            # theta in [0, pi]
N_PSI = 2 * PAS          # psi   in [0, 2 pi)
PI = math.pi


# --------------------------------------------------------------------------
# Small helpers (identical to SHS.py)
# --------------------------------------------------------------------------
# number format of the intensities: historical '%f' when the multiplier is 1,
# scientific notation otherwise (values can be ~1e-9 with PYSHS_UNITS=rad)
FMT = "%f"


class TokenReader:
    """Reads whitespace separated values from stdin like C scanf."""

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
    """ZYZ Euler matrices, identical to Txx..Tzz of the C code.
    phi, theta, psi : arrays of shape (n,)  ->  T : (n, 3, 3)"""
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


_JK = ((0, 0), (1, 1), (0, 1))


def rotate_needed(T, S):
    """L[n,i,q] = sum_a T[n,i,a] * ( sum_bc T[n,j,b] T[n,k,c] S[n,a,b,c] ),
    q <-> (j,k) in ((0,0),(1,1),(0,1))."""
    L = np.empty(S.shape[:1] + (3, 3), dtype=S.dtype)
    for q, (j, k) in enumerate(_JK):
        w = np.einsum("nb,nc,nabc->na", T[:, j, :], T[:, k, :], S)
        L[:, :, q] = np.einsum("nia,na->ni", T, w)
    return L


def pol_sums(w, xx, yy, xy):
    """Weighted sums (a, b, c, b1, b2) of the intensity coefficients."""
    xyc = np.conj(xy)
    return np.array([
        np.sum(w * (xx * np.conj(xx)).real),
        np.sum(w * (4.0 * (xy * xyc).real + 2.0 * (xx * np.conj(yy)).real)),
        np.sum(w * (yy * np.conj(yy)).real),
        np.sum(w * 4.0 * (xx * xyc).real),
        np.sum(w * 4.0 * (yy * xyc).real),
    ])


# --------------------------------------------------------------------------
# Inputs
# --------------------------------------------------------------------------
def read_beta(name):
    """1 header line + 27 lines "name value"."""
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
# Core calculation (HRS = single dipole at origin, identity orientation)
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



class Orientations:
    """Orientation table (Euler matrices and weights), computed ONCE."""

    def __init__(self):
        hx, hz = 2 * PI / N_PHI, 2 * PI / N_PSI
        th_nodes, th_w = theta_rule(N_THETA)
        i, j, k = np.meshgrid(np.arange(N_PHI), np.arange(N_THETA),
                              np.arange(N_PSI), indexing="ij")
        i, j, k = i.ravel(), j.ravel(), k.ravel()
        self.T = euler_matrices(0.5 * hx + i * hx, th_nodes[j], 0.5 * hz + k * hz)
        self.w = hx * hz * th_w[j]
        self.n = i.size


class AngleTable:
    """Angle dependent constants, computed ONCE."""

    def __init__(self, theta_s, petitk, grandk):
        theta_s = np.asarray(theta_s, dtype=np.float64)
        self.cs, self.sn = np.cos(theta_s), np.sin(theta_s)
        # q is not needed for HRS (r = 0 => phase = 1), but cs/sn are used
        # in finalize for the H polarisation combination.


def lab_amplitudes(S, ori):
    """Lab-frame amplitudes A[n, i, q] (q = xx, yy, xy) for every orientation.
    They do not depend on the scattering angle, so they are computed once.
    S is real for HRS (phase = 1; the imaginary half S[:, 27:] is zero)."""
    n = ori.n
    Lr = rotate_needed(ori.T, S[:, :27].reshape(n, 3, 3, 3))
    return Lr + 0j


def combine(A, w, c, s, legacy):
    """Polarisation sums and orientational average for one angle (cos c, sin s).
    Returns the (2, 5) array [V; H] of coefficients (a, b, c, b1, b2), unnormalised."""
    V = pol_sums(w, A[:, 0, 0], A[:, 0, 1], A[:, 0, 2])
    Hxx = c * A[:, 1, 0] - s * A[:, 2, 0]
    Hyy = c * A[:, 1, 1] - s * A[:, 2, 1]
    Hxy = c * A[:, 1, 2] - s * A[:, 2, 2]
    H = pol_sums(w, Hxx, Hyy, Hxy)
    if legacy:
        gyx, gzx, gyy, gzy = A[:, 1, 0], A[:, 2, 0], A[:, 1, 1], A[:, 2, 1]
        gyxy, gzxy = A[:, 1, 2], A[:, 2, 2]

        def r4(u, v):
            return 4.0 * (u * np.conj(v)).real

        H[3] += np.sum(w * (-r4(Hxx, Hxy)))
        H[3] += np.sum(w * (c * c * r4(gyx, gyxy) + s * s * r4(gzx, gzxy)
                            - s * c * (r4(gzx, gzxy) + r4(gyx, gyxy))))
        H[4] += np.sum(w * (-r4(Hyy, Hxy)))
        H[4] += np.sum(w * (c * c * r4(gyy, gyxy) + s * s * r4(gzy, gzxy)
                            - s * c * (r4(gzy, gzxy) + r4(gyy, gyxy))))
    return np.vstack((V, H))


def finalize(S, ori, ang, ag, legacy):
    """Rotation to the lab frame and orientational average for angle ag
    (kept as one call for the tests; compute() rotates only once)."""
    return combine(lab_amplitudes(S, ori), ori.w, ang.cs[ag], ang.sn[ag], legacy)


def compute(beta, petitk, grandk, theta_s, legacy):
    """HRS calculation for a single molecule.
    beta is already the molecular tensor (meso = mol for identity orientation).
    Returns acc[nang, 2, 5]."""
    nang = len(theta_s)
    ori = Orientations()
    ang = AngleTable(theta_s, petitk, grandk)
    acc = np.zeros((nang, 2, 5))

    # Single dipole at origin with identity orientation: the phase factor is 1, so the
    # lab-frame amplitudes are the same for every angle; only the H combination
    # (cos Theta * Y - sin Theta * Z) depends on the angle.
    S = np.zeros((ori.n, 54))
    S[:, :27] = beta.reshape(27)        # same tensor for all orientations
    A = lab_amplitudes(S, ori)
    next_print = 0
    for ag in range(nang):
        acc[ag] = combine(A, ori.w, ang.cs[ag], ang.sn[ag], legacy)
        pct = int(100.0 * (ag + 1) / nang)
        if pct >= next_print:
            say("calculation in progress : %f pourcent " % pct)
            next_print = pct + 5
    say("calculation in progress : %f pourcent " % 100.0)
    return acc


# --------------------------------------------------------------------------
# Outputs (identical format to SHS.py)
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
    f.write("******          SIMULATION RESULTS               *******\n")
    f.write("********************************************************\n")
    f.write("******     SURFACE HYPERPOLARISABILITY TENSOR     ******\n")
    f.write("\n")
    f.write("    xxx=%f   xxy=%f   xxz= %f            \n" % (b[0, 0, 0], b[0, 0, 1], b[0, 0, 2]))
    f.write("    xyx=%f   xyy=%f   xyz= %f  \n" % (b[0, 1, 0], b[0, 1, 1], b[0, 1, 2]))
    f.write("    xzx=%f   xzy=%f   xzz= %f  \n" % (b[0, 2, 0], b[0, 2, 1], b[0, 2, 2]))
    f.write("\n")
    f.write("    yxx=%f   yxy=%f   yxz= %f  \n" % (b[1, 0, 0], b[1, 0, 1], b[1, 0, 2]))
    f.write("    yyx=%f   yyy=%f   yyz= %f  \n" % (b[1, 1, 0], b[1, 1, 1], b[1, 1, 2]))
    f.write("    yzx=%f   yzy=%f   yzz= %f  \n" % (b[1, 2, 0], b[1, 2, 1], b[1, 2, 2]))
    f.write("\n")
    f.write("    zxx=%f   zxy=%f   zxz= %f  \n" % (b[2, 0, 0], b[2, 0, 1], b[2, 0, 2]))
    f.write("    zyx=%f   zyy=%f   zyz= %f  \n" % (b[2, 1, 0], b[2, 1, 1], b[2, 1, 2]))
    f.write("    zzx=%f   zzy=%f   zzz= %f  \n" % (b[2, 2, 0], b[2, 2, 1], b[2, 2, 2]))
    f.write("\n")
    f.write("*******************************************************\n")


def report_polar(fout, deg, ndip, V, H):
    """Final report of polarplot_single / polarplot_integrate."""
    av, bv, cv, bv1, bv2 = V[0], V[1], V[2], V[3], V[4]
    ah, bh, ch, bh1, bh2 = H[0], H[1], H[2], H[3], H[4]
    with np.errstate(all="ignore"):
        I4v = (av - bv + cv) / (3 * av + bv + 3 * cv)
        I2v = 4.0 * (av - cv) / (3 * av + bv + 3 * cv)
        I4h = (ah - bh + ch) / (3 * ah + bh + 3 * ch)
        I2h = 4.0 * (ah - ch) / (3 * ah + bh + 3 * ch)

    say("")
    say("**********calculation complete ************")
    for name, val in (("av", av), ("bv", bv), ("bv1", bv1), ("bv2", bv2), ("cv", cv),
                      ("ah", ah), ("bh", bh), ("bh1", bh1), ("bh2", bh2), ("ch", ch)):
        say(("%s = " + FMT) % (name, val))
    say("********************************************")
    say("I2v = %f " % I2v)
    say("I4v = %f " % I4v)
    say("I2h = %f " % I2h)
    say("I4h = %f " % I4h)

    fout.write("                                      \n")
    fout.write(" POLAR PLOT - CONFIGURATION AT %f °  \n" % deg)
    fout.write("                                      \n")
    fout.write("nombre de dipole = %d  \n" % ndip)
    fout.write("*************************\n")
    fout.write("I2v=%f  I4v=%f \n" % (I2v, I4v))
    fout.write("I2h=%f  I4h=%f \n" % (I2h, I4h))
    fout.write("*************************\n")
    fout.write("av=%f   bv=%f   bv1=%f  bv2=%f  cv= %f \n".replace("%f", FMT) % (av, bv, bv1, bv2, cv))
    fout.write("ah=%f   bh=%f   bh1=%f  bh2=%f  ch= %f \n".replace("%f", FMT) % (ah, bh, bh1, bh2, ch))

    try:
        fp = open(shs_io.path("out_plot"), "w", encoding="utf-8")
    except OSError:
        say("impossible d'ouvrir le fichier pour ecrire gnuplot!", end="")
        return
    with fp:
        fp.write("gamma IV(gamma) IH(gamma)\n")
        g = 2.0 * PI * np.arange(100) / 100
        c, s = np.cos(g), np.sin(g)
        Iv = av * c**4 + bv * c**2 * s**2 + cv * s**4 + bv1 * c**3 * s + bv2 * c * s**3
        Ih = ah * c**4 + bh * c**2 * s**2 + ch * s**4 + bh1 * c**3 * s + bh2 * c * s**3
        for gg, a, b in zip(g, Iv, Ih):
            fp.write(("%f " + FMT + " " + FMT + "\n") % (gg, a, b))


# --------------------------------------------------------------------------
def main(argv):
    legacy = os.environ.get("PYSHS_LEGACY", "") not in ("", "0")
    if "--legacy" in argv:
        legacy = True
        argv = [a for a in argv if a != "--legacy"]
    if legacy:
        print("NOTE: legacy mode (old C quirks) only affects bv1, bv2, bh1, bh2, which are "
              "zero after the isotropic average; av, bv, cv, ah, bh, ch are unchanged.")
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
    say("                                                        ")
    say("           Don't forget to cite this work:              ")
    say("    J. Chem. Inf. Model. 2020, 60, 12, pp 5912–5917     ")
    say("                                                        ")
    say("********************************************************")
    say("********************************************************")
    say("********       Enter the input parameters     **********")
    say("********************************************************")

    # HRS: 3 arguments (keyword, beta_file, output_file) — no orientation file
    if len(argv) != 4:
        say("Computation failed: ")
        say("3 arguments expected : the keyword of the calculation, "
            "the beta inputfile and the outputfile. ")
        say("For more details see the documentation ")
        return 1

    keyword = argv[1]
    mode_single = keyword == "polarplot_single"
    mode_integ = keyword == "polarplot_integrate"
    mode_angle = keyword == "angle_scattering"
    if not (mode_single or mode_integ or mode_angle):
        say("--------------------- ----------------------------------------------------")
        say("             ERROR in the name of the keyword                             ")
        say(" USE ONLY polarplot_single OR polarplot_integrate OR angle_scattering  !! ")
        say(" ! no calculation has been done ! Retry with the correct keyword          ")
        say("--------------------------------------------------------------------------")
        return 1

    rd = TokenReader()
    try:
        # HRS: no "Number of dipoles" — always N = 1
        ndip = 1
        say("HRS mode: single dipole at origin (N = 1)")
        say("Enter the wavelength in nm? ", end="")
        lam = rd.real()
        say("The wavelength is: %f " % lam)
        say("Enter the real part of the refractive index? ", end="")
        nr = rd.real()
        say("Enter the imaginary part of the refractive index? ", end="")
        ni = rd.real()
    except (EOFError, ValueError):
        return 1

    n = complex(nr, ni)
    say("The refractive index is: %f + I*%f " % (n.real, n.imag))
    say("The input beta file is %s and the output file is %s "
        % (argv[2], argv[3]))
    grandk = 4.0 * PI * n / lam
    petitk = 2.0 * PI * n / lam

    beta = read_beta(argv[2])
    if beta is None:
        say("Problem to read the input beta file \n!")
        return 1

    shs_io.set_output(argv[3])

    try:
        fout = open(argv[3], "w", encoding="utf-8")
    except OSError:
        say("problem to open the output file to write results!")
        return 1

    unit = shs_units.settings(grandk, "raw")
    for _l in unit.lines():
        say(_l)
    write_header(fout, beta)
    shs_units.write_header_lines(fout, unit)

    # ---------- angles ----------
    t1 = t2 = 0.0
    nang = 1
    try:
        if mode_single:
            say("Enter the scattered angle in degree (0°=transmission) ? ", end="")
            t1 = rd.real()
            say("      ******************************     ")
            say("********** Begin calculation polar plot at  %f ° ********" % t1)
        elif mode_integ:
            say("Enter the beginning scattered angle in degree (0°=transmission) ? ", end="")
            t1 = rd.real()
            say("      ******************************     ")
            say("Enter the final scattered angle in degree (0°=transmission) ? ", end="")
            t2 = rd.real()
            say("      ******************************     ")
            say(" Number of integrale points ? ", end="")
            nang = rd.integer()
            say("      ******************************     ")
            if nang < 1:
                return 1
        else:
            fout.write("                                      \n")
            fout.write("* HRS angular distribution            \n")
            fout.write("                                      \n")
            fout.write("nombre de dipole = %d  \n" % ndip)
            fout.write(" ****************OUTPUT****************************************\n")
            fout.write(" grandtheta  IPSS-IH(0°)   IPPP-IH(90°)   ISSS-IV(0°)   ISPP-IV(90°) \n")
            say("Enter the number of points between 0 and 180° for the angular HRS analysis ?", end="")
            pas2 = rd.integer()
            if pas2 < 1:
                return 1
            nang = shs_angles.count(pas2)
            say("********** Begin HRS Angular Distribution calculation**************")
    except (EOFError, ValueError):
        return 1

    ia = np.arange(nang)
    if mode_angle:
        th = shs_angles.scattering_angles(pas2)
    elif mode_integ:
        th = (np.full(nang, t1) if nang == 1
              else t1 + ia * (t2 - t1) / (nang - 1)) * 2.0 * PI / 360.0
    else:
        th = np.full(1, t1 * 2.0 * PI / 360.0)

    acc = compute(beta, petitk, grandk, th, legacy)

    norm = unit.multiplier / (8.0 * PI * PI)      # units: see shs_units.py
    global FMT
    if unit.multiplier != 1.0:
        FMT = "%.6e"
    if mode_angle:
        with open(shs_io.path("out_plot"), "w", encoding="utf-8") as fgr:
            fgr.write(" grandtheta   IH(0°)  IH(90°)  IV(0°)  IV(90°) \n")
            for k in range(nang):
                deg = th[k] * 360.0 / (2.0 * PI)
                IH_0, IH_90 = acc[k, 1, 0] * norm, acc[k, 1, 2] * norm
                IV_0, IV_90 = acc[k, 0, 0] * norm, acc[k, 0, 2] * norm
                say("angle value: %f " % deg)
                say("IH_0=%f    IH_90=%f   IV_0=%f    IV_90=%f ".replace("%f", FMT)
                    % (IH_0, IH_90, IV_0, IV_90))
                row = (" %f    %f   %f   %f   %f \n".replace("%f", FMT).replace(FMT, "%f", 1)) % (deg, IH_0, IH_90, IV_0, IV_90)
                fout.write(row)
                fgr.write(row)
        say("calculation finished ")
    else:
        tot = acc.sum(axis=0)                      # (2, 5)
        f = norm / nang
        if legacy:
            tot[:, :3] *= f
            tot[:, 3:] *= unit.multiplier / nang
        else:
            tot *= f
        report_polar(fout, th[nang - 1] * 360.0 / (2.0 * PI), ndip,
                     tot[0], tot[1])
    fout.close()
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))

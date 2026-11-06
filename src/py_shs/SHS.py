#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
SHS.py  --  Py_SHS : SHS of an assembly of N oriented dipoles (Python / NumPy)

Python version of SHS.c: same physics and same input protocol (see the README for the
differences in numerical defaults: Gauss-Legendre orientation grid, automatic grid size,
units, angle grid).  All the dipoles are loaded in memory (about 240 bytes per dipole), so
the memory grows linearly with N.

Method
------
The rotation of the tensors by the orientation (phi, theta, psi) is LINEAR,
so it commutes with the sum over the dipoles:

    A_lab[i][jk] = sum_abc T_ia T_jb T_kc  S_abc
    S_abc        = sum_dipoles  exp(i q.r_d) * P_d[abc]

In NumPy the sum over the dipoles is a matrix product
    S[orientation, 27] = phase[orientation, dipole] @ P[dipole, 27]
(done by BLAS), and the rotation to the laboratory frame is done once per
(orientation, angle).  The orientations are independent: they are cut in
chunks that are processed by a pool of threads (NumPy releases the GIL).

Keywords : polarplot_single | polarplot_integrate | angle_scattering
Usage    : python SHS.py keyword beta_file orientation_file output_file [--legacy]
stdin    : N, wavelength (nm), Re(n), Im(n), then the angle(s) of the mode

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
PYSHS_THREADS=k  number of threads (default: number of CPU cores).
"""
import os

# Our own thread pool provides the parallelism: avoid BLAS oversubscription.
for _v in ("OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "OMP_NUM_THREADS",
           "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")

import sys
import math
from concurrent.futures import ThreadPoolExecutor, as_completed

import numpy as np

from py_shs import shs_grid
from py_shs import shs_units
from py_shs import shs_angles
from py_shs import shs_io

# ---- orientation grid of the whole assembly (same values as the C code) ----
PAS = 10
N_PHI = 2 * PAS          # phi   in [0, 2 pi)
N_THETA = PAS            # theta in [0, pi]
N_PSI = 2 * PAS          # psi   in [0, 2 pi)
PI = math.pi

BLOCK = 4096             # dipoles per block of the phase matrix
TARGET_ELEMS = 1 << 20   # elements of one (orientations x dipoles) work array


# --------------------------------------------------------------------------
# Small helpers
# --------------------------------------------------------------------------
# number format of the intensities: historical '%f' when the multiplier is 1,
# scientific notation otherwise (values can be ~1e-9 with PYSHS_UNITS=rad)
FMT = "%f"


class TokenReader:
    """Reads whitespace separated values from stdin like C scanf (works for
    piped input as well as for interactive use)."""

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


def rotate_tensors(T, beta):
    """out[n,i,j,k] = sum_abc T[n,i,a] T[n,j,b] T[n,k,c] beta[a,b,c]
    (one rotated copy of the 3x3x3 tensor beta per Euler matrix)."""
    t1 = np.einsum("nia,abc->nibc", T, beta)
    t2 = np.einsum("njb,nibc->nijc", T, t1)
    return np.einsum("nkc,nijc->nijk", T, t2)


# the 3 lab polarisation components that are used: (j,k) = xx, yy, xy
_JK = ((0, 0), (1, 1), (0, 1))


def rotate_needed(T, S):
    """L[n,i,q] = sum_a T[n,i,a] * ( sum_bc T[n,j,b] T[n,k,c] S[n,a,b,c] ),
    q <-> (j,k) in ((0,0),(1,1),(0,1)).  Works for real or complex S."""
    L = np.empty(S.shape[:1] + (3, 3), dtype=S.dtype)
    for q, (j, k) in enumerate(_JK):
        w = np.einsum("nb,nc,nabc->na", T[:, j, :], T[:, k, :], S)
        L[:, :, q] = np.einsum("nia,na->ni", T, w)
    return L


def pol_sums(w, xx, yy, xy):
    """Weighted sums (a, b, c, b1, b2) of the intensity coefficients
        I(gamma) = a c^4 + b c^2 s^2 + c' s^4 + b1 c^3 s + b2 c s^3
    for one polarisation, from the complex amplitudes (Xxx, Xyy, Xxy)."""
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


def parse_lines(lines, first_index):
    """Parses dipole lines (phi theta psi x y z ...) -> (n, 6) array.
    Raises ValueError("<dipole number>") on a malformed line."""
    n = len(lines)
    try:                                    # fast path: exactly 6 numbers per line
        tok = " ".join(lines).split()
        if len(tok) == 6 * n:
            return np.array(tok, dtype=np.float64).reshape(n, 6)
    except ValueError:
        pass
    out = np.empty((n, 6))                  # slow path: like sscanf, extra columns ignored
    for i, line in enumerate(lines):
        try:
            tok = line.split()
            if len(tok) < 6:
                raise ValueError
            out[i] = [float(t) for t in tok[:6]]
        except ValueError:
            raise ValueError(str(first_index + i + 1))
    return out


def write_arrows(farr, data):
    """out_arrows_gnuplot : x y z and the unit vector of each dipole."""
    phi, theta = data[:, 0], data[:, 1]
    arr = np.column_stack((data[:, 3], data[:, 4], data[:, 5],
                           np.cos(phi) * np.sin(theta),
                           np.sin(phi) * np.sin(theta),
                           np.cos(theta)))
    np.savetxt(farr, arr, fmt="%f", delimiter=" ", newline=" \n")


def load_dipoles(name, ndip, beta, farr):
    """Reads the N dipoles; returns the tensors P (N,27) in the assembly frame
    and the positions (N,3)."""
    P = np.empty((ndip, 27))
    pos = np.empty((ndip, 3))
    try:
        f = open(name, "rt")
    except OSError:
        say("Problem to open the input orientation file \n!")
        return None
    with f:
        d = 0
        while d < ndip:
            n = min(BLOCK, ndip - d)
            lines = []
            for _ in range(n):
                line = f.readline()
                if line == "":
                    break
                lines.append(line)
            try:
                if len(lines) < n:
                    raise ValueError(str(d + len(lines) + 1))
                data = parse_lines(lines, d)
            except ValueError as e:
                sys.stderr.write("Problem reading dipole %s of the orientation file\n" % e)
                return None
            if farr is not None:
                write_arrows(farr, data)
            Td = euler_matrices(data[:, 0], data[:, 1], data[:, 2])
            P[d:d + n] = rotate_tensors(Td, beta).reshape(n, 27)
            pos[d:d + n] = data[:, 3:6]
            d += n
        # the N entered by the user is authoritative: warn if the file holds more dipoles
        extra = 0
        for line in f:
            if line.strip():
                extra += 1
        if extra:
            say("WARNING: the orientation file contains %d more line(s) than the %d dipole(s) "
                "requested: they are ignored. Check the number of dipoles you entered." % (extra, ndip))
    return P, pos


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



class Orientations:
    """Orientation table (Euler matrices and weights), computed ONCE."""

    def __init__(self, grid=None):
        n_phi, n_theta, n_psi = grid if grid is not None else (N_PHI, N_THETA, N_PSI)
        hx, hz = 2 * PI / n_phi, 2 * PI / n_psi
        th_nodes, th_w = theta_rule(n_theta)
        i, j, k = np.meshgrid(np.arange(n_phi), np.arange(n_theta),
                              np.arange(n_psi), indexing="ij")
        i, j, k = i.ravel(), j.ravel(), k.ravel()
        self.T = euler_matrices(0.5 * hx + i * hx, th_nodes[j], 0.5 * hz + k * hz)
        self.w = hx * hz * th_w[j]
        self.n = i.size


class AngleTable:
    """Angle dependent constants, computed ONCE."""

    def __init__(self, theta_s, petitk, grandk):
        theta_s = np.asarray(theta_s, dtype=np.float64)
        self.cs, self.sn = np.cos(theta_s), np.sin(theta_s)
        qz = 2.0 * petitk - grandk * self.cs          # z-term of the phase
        qy = -grandk * self.sn                        # y-term of the phase
        self.qzr, self.qzi = qz.real.copy(), qz.imag.copy()
        self.qyr, self.qyi = qy.real.copy(), qy.imag.copy()


def accumulate(S, T, ang, ag, pos_blk, P_blk):
    """S (n, 54) += sum over the dipoles of the block of
    exp(i kv.r_d) * P_d,  kv = qz*T[2,:] + qy*T[1,:]  (complex 3-vector).
    S[:, :27] is the real part, S[:, 27:] the imaginary part."""
    kr = ang.qzr[ag] * T[:, 2, :] + ang.qyr[ag] * T[:, 1, :]
    ar = kr @ pos_blk.T                                   # (n, B)
    pr = np.cos(ar)
    pi_ = np.sin(ar)
    if ang.qzi[ag] != 0.0 or ang.qyi[ag] != 0.0:          # absorbing medium
        ki = ang.qzi[ag] * T[:, 2, :] + ang.qyi[ag] * T[:, 1, :]
        e = np.exp(-(ki @ pos_blk.T))
        pr *= e
        pi_ *= e
    S[:, :27] += pr @ P_blk
    S[:, 27:] += pi_ @ P_blk


def finalize(S, ori, ang, ag, legacy):
    """Rotation to the lab frame (once per orientation) and orientation
    average for the angle ag.  Returns the (2, 5) array [V; H] of the
    coefficients (a, b, c, b1, b2), still unnormalised."""
    n = ori.n
    T, w = ori.T, ori.w
    Lr = rotate_needed(T, S[:, :27].reshape(n, 3, 3, 3))
    Li = rotate_needed(T, S[:, 27:].reshape(n, 3, 3, 3))
    A = Lr + 1j * Li                                      # A[n, i, q]
    c, s = ang.cs[ag], ang.sn[ag]
    V = pol_sums(w, A[:, 0, 0], A[:, 0, 1], A[:, 0, 2])   # lab x component
    Hxx = c * A[:, 1, 0] - s * A[:, 2, 0]                 # cos(Theta) y - sin(Theta) z
    Hyy = c * A[:, 1, 1] - s * A[:, 2, 1]
    Hxy = c * A[:, 1, 2] - s * A[:, 2, 2]
    H = pol_sums(w, Hxx, Hyy, Hxy)
    if legacy:
        # cross terms of bh1/bh2 as written in the old code (they used
        # gyxx*gyxy and gzxx*gzxy instead of the mixed products)
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


def compute(P, pos, petitk, grandk, theta_s, nthreads, legacy, grid=None):
    """Full calculation.  Returns acc[nang, 2, 5]."""
    ndip, nang = P.shape[0], len(theta_s)
    ori = Orientations(grid)
    ang = AngleTable(theta_s, petitk, grandk)
    acc = np.zeros((nang, 2, 5))

    blk = max(1, min(BLOCK, ndip))
    nochunk = max(1, min(ori.n, TARGET_ELEMS // blk))
    chunks = [(o0, min(o0 + nochunk, ori.n)) for o0 in range(0, ori.n, nochunk)]
    say("Python threads: %d" % nthreads)

    def work(S, o0, o1, ag):
        T = ori.T[o0:o1]
        for d0 in range(0, ndip, blk):
            d1 = min(d0 + blk, ndip)
            accumulate(S[o0:o1], T, ang, ag, pos[d0:d1], P[d0:d1])

    total = nang * len(chunks)
    done, next_print = 0, 0
    with ThreadPoolExecutor(max_workers=nthreads) as pool:
        for ag in range(nang):
            S = np.zeros((ori.n, 54))
            futures = [pool.submit(work, S, o0, o1, ag) for o0, o1 in chunks]
            for fut in as_completed(futures):
                fut.result()
                pct = int(100.0 * done / total)
                done += 1
                if pct >= next_print:
                    say("calculation in progress : %f pourcent " % pct)
                    next_print = pct + 5
            acc[ag] = finalize(S, ori, ang, ag, legacy)
    say("calculation in progress : %f pourcent " % 100.0)
    return acc


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

    if len(argv) != 5:
        say("Computation failed: ")
        say("4 arguments expected : the keyword of the calculation, the beta and orientation inputfiles and the outputfile. ")
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
        say("Number of dipoles ? ", end="")
        ndip = rd.integer()
        if ndip <= 0:
            sys.stderr.write("invalid number of dipoles\n")
            return 1
        say("The number of dipoles is: %d " % ndip)
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
    say("The input beta file is %s, the input orientation file is %s and the output file is %s "
        % (argv[2], argv[3], argv[4]))
    grandk = 4.0 * PI * n / lam
    petitk = 2.0 * PI * n / lam

    beta = read_beta(argv[2])
    if beta is None:
        say("Problem to read the input beta file \n!")
        return 1
    shs_io.set_output(argv[4])
    try:
        fout = open(argv[4], "w", encoding="utf-8")
    except OSError:
        say("problem to open the output file to write results!")
        return 1
    try:
        farr = open(shs_io.path("out_arrows_gnuplot"), "w", encoding="utf-8")
    except OSError:
        farr = None
        say("impossible d'ouvrir le fichier pour ecrire fleche gnuplot!", end="")
    loaded = load_dipoles(argv[3], ndip, beta, farr)
    if farr is not None:
        farr.close()
    if loaded is None:
        return 1
    P, pos = loaded
    unit = shs_units.settings(grandk, "raw")
    for _l in unit.lines():
        say(_l)
    write_header(fout, beta)
    shs_units.write_header_lines(fout, unit)

    # ---------- angles ----------
    nthreads = int(os.environ.get("PYSHS_THREADS", "0") or 0)
    if nthreads <= 0:
        try:
            nthreads = len(os.sched_getaffinity(0))
        except AttributeError:
            nthreads = os.cpu_count() or 1
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
            fout.write("* SHS angular distribution            \n")
            fout.write("                                      \n")
            fout.write("nombre de dipole = %d  \n" % ndip)
            fout.write(" ****************OUTPUT****************************************\n")
            fout.write(" grandtheta  IPSS-IH(0°)   IPPP-IH(90°)   ISSS-IV(0°)   ISPP-IV(90°) \n")
            say("Enter the number of points between 0 and 180° for the angular SHS analysis ?", end="")
            pas2 = rd.integer()
            if pas2 < 1:
                return 1
            nang = shs_angles.count(pas2)
            say("********** Begin SHS Angular Distribution calculation**************")
    except (EOFError, ValueError):
        return 1

    ia = np.arange(nang)
    if mode_angle:
        th = shs_angles.scattering_angles(pas2)
    elif mode_integ:
        th = (np.full(nang, t1) if nang == 1 else t1 + ia * (t2 - t1) / (nang - 1)) * 2.0 * PI / 360.0
    else:
        th = np.full(1, t1 * 2.0 * PI / 360.0)

    rmax = float(np.max(np.linalg.norm(pos, axis=1))) if len(pos) else 0.0
    n_ph, n_th, n_ps, grid_msg = shs_grid.dipole_grid(
        rmax, petitk, grandk, (N_PHI, N_THETA, N_PSI))
    say("Orientation grid: %d x %d x %d  (max |r| = %.1f, q*rmax = %.1f)"
        % (n_ph, n_th, n_ps, rmax, shs_grid.qmax(petitk, grandk) * rmax))
    if grid_msg:
        say(grid_msg)
    acc = compute(P, pos, petitk, grandk, th, nthreads, legacy, (n_ph, n_th, n_ps))

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
                say("IH_0=%f    IH_90=%f   IV_0=%f    IV_90=%f ".replace("%f", FMT) % (IH_0, IH_90, IV_0, IV_90))
                row = (" %f    %f   %f   %f   %f \n".replace("%f", FMT).replace(FMT, "%f", 1)) % (deg, IH_0, IH_90, IV_0, IV_90)
                fout.write(row)
                fgr.write(row)
        say("calculation finished ")
    else:
        # average over the scattering angles, then 1/(8 pi^2) for the orientations
        tot = acc.sum(axis=0)                      # (2, 5)
        f = norm / nang
        if legacy:
            # old code: bv1, bv2, bh1, bh2 were NOT divided by 8 pi^2
            tot[:, :3] *= f
            tot[:, 3:] *= unit.multiplier / nang
        else:
            tot *= f
        report_polar(fout, th[nang - 1] * 360.0 / (2.0 * PI), ndip, tot[0], tot[1])
    fout.close()
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))

"""Automatic quadrature-grid sizing for the PySHS modules.

All the intensity integrals contain the phase factor exp(i q.r) with
|q| <= qmax = 2 Re(k_w) + Re(K_2w).  A grid that is fine enough for a small object
silently becomes under-resolved for a large one, so the grids are sized from the
dimensionless "size parameter"  a = qmax * (largest half-extent of the object).

Environment variables
---------------------
PYSHS_GRID        auto (default) | fixed
                  fixed: use the historical hard-coded grid sizes of each module.
PYSHS_GRID_SCALE  float (default 1.0): multiplies the size parameter used by the rules
                  (e.g. 2 for a safer / slower grid, 0.5 for a faster one).
PYSHS_GRID_MAX    int (default depends on module): upper bound on the grid size, to
                  keep memory and run time under control. A warning is issued when the
                  bound prevents the rule from being satisfied.
"""
import math
import os


def auto_enabled():
    return os.environ.get("PYSHS_GRID", "auto").strip().lower() != "fixed"


def scale():
    try:
        s = float(os.environ.get("PYSHS_GRID_SCALE", "1"))
        return s if s > 0 else 1.0
    except ValueError:
        return 1.0


def grid_max(default):
    try:
        return int(os.environ.get("PYSHS_GRID_MAX", default))
    except ValueError:
        return default


def qmax(petitk, grandk):
    """Largest modulus of the scattering vector q = 2 k_w - K_2w (real part)."""
    return 2.0 * abs(complex(petitk).real) + abs(complex(grandk).real)


def size_parameter(petitk, grandk, extent):
    return qmax(petitk, grandk) * extent * scale()


# --------------------------------------------------------------------------
# Sphere: Gauss-Legendre in cos(theta_c) (ny points), trapezoid in phi_c (nx = 2 ny)
# calibrated: error ~1e-12 once ny >= a/2 + 20  (a = qmax R)
# --------------------------------------------------------------------------
def sphere_pas(R, petitk, grandk, base=40):
    """Return (pas, message). pas = ny ; nx = 2*pas."""
    if not auto_enabled():
        return base, ""
    a = size_parameter(petitk, grandk, R)
    pas = max(base, int(math.ceil(0.5 * a)) + 20)
    cap = grid_max(400)
    msg = ""
    if pas > cap:
        msg = ("WARNING: surface grid limited to %d (size parameter q*R = %.1f needs %d): "
               "results may be inaccurate. Raise PYSHS_GRID_MAX to increase it." % (cap, a, pas))
        pas = cap
    return pas, msg


# --------------------------------------------------------------------------
# Cylinder: periodic trapezoid around the axis (n_beta), analytic along the axis,
# trapezoid(phi) x Gauss-Legendre(cos theta) orientation grid (n_phi = n_theta = N).
# calibrated: n_beta >= q R + 14 ; N >= q sqrt(R^2 + L^2/4) + 12   (error <~ 1e-4)
# --------------------------------------------------------------------------
def cylinder_grid(R, length, petitk, grandk, base=(10, 10, 10)):
    """Return (n_beta, n_phi, n_theta, message)."""
    if not auto_enabled():
        return base[0], base[1], base[2], ""
    aR = size_parameter(petitk, grandk, R)
    aO = size_parameter(petitk, grandk, math.sqrt(R * R + 0.25 * length * length))
    n_beta = max(12, int(math.ceil(aR)) + 14, base[0])
    n_ori = max(12, int(math.ceil(aO)) + 12, base[1], base[2])
    msgs = []
    cap_b = grid_max(128)
    cap_o = grid_max(160)
    if n_beta > cap_b:
        msgs.append("WARNING: azimuthal grid limited to %d (needs %d for q*R = %.1f)." % (cap_b, n_beta, aR))
        n_beta = cap_b
    if n_ori > cap_o:
        msgs.append("WARNING: orientation grid limited to %d x %d (needs %d for size parameter %.1f)."
                    % (cap_o, cap_o, n_ori, aO))
        n_ori = cap_o
    if msgs:
        msgs.append("Results may be inaccurate; raise PYSHS_GRID_MAX to increase the limits.")
    return n_beta, n_ori, n_ori, "\n".join(msgs)


# --------------------------------------------------------------------------
# Dipole assemblies (SHS): trapezoid in phi and psi (n_phi = n_psi = N),
# Gauss-Legendre in cos(theta) (n_theta = N/2).  a = qmax * max|r| (distance to origin).
# calibrated on coherent shells and chains: N >= 2 a + 8   (error <~ 1e-5)
# --------------------------------------------------------------------------
def dipole_grid(rmax, petitk, grandk, base=(20, 10, 20)):
    """Return (n_phi, n_theta, n_psi, message)."""
    if not auto_enabled():
        return base[0], base[1], base[2], ""
    a = size_parameter(petitk, grandk, rmax)
    n = max(base[0], base[2], int(math.ceil(2.0 * a)) + 8)
    msg = ""
    cap = grid_max(80)
    if n > cap:
        msg = ("WARNING: orientation grid limited to %d x %d x %d (needs N = %d for size "
               "parameter q*rmax = %.1f): results may be inaccurate. Raise PYSHS_GRID_MAX "
               "to increase it (cost grows as N^3)." % (cap, cap // 2, cap, n, a))
        n = cap
    return n, max(base[1], n // 2), n, msg


# --------------------------------------------------------------------------
# Ellipsoid of revolution: surface = Gauss-Legendre in cos(alpha) x periodic trapezoid in
# beta; orientations = trapezoid(phi) x Gauss-Legendre(cos theta).  Half extent
# e = max(a, c) and size parameter aE = qmax * e.
# rules (same calibration as sphere / cylinder):
#   n_alpha >= aE/2 + 20, n_beta >= aE + 16, N_orientation >= aE + 12
# --------------------------------------------------------------------------
def ellipsoid_grid(a, c, petitk, grandk, base=(30, 20, 20, 20)):
    """Return (n_alpha, n_beta, n_phi, n_theta, message)."""
    if not auto_enabled():
        return base[0], base[1], base[2], base[3], ""
    aE = size_parameter(petitk, grandk, max(a, c))
    n_alpha = max(base[0], int(math.ceil(0.5 * aE)) + 20)
    n_beta = max(base[1], int(math.ceil(aE)) + 16)
    n_ori = max(base[2], base[3], int(math.ceil(aE)) + 12)
    caps = (grid_max(200), grid_max(300), grid_max(160))
    msgs = []
    if n_alpha > caps[0] or n_beta > caps[1]:
        msgs.append("WARNING: surface grid limited to %d x %d (needs %d x %d for size parameter %.1f)."
                    % (min(n_alpha, caps[0]), min(n_beta, caps[1]), n_alpha, n_beta, aE))
        n_alpha, n_beta = min(n_alpha, caps[0]), min(n_beta, caps[1])
    if n_ori > caps[2]:
        msgs.append("WARNING: orientation grid limited to %d x %d (needs %d)." % (caps[2], caps[2], n_ori))
        n_ori = caps[2]
    if msgs:
        msgs.append("Results may be inaccurate; raise PYSHS_GRID_MAX to increase the limits.")
    return n_alpha, n_beta, n_ori, n_ori, "\n".join(msgs)


# --------------------------------------------------------------------------
# Box (cube): faces are integrated analytically (sinc factors), so only the orientation
# grid (N x N x N : trapezoid phi, Gauss-Legendre cos theta, trapezoid psi) matters.
# Half diagonal d = sqrt(a^2+b^2+c^2)/2, size parameter aD = qmax * d.
# calibrated: N >= aD + 12
# --------------------------------------------------------------------------
def box_grid(a, b, c, petitk, grandk, base=(24, 24, 24)):
    """Return (n_phi, n_theta, n_psi, message)."""
    if not auto_enabled():
        return base[0], base[1], base[2], ""
    d = 0.5 * math.sqrt(a * a + b * b + c * c)
    aD = size_parameter(petitk, grandk, d)
    n = max(max(base), int(math.ceil(aD)) + 12)
    msg = ""
    cap = grid_max(80)
    if n > cap:
        msg = ("WARNING: orientation grid limited to %d x %d x %d (needs N = %d for size "
               "parameter q*d = %.1f): results may be inaccurate. Raise PYSHS_GRID_MAX to "
               "increase it (cost grows as N^3)." % (cap, cap, cap, n, aD))
        n = cap
    return n, n, n, msg

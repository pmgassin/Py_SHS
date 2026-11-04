"""Scattering-angle grid of the `angle_scattering` mode, common to all PySHS modules.

The user enters N = "number of points between 0 and 180 deg".  The grid has N intervals of
180/N degrees.

Because I_{V,H}(Theta) = I_{V,H}(-Theta) for the four intensities that are reported
(IH(0), IH(90), IV(0), IV(90)) -- checked numerically in every module, chiral tensors included --
a half range is sufficient.  Historically the modules did not agree: HRS, SHS and the cylinder
covered [-180, 0] (N+1 points) whereas the sphere, ellipsoid, cube and Mie covered [-180, +180]
(2N+1 points).

PYSHS_ANGLE_RANGE
    half      (default) Theta = 0 ... 180 deg, N+1 points     <- what the prompts announce
    negative  Theta = -180 ... 0 deg, N+1 points              (former HRS/SHS/cylinder output)
    full      Theta = -180 ... +180 deg, 2N+1 points          (former sphere/ellipsoid/cube/Mie output)
"""
import math
import os

import numpy as np

PI = math.pi


def _mode():
    m = os.environ.get("PYSHS_ANGLE_RANGE", "half").strip().lower()
    return m if m in ("half", "negative", "full") else "half"


def count(n):
    """Number of angles for N intervals."""
    return 2 * n + 1 if _mode() == "full" else n + 1


def scattering_angles(n):
    """Scattering angles Theta in radians (numpy array) for N intervals of 180/N deg."""
    k = np.arange(count(n))
    m = _mode()
    if m == "half":
        return k * (PI / n)
    return -PI + k * (PI / n)          # negative: -pi..0 ; full: -pi..+pi

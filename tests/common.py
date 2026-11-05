"""Shared helpers for the PySHS test-suite."""
import contextlib
import io
import os
import re
import subprocess
import sys
import tempfile

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.normpath(os.path.join(HERE, "..", "Src"))
if SRC not in sys.path:
    sys.path.insert(0, SRC)

PI = np.pi
NAMES = [a + b + c for a in "xyz" for b in "xyz" for c in "xyz"]

# Wavelength / indices used by the in-process physics tests
LAM, N_P, N_M = 800.0, 1.5, 1.33
K_FUND = 2 * PI * N_P / LAM          # k_omega (inside the particle)
K_SH = 4 * PI * N_M / LAM            # K_2omega (in the medium)
PREF2 = (K_SH * K_SH / (4 * PI)) ** 2


def quiet():
    """Context manager swallowing the progress messages printed by the modules."""
    return contextlib.redirect_stdout(io.StringIO())


def tensor_surface():
    """Achiral surface tensor, symmetric in (j,k), isotropic in the surface plane."""
    b = np.zeros((3, 3, 3))
    b[2, 2, 2] = 1.0
    b[2, 0, 0] = b[2, 1, 1] = 0.4
    b[0, 0, 2] = b[0, 2, 0] = b[1, 1, 2] = b[1, 2, 1] = 0.3
    return b


def tensor_chiral(seed=11):
    """Random chiral tensor, symmetric in (j,k), not Kleinman symmetric."""
    b = np.random.default_rng(seed).normal(size=(3, 3, 3))
    return 0.5 * (b + b.transpose(0, 2, 1))


def tensor_asym():
    """Isotropic in the plane but NOT symmetric in (j,k): only xxz and yyz are set."""
    b = np.zeros((3, 3, 3))
    b[2, 2, 2] = 1.0
    b[2, 0, 0] = b[2, 1, 1] = 0.5
    b[0, 0, 2] = b[1, 1, 2] = 1.0
    return b


def write_tensor(path, b):
    with open(path, "w") as f:
        f.write("# test tensor\n")
        for n, v in zip(NAMES, b.reshape(27)):
            f.write("%s %.12f\n" % (n, v))


ORIENT2 = "0 0 0 0 0 0\n0.3 0.5 0.2 4 0 2\n"           # two dipoles
ORIENT3 = "0.3 0.5 0.2 0 0 0\n1.1 0.7 2.0 40 10 -30\n2.5 1.9 0.4 -20 35 15\n"


def cli(module, keyword, stdin, files, env=None, workdir=None):
    """Run `python Src/<module>.py keyword files...` with `stdin`; return (stdout, workdir).

    Runs in `workdir` (a fresh temporary directory by default) so that `out_plot` and the
    other side files of the modules never collide between tests.
    """
    wd = workdir or tempfile.mkdtemp(prefix="pyshs_test_")
    e = dict(os.environ)
    for k in ("PYSHS_UNITS", "PYSHS_I_SCALE", "PYSHS_GRID", "PYSHS_GRID_SCALE", "PYSHS_GRID_MAX",
              "PYSHS_ANGLE_RANGE", "PYSHS_THETA_RULE", "PYSHS_LEGACY"):
        e.pop(k, None)
    e.update(env or {})
    r = subprocess.run([sys.executable, os.path.join(SRC, module + ".py"), keyword] + list(files),
                       input=stdin, capture_output=True, text=True, cwd=wd, env=e)
    if r.returncode != 0:
        raise RuntimeError("%s failed:\n%s\n%s" % (module, r.stdout[-500:], r.stderr[-800:]))
    return r.stdout, wd


def values(out, names=("av", "bv", "cv", "ah", "bh", "ch")):
    """Parse 'name = value' lines of a module's stdout."""
    res = []
    for n in names:
        m = re.search(r"^%s = (\S+)" % n, out, re.M)
        if m is None:
            raise AssertionError("value %r not found in output:\n%s" % (n, out[-600:]))
        res.append(float(m.group(1)))
    return np.array(res)


def read_table(path):
    """Numeric rows of an `out_plot` file (header skipped)."""
    rows = []
    with open(path) as f:
        for line in f.read().splitlines()[1:]:
            if line.strip():
                rows.append([float(x) for x in line.split()])
    return np.array(rows)


def rel(a, b, floor=1e-12):
    """max relative difference, ignoring entries that are ~0 in the reference."""
    a, b = np.asarray(a, float), np.asarray(b, float)
    if b.size == 0:
        return 0.0
    m = np.abs(b) > floor * max(np.abs(b).max(), 1e-300)
    return float(np.abs(a[m] / b[m] - 1).max()) if m.any() else 0.0


def fibonacci_sphere(n, radius):
    """n quasi-uniform points on a sphere: returns (phi, theta, positions)."""
    i = np.arange(n) + 0.5
    t = np.arccos(1 - 2 * i / n)
    p = PI * (1 + 5 ** 0.5) * i
    pos = radius * np.column_stack([np.sin(t) * np.cos(p), np.sin(t) * np.sin(p), np.cos(t)])
    return p % (2 * PI), t, pos


def write_text(path, text):
    with open(path, "w") as f:
        f.write(text)

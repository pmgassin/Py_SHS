"""Definition of the frozen regression cases (command-line runs, fixed inputs).

The snapshot in golden/golden.json was produced by *this* version of the Python code
(`python make_golden.py`).  It is NOT a reference from the former C programs: it only detects
unintended changes.  To validate against the C programs, see compare_outputs.py.
"""
import importlib.util
import os
import tempfile

import numpy as np

import common as C

HAVE_SCIPY = importlib.util.find_spec("scipy") is not None

# name: (module, keyword, stdin, tensor kind, needs orientation file, quantities)
CASES = {
    "HRS_polar_90": ("HRS", "polarplot_single", "800\n1.33\n0\n90\n", "chiral", False),
    "SHS_polar_90": ("SHS", "polarplot_single", "2\n800\n1.33\n0\n90\n", "chiral", True),
    "sphere_polar_90": ("sphere_SHS", "polarplot_single", "1\n800\n1.5\n0\n1.33\n0\n150\n90\n", "surface", False),
    "sphere_polar_150": ("sphere_SHS", "polarplot_single", "1\n800\n1.5\n0\n1.33\n0\n150\n150\n", "surface", False),
    "cylinder_polar_90": ("cylinder_SHS", "polarplot_single", "800\n1.5\n0\n1.33\n0\n50\n300\n90\n", "surface", False),
    "ellipsoid_polar_90": ("ellipsoid_SHS", "polarplot_single", "1\n800\n1.5\n0\n1.33\n0\n80\n400\n90\n", "surface", False),
    "cube_polar_90": ("cube_SHS", "polarplot_single", "1\n800\n1.5\n0\n1.33\n0\n200\n200\n200\n90\n", "surface", False),
    "HRS_angles_N2": ("HRS", "angle_scattering", "800\n1.33\n0\n2\n", "chiral", False),
    "sphere_integrate": ("sphere_SHS", "polarplot_integrate", "1\n800\n1.5\n0\n1.33\n0\n150\n60\n120\n3\n", "surface", False),
}
if HAVE_SCIPY:
    CASES["Mie_polar_90"] = ("shs_mie", "polarplot_single", "800\n1.5\n0\n1.33\n0\n100\n90\n", "surface", False)

QUANTITIES = ("av", "bv", "cv", "ah", "bh", "ch")


def run_case(name, workdir=None):
    mod, kw, stdin, tensor, needs_orient = CASES[name]
    d = workdir or tempfile.mkdtemp(prefix="pyshs_golden_")
    bf = os.path.join(d, "beta.txt")
    C.write_tensor(bf, C.tensor_chiral() if tensor == "chiral" else C.tensor_surface())
    of = os.path.join(d, "orient.txt")
    C.write_text(of, C.ORIENT2)
    files = (bf, of, "out.txt") if needs_orient else (bf, "out.txt")
    out, _ = C.cli(mod, kw, stdin, files, workdir=d)
    if kw == "angle_scattering" or mod == "shs_mie":
        rows = C.read_table(os.path.join(d, "out_plot"))
        if mod == "shs_mie":
            rows = rows[[0, 25, 50]]
        return [float(x) for x in rows.ravel()]
    return [float(x) for x in C.values(out, QUANTITIES)]


def compute_all():
    return {name: run_case(name) for name in CASES}

#!/usr/bin/env python3
"""Convert a Gaussian .log hyperpolarizability block to the PySHS beta input format.

    python script_gauss2PySHS.py job.log beta_output [--marker TEXT] [--occurrence K] [--offset N]

The script looks for the line whose first word is the marker (default ``Beta(-2w;w,w)``),
then reads 18 consecutive lines, ``offset`` lines after the marker (default 8), taking the
SECOND word of each line.  The 18 values are the components beta_{i,jk} with (jk) in the order
xx, xy, yy, xz, yz, zz and, for each, i = x, y, z.  They are expanded into the full 27-component
tensor with beta_ijk = beta_ikj (the intrinsic permutation symmetry of the last two indices
that PySHS assumes -- NOT Kleinman symmetry, which would also symmetrise the first index).

Notes
-----
* The layout of the block (marker, 8-line offset, 18 rows) is the one assumed by the original
  script; it could not be checked against a Gaussian version here.  The values are therefore
  validated (18 numeric rows, finite) and the block that was used is reported; if your Gaussian
  version prints it differently use --offset, or --marker for another block (for instance a
  static ``Beta(0;0,0)``, if your log has one in the same layout).
* If the log contains several blocks (several frequencies), the first one is used (the former
  behaviour) and the others are listed: choose with --occurrence K (1 = first).
"""
import argparse
import math
import sys

import numpy as np

LABELS = [
    "xxx", "xxy", "xxz", "xyx", "xyy", "xyz", "xzx", "xzy", "xzz",
    "yxx", "yxy", "yxz", "yyx", "yyy", "yyz", "yzx", "yzy", "yzz",
    "zxx", "zxy", "zxz", "zyx", "zyy", "zyz", "zzx", "zzy", "zzz",
]

# index of beta_{i,jk} in the 18 values : 3 * (position of jk in xx,xy,yy,xz,yz,zz) + i
_JK = {(0, 0): 0, (0, 1): 1, (1, 1): 2, (0, 2): 3, (1, 2): 4, (2, 2): 5}


def parse_gauss_float(s: str) -> float:
    """Parse Gaussian-style floats (0.123D-04, 1.23E+02, etc.)."""
    return float(s.strip().replace("D", "E").replace("d", "E"))


def expand(values):
    """18 contracted values -> full (3,3,3) tensor, symmetric in the last two indices."""
    beta = np.zeros((3, 3, 3))
    for i in range(3):
        for j in range(3):
            for k in range(3):
                beta[i, j, k] = values[3 * _JK[(min(j, k), max(j, k))] + i]
    return beta


def read_block(lines, marker, occurrence, offset):
    """Return (values[18], marker_line_index, all_marker_line_indices)."""
    found = [i for i, line in enumerate(lines) if line.split() and line.split()[0] == marker]
    if not found:
        raise ValueError("%s block not found in the Gaussian log." % marker)
    if not 1 <= occurrence <= len(found):
        raise ValueError("--occurrence %d requested but only %d block(s) found." % (occurrence, len(found)))
    index = found[occurrence - 1]
    start = index + offset
    values = np.zeros(18)
    for n in range(18):
        if start + n >= len(lines):
            raise ValueError("unexpected end of file while reading the beta block.")
        parts = lines[start + n].split()
        if len(parts) < 2:
            raise ValueError("cannot parse beta line %d: %r" % (start + n + 1, lines[start + n]))
        try:
            values[n] = parse_gauss_float(parts[1])
        except ValueError:
            raise ValueError("line %d: %r is not a number (is --offset right for this log?)"
                             % (start + n + 1, parts[1]))
    if not np.all(np.isfinite(values)):
        raise ValueError("non-finite value in the beta block.")
    return values, index, found


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("log")
    ap.add_argument("output")
    ap.add_argument("--marker", default="Beta(-2w;w,w)")
    ap.add_argument("--occurrence", type=int, default=1)
    ap.add_argument("--offset", type=int, default=8)
    a = ap.parse_args(argv)
    with open(a.log, "r", encoding="utf-8", errors="replace") as f:
        lines = f.readlines()
    try:
        values, index, found = read_block(lines, a.marker, a.occurrence, a.offset)
    except ValueError as exc:
        print("ERROR: %s" % exc, file=sys.stderr)
        return 1
    if len(found) > 1:
        print("WARNING: %d '%s' blocks found (lines %s); using number %d."
              % (len(found), a.marker, ", ".join(str(i + 1) for i in found), a.occurrence), file=sys.stderr)
    beta = expand(values)
    with open(a.output, "w", encoding="utf-8") as f:
        f.write("# beta in atomic units (a.u.) from Gaussian (%s, line %d of the log)\n" % (a.marker, index + 1))
        for idx, (i, j, k) in enumerate((i, j, k) for i in range(3) for j in range(3) for k in range(3)):
            f.write("%s %12.6e\n" % (LABELS[idx], beta[i, j, k]))
    print("Gaussian beta written to", a.output)
    return 0


if __name__ == "__main__":
    sys.exit(main())

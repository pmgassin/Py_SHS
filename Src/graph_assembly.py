#!/usr/bin/env python3
"""Plot an assembly of oriented dipoles (the orientation file read by SHS.py).

    python graph_assembly.py orientation_file output_image.png [--positions-first] [--arrow L]

Each line of the file is  phi theta psi x y z  (angles in radians, positions in nm), the
format documented in the manual.  With --positions-first the columns are read as
x y z phi theta psi instead.  (The former version tried to guess the order from the values
and got it wrong for valid files with an angle larger than 3.2 rad.)
The arrows show the dipole z axis.  Their length is 15 % of the size of the assembly unless
--arrow gives it (nm).
"""
import argparse
import sys

import numpy as np


def load(path, positions_first=False):
    rows = []
    with open(path, "r") as f:
        for line in f:
            parts = line.split()
            if len(parts) >= 6:
                try:
                    rows.append([float(p) for p in parts[:6]])
                except ValueError:
                    continue
    d = np.array(rows)
    if d.size == 0:
        return None
    if positions_first:
        d = d[:, [3, 4, 5, 0, 1, 2]]
    phi, theta, psi, x, y, z = d.T
    return phi, theta, psi, x, y, z


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("orientation_file")
    ap.add_argument("output_image")
    ap.add_argument("--positions-first", action="store_true")
    ap.add_argument("--arrow", type=float, default=None)
    a = ap.parse_args(argv)
    data = load(a.orientation_file, a.positions_first)
    if data is None:
        print("No dipole found in the orientation file.")
        return 1
    phi, theta, psi, x, y, z = data
    import matplotlib
    if not matplotlib.get_backend():
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    extent = max(np.ptp(x), np.ptp(y), np.ptp(z))
    bip = a.arrow if a.arrow else (0.15 * extent if extent > 0 else 5.0)
    fig = plt.figure(figsize=(8, 7))
    ax = fig.add_subplot(111, projection="3d")
    ax.quiver(x, y, z,
              bip * np.cos(phi) * np.sin(theta),
              bip * np.sin(phi) * np.sin(theta),
              bip * np.cos(theta),
              color="crimson", arrow_length_ratio=0.3)
    ax.set_xlabel("X (nm)")
    ax.set_ylabel("Y (nm)")
    ax.set_zlabel("Z (nm)")
    ax.set_title("Dipole assembly (%d dipoles)" % len(x))
    plt.tight_layout()
    plt.savefig(a.output_image, dpi=150)
    print("Assembly picture saved: %s" % a.output_image)
    return 0


if __name__ == "__main__":
    sys.exit(main())

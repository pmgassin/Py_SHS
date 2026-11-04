"""Dipoles on a sphere (no disorder), for SHS.py.

Usage:  python script_sphere.py <output_file>      (answers are read from stdin)

Output (one line per dipole):   phi theta psi x y z
The dipoles are quasi-uniformly distributed (Fibonacci spiral, random rotation of the spiral)
and oriented along the outward normal (psi = 0), the frame used by sphere_SHS.py.
A 3D picture <output_file>.png is written when matplotlib is available (the former version
overwrote the dipole file with the picture and crashed with matplotlib >= 3.6).
"""
import math
import random
import sys

import numpy as np


def fibonacci_sphere(samples=1, randomize=True):
    rnd = random.random() * samples if randomize else 1.0
    points = []
    offset = 2.0 / samples
    increment = math.pi * (3.0 - math.sqrt(5.0))
    for i in range(samples):
        y = ((i * offset) - 1) + (offset / 2)
        r = math.sqrt(1 - y * y)
        phi = ((i + rnd) % samples) * increment
        points.append([math.cos(phi) * r, y, math.sin(phi) * r])
    return points


def main():
    if len(sys.argv) < 2:
        print("usage: python script_sphere.py <output_file>")
        return 1
    uuv = sys.argv[1]
    a = float(input("radius of the sphere in nm? "))
    nombre = int(input("number of dipoles onto the sphere? "))
    uu = np.array(fibonacci_sphere(nombre))
    x, y, z = a * uu[:, 0], a * uu[:, 1], a * uu[:, 2]
    theta = np.arccos(np.clip(z / a, -1.0, 1.0))
    phi = np.arctan2(y, x)                       # azimuth, robust at the poles
    psi = np.zeros(nombre)
    with open(uuv, "w") as f:
        for i in range(nombre):
            f.write("%12.6f %12.6f %12.6f %12.6f %12.6f %12.6f\n" % (phi[i], theta[i], psi[i], x[i], y[i], z[i]))
    try:
        import matplotlib
        matplotlib.use("Agg", force=False)
        import matplotlib.pyplot as plt
        fig = plt.figure()
        ax = fig.add_subplot(111, projection="3d")
        ax.quiver(x, y, z, (a / 3) * np.cos(phi) * np.sin(theta), (a / 3) * np.sin(phi) * np.sin(theta),
                  (a / 3) * np.cos(theta), color="b")
        plt.savefig(uuv if uuv.lower().endswith((".png", ".pdf", ".jpg")) else uuv + ".png")
    except Exception as exc:                     # the picture is optional
        print("(no 3D picture: %s)" % exc)
    return 0


if __name__ == "__main__":
    sys.exit(main())

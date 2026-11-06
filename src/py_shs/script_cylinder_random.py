# Dipoles on the lateral surface of a cylinder (axis z), for SHS.py.
# Output lines: phi theta psi x y z ; dipole along the outward radial direction
# (phi = azimuth, theta = pi/2, psi = 0 -- the frame used by cylinder_SHS.py) + Gaussian noise
# sigma on phi and theta.  The cylinder has no end caps (as in cylinder_SHS.py).
import sys

import numpy as np
try:
    import matplotlib.pyplot as plt
except ImportError:                       # the picture is optional
    plt = None

uu = sys.argv[1]

a = float(input("cylinder radius in nm?"))
l = float(input("length of the cylinder in nm?"))
nombre1 = int(input("Number of dipoles around the circumference?"))
nombre2 = int(input("Number of dipoles on the length ?"))
ind = float(input("indice de desordre, sigma distribution ?"))
nombrefinal = nombre1 * nombre2
cc = float(l / nombre2) if nombre2 > 0 else 0.0

phi = np.zeros(nombrefinal)
theta = np.zeros(nombrefinal)
psi = np.zeros(nombrefinal)
x = np.zeros(nombrefinal)
y = np.zeros(nombrefinal)
z = np.zeros(nombrefinal)

kk = 0
for i in range(nombre1):
    for j in range(nombre2):
        phi[kk] = (np.random.randn(1) * ind + (i * 2 * np.pi / nombre1))[0]
        theta[kk] = (np.random.randn(1) * ind + np.pi / 2)[0]
        psi[kk] = 0
        x[kk] = a * np.cos(i * 2 * np.pi / nombre1)
        y[kk] = a * np.sin(i * 2 * np.pi / nombre1)
        z[kk] = -(l / 2) + cc * j
        kk += 1

b = 0.4 * a if a > 0 else 1.0          # arrow length of the picture

lines = []
for j in range(nombrefinal):
    lines.append("%12.6f %12.6f %12.6f %12.6f %12.6f %12.6f\n" % (phi[j], theta[j], psi[j], x[j], y[j], z[j]))

with open(uu, "w") as f:
    f.writelines(lines)

print("total number of dipoles on the cylinder:%d" % nombrefinal)
try:
    if plt is None:
        raise RuntimeError('matplotlib not installed')
    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")
    ax.quiver(x, y, z, b * np.cos(phi) * np.sin(theta), b * np.sin(phi) * np.sin(theta), b * np.cos(theta))
    plot_name = uu if uu.lower().endswith((".png", ".pdf", ".jpg")) else (uu + ".png")
    plt.savefig(plot_name)
except Exception as exc:                 # the picture is optional
    print("(no 3D picture: %s)" % exc)

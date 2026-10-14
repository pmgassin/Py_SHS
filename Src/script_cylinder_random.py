import matplotlib.pyplot as plt
import numpy as np
import sys

uu = sys.argv[1]

a = float(input("cylinder radius in nm?"))
l = float(input("length of the cylinder in nm?"))
nombre1 = int(input("Number of dipoles per radius?"))
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

xmax = np.amax(x) if nombrefinal else 0
ymax = np.amax(y) if nombrefinal else 0
zmax = np.amax(z) if nombrefinal else 0
b = (xmax + ymax + zmax) / 10 if (xmax + ymax + zmax) != 0 else 1.0

lines = []
for j in range(nombrefinal):
    lines.append("%12.6f %12.6f %12.6f %12.6f %12.6f %12.6f\n" % (phi[j], theta[j], psi[j], x[j], y[j], z[j]))

with open(uu, "w") as f:
    f.writelines(lines)

print("total number of dipoles on the cylinder:%d" % nombrefinal)
fig = plt.figure()
ax = fig.add_subplot(111, projection="3d")
ax.quiver(x, y, z, b * np.cos(phi) * np.sin(theta), b * np.sin(phi) * np.sin(theta), b * np.cos(theta))
plot_name = uu if uu.lower().endswith((".png", ".pdf", ".jpg")) else (uu + ".png")
plt.savefig(plot_name)

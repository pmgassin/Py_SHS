import math, random
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import csv
import sys


uuv=sys.argv[1]

a=float(input("radius of the sphere in nm?"))
nombre=int(input("number of dipoles onto the sphere?"))
ind=float(input("indice de desordre, sigma distribution ?"))
phi=np.zeros(nombre)
theta=np.zeros(nombre)
psi=np.zeros(nombre)
x=np.zeros(nombre)
y=np.zeros(nombre)
z=np.zeros(nombre)

def fibonacci_sphere(samples=1,randomize=True):
    rnd = 1.
    if randomize:
        rnd = random.random() * samples
    points = []
    offset = 2./samples
    increment = math.pi * (3. - math.sqrt(5.));

    for i in range(samples):
        y = ((i * offset) - 1) + (offset / 2);
        r = math.sqrt(1 - pow(y,2))
        phi = ((i + rnd) % samples) * increment
        x = math.cos(phi) * r
        z = math.sin(phi) * r
        points.append([x,y,z])

    return points
uu=fibonacci_sphere(nombre)

lines = []
for i in range(nombre):
    x[i]=a*uu[i][0]
    y[i]=a*uu[i][1]
    z[i]=a*uu[i][2]
    theta[i]=np.arccos(z[i]/a)
    psi[i]=0
    if y[i]>0:
        phi[i]=(np.arccos(x[i]/(a*np.sin(theta[i]))))
    else:
        phi[i]=-(np.arccos(x[i]/(a*np.sin(theta[i]))))
    theta[i]=theta[i]+(np.random.randn(1)*ind)[0]
    phi[i]=phi[i]+(np.random.randn(1)*ind)[0]
    lines.append("%12.6f %12.6f %12.6f %12.6f %12.6f %12.6f\n" % (phi[i],theta[i],psi[i],x[i],y[i],z[i]))

with open(uuv, "w") as f:
    f.writelines(lines)

# Sauvegarde du graphique sous un nom distinct (évite d'écraser le .txt)
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.quiver(x,y,z,(a/3)*np.cos(phi)*np.sin(theta),(a/3)*np.sin(phi)*np.sin(theta),(a/3)*np.cos(theta),color='b')
plot_name = uuv if uuv.lower().endswith(('.png', '.pdf', '.jpg')) else (uuv + ".png")
plt.savefig(plot_name)
# plt.show()  # désactivé pour usage non-interactif (GUI)



#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import sys

if len(sys.argv) < 3:
    print("Usage: python graph_assembly.py <orientation_file> <output_image.png>")
    sys.exit(1)

orient_file = sys.argv[1]
output_image = sys.argv[2]

x, y, z, phi, theta, psi = [], [], [], [], [], []

with open(orient_file, 'r') as f:
    for line in f:
        parts = line.split()
        if len(parts) >= 6:
            try:
                # Format du manuel : φ' θ' ψ' x' y' z'
                # ou parfois x y z φ θ ψ → on détecte
                vals = [float(p) for p in parts[:6]]
                # Heuristique simple : si les 3 premiers sont entre -π et π → angles d'abord
                if abs(vals[0]) <= 3.2 and abs(vals[1]) <= 3.2:
                    phi.append(vals[0])
                    theta.append(vals[1])
                    psi.append(vals[2])
                    x.append(vals[3])
                    y.append(vals[4])
                    z.append(vals[5])
                else:
                    x.append(vals[0])
                    y.append(vals[1])
                    z.append(vals[2])
                    phi.append(vals[3])
                    theta.append(vals[4])
                    psi.append(vals[5])
            except ValueError:
                continue

if not x:
    print("Aucune donnée d'orientation trouvée")
    sys.exit(1)

x, y, z = np.array(x), np.array(y), np.array(z)
phi, theta = np.array(phi), np.array(theta)

bip = 5.0
fig = plt.figure(figsize=(8, 7))
ax = fig.add_subplot(111, projection='3d')
ax.quiver(x, y, z,
          bip * np.cos(phi) * np.sin(theta),
          bip * np.sin(phi) * np.sin(theta),
          bip * np.cos(theta),
          color='crimson', arrow_length_ratio=0.3)
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.set_title('Assemblage des dipôles')
plt.tight_layout()
plt.savefig(output_image, dpi=150)
print(f"Graphique assemblage sauvegardé : {output_image}")

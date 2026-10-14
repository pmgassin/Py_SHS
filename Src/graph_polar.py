#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import sys
from pathlib import Path

if len(sys.argv) < 3:
    print("Usage: python graph_polar.py <out_plot_file> <output_image.png>")
    sys.exit(1)

input_file = sys.argv[1]
output_image = sys.argv[2]

gamma, IV, IH = [], [], []

with open(input_file, 'r') as f:
    lines = f.readlines()
    # On ignore la première ligne si c'est un en-tête
    start = 1 if lines and not lines[0].strip()[0].isdigit() else 0
    for line in lines[start:]:
        parts = line.split()
        if len(parts) >= 3:
            try:
                gamma.append(float(parts[0]))
                IV.append(float(parts[1]))
                IH.append(float(parts[2]))
            except ValueError:
                continue

if not gamma:
    print("Aucune donnée trouvée dans", input_file)
    sys.exit(1)

plt.figure(figsize=(8, 5))
plt.plot(gamma, IV, label='IV', linewidth=2)
plt.plot(gamma, IH, label='IH', linewidth=2)
plt.legend()
plt.xlabel('Gamma (°)')
plt.ylabel('SHS Intensity')
plt.title('Polarization-resolved SHS')
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(output_image, dpi=150)
print(f"Graphique polar sauvegardé : {output_image}")

#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import sys

if len(sys.argv) < 3:
    print("Usage: python graph_angle.py <out_plot_file> <output_image.png>")
    sys.exit(1)

input_file = sys.argv[1]
output_image = sys.argv[2]

gamma, IH_0, IH_90, IV_0, IV_90 = [], [], [], [], []

with open(input_file, 'r') as f:
    lines = f.readlines()
    start = 1 if lines and not lines[0].strip()[0].isdigit() else 0
    for line in lines[start:]:
        parts = line.split()
        if len(parts) >= 5:
            try:
                gamma.append(float(parts[0]))
                IH_0.append(float(parts[1]))
                IH_90.append(float(parts[2]))
                IV_0.append(float(parts[3]))
                IV_90.append(float(parts[4]))
            except ValueError:
                continue

if not gamma:
    print("Aucune donnée trouvée")
    sys.exit(1)

# Symétrisation (comme dans le script original)
n = len(gamma)
gamma_full = gamma + [-g for g in reversed(gamma[:-1])]
IH_0_full = IH_0 + list(reversed(IH_0[:-1]))
IH_90_full = IH_90 + list(reversed(IH_90[:-1]))
IV_0_full = IV_0 + list(reversed(IV_0[:-1]))
IV_90_full = IV_90 + list(reversed(IV_90[:-1]))

plt.figure(figsize=(9, 5))
plt.plot(gamma_full, IV_0_full, label='IV(0°)')
plt.plot(gamma_full, IV_90_full, label='IV(90°)')
plt.plot(gamma_full, IH_0_full, label='IH(0°)')
plt.plot(gamma_full, IH_90_full, label='IH(90°)')
plt.legend()
plt.xlabel('Angle (°)')
plt.ylabel('SHS Intensity')
plt.title('Angle-resolved SHS')
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(output_image, dpi=150)
print(f"Graphique angle sauvegardé : {output_image}")

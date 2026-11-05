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
    for line in f:
        parts = line.split()
        if len(parts) >= 5:
            try:
                row = [float(p) for p in parts[:5]]
            except ValueError:          # header line
                continue
            gamma.append(row[0])
            IH_0.append(row[1]); IH_90.append(row[2])
            IV_0.append(row[3]); IV_90.append(row[4])

if not gamma:
    print("Aucune donnée trouvée")
    sys.exit(1)

# The four plotted intensities satisfy I(Theta) = I(-Theta).  The modules now write
# Theta in [0, 180] by default (older files: [-180, 0] or [-180, +180]).  When the data
# cover only one sign, mirror them to obtain the full -180 ... +180 curve.
order = np.argsort(gamma)
g = np.asarray(gamma)[order]
curves = [np.asarray(c)[order] for c in (IH_0, IH_90, IV_0, IV_90)]
if g.min() >= -1e-12 or g.max() <= 1e-12:
    keep = np.abs(g) > 1e-12                    # do not duplicate Theta = 0
    gamma_full = np.concatenate([-g[keep], g])
    srt = np.argsort(gamma_full)
    gamma_full = gamma_full[srt]
    IH_0_full, IH_90_full, IV_0_full, IV_90_full = [
        np.concatenate([c[keep], c])[srt] for c in curves]
else:
    gamma_full = g
    IH_0_full, IH_90_full, IV_0_full, IV_90_full = curves

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

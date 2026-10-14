#!/usr/bin/env python3
"""Convert Gaussian .log hyperpolarizability block to PySHS beta input format."""
import numpy as np
import sys
import re

uv = sys.argv[1]
ut = sys.argv[2]

LABELS = [
    "xxx", "xxy", "xxz", "xyx", "xyy", "xyz", "xzx", "xzy", "xzz",
    "yxx", "yxy", "yxz", "yyx", "yyy", "yyz", "yzx", "yzy", "yzz",
    "zxx", "zxy", "zxz", "zyx", "zyy", "zyz", "zzx", "zzy", "zzz",
]

def parse_gauss_float(s: str) -> float:
    """Parse Gaussian-style floats (0.123D-04, 1.23E+02, etc.)."""
    s = s.strip().replace("D", "E").replace("d", "E")
    return float(s)

with open(uv, "r", encoding="utf-8", errors="replace") as f:
    lines = f.readlines()

index = None
for i, line in enumerate(lines):
    if line.split() and line.split()[0] == "Beta(-2w;w,w)":
        index = i
        break

if index is None:
    print("ERROR: Beta(-2w;w,w) block not found in Gaussian log.", file=sys.stderr)
    sys.exit(1)

# Gaussian contracted beta table starts 8 lines after the marker
kk = index + 8
beta_contract_valeur = np.zeros(18)
for i in range(18):
    if kk + i >= len(lines):
        print("ERROR: unexpected end of file while reading beta block.", file=sys.stderr)
        sys.exit(1)
    parts = lines[kk + i].split()
    if len(parts) < 2:
        print(f"ERROR: cannot parse beta line: {lines[kk + i]!r}", file=sys.stderr)
        sys.exit(1)
    beta_contract_valeur[i] = parse_gauss_float(parts[1])

beta = np.zeros((3, 3, 3))
# Map contracted components to full tensor (Kleinman-type filling used by original script)
beta[0, 0, 0] = beta_contract_valeur[0]
beta[0, 0, 1] = beta_contract_valeur[3]
beta[0, 0, 2] = beta_contract_valeur[9]
beta[0, 1, 0] = beta_contract_valeur[3]
beta[0, 1, 1] = beta_contract_valeur[6]
beta[0, 1, 2] = beta_contract_valeur[12]
beta[0, 2, 0] = beta_contract_valeur[9]
beta[0, 2, 1] = beta_contract_valeur[12]
beta[0, 2, 2] = beta_contract_valeur[15]

beta[1, 0, 0] = beta_contract_valeur[1]
beta[1, 0, 1] = beta_contract_valeur[4]
beta[1, 0, 2] = beta_contract_valeur[10]
beta[1, 1, 0] = beta_contract_valeur[4]
beta[1, 1, 1] = beta_contract_valeur[7]
beta[1, 1, 2] = beta_contract_valeur[13]
beta[1, 2, 0] = beta_contract_valeur[10]
beta[1, 2, 1] = beta_contract_valeur[13]
beta[1, 2, 2] = beta_contract_valeur[16]

beta[2, 0, 0] = beta_contract_valeur[2]
beta[2, 0, 1] = beta_contract_valeur[5]
beta[2, 0, 2] = beta_contract_valeur[11]
beta[2, 1, 0] = beta_contract_valeur[5]
beta[2, 1, 1] = beta_contract_valeur[8]
beta[2, 1, 2] = beta_contract_valeur[14]
beta[2, 2, 0] = beta_contract_valeur[11]
beta[2, 2, 1] = beta_contract_valeur[14]
beta[2, 2, 2] = beta_contract_valeur[17]

with open(ut, "w", encoding="utf-8") as f:
    f.write("# beta in atomic units (a.u.) from Gaussian\n")
    idx = 0
    for i in range(3):
        for j in range(3):
            for k in range(3):
                f.write("%s %12.6e\n" % (LABELS[idx], beta[i, j, k]))
                idx += 1

print("Gaussian beta written to", ut)

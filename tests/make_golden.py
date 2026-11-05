"""Regenerate golden/golden.json.  Run only after an intentional change of the numerics:

    python tests/make_golden.py
"""
import json
import os

import golden_cases as G

path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "golden", "golden.json")
data = G.compute_all()
os.makedirs(os.path.dirname(path), exist_ok=True)
with open(path, "w") as f:
    json.dump(data, f, indent=1, sort_keys=True)
print("wrote %s (%d cases)" % (path, len(data)))

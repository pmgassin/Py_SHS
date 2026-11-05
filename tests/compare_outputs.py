#!/usr/bin/env python3
"""Compare two PySHS output files (e.g. one produced by the former C program and one by this
Python version) quantity by quantity.

    python tests/compare_outputs.py reference.txt new.txt [--rtol 1e-6] [--all-keys]

Compared:
  * the summary quantities  av bv bv1 bv2 cv ah bh bh1 bh2 ch I2v I4v I2h I4h  (lines such as
    "av = 1.007732" or "av=1.007732"), in order of appearance; with --all-keys every
    "name = number" pair is compared;
  * numeric tables (out_plot style: rows of >= 3 numbers), row by row.
Exit status 0 if every difference is within tolerance, 1 otherwise.

Tolerance: |a-b| <= rtol * max(|a|,|b|,  scale*1e-12)  where scale is the largest magnitude
of the quantities of the same kind (so that values that are zero up to rounding noise do not
raise false alarms).
"""
import argparse
import re
import sys

KEYS = ("av", "bv", "bv1", "bv2", "cv", "ah", "bh", "bh1", "bh2", "ch", "I2v", "I4v", "I2h", "I4h")
PAIR = re.compile(r"\b([A-Za-z][A-Za-z0-9_]*)\s*=\s*([-+]?(?:\d+\.?\d*|\.\d+)(?:[eEdD][-+]?\d+)?)")
NUM = re.compile(r"^[-+]?(?:\d+\.?\d*|\.\d+)(?:[eEdD][-+]?\d+)?$")


def _f(s):
    return float(s.replace("D", "E").replace("d", "e"))


def parse(path, all_keys=False):
    pairs, rows = [], []
    with open(path, encoding="utf-8", errors="replace") as f:
        for line in f:
            for k, v in PAIR.findall(line):
                if all_keys or k in KEYS:
                    pairs.append((k, _f(v)))
            toks = line.split()
            if len(toks) >= 3 and all(NUM.match(t) for t in toks):
                rows.append([_f(t) for t in toks])
    return pairs, rows


def close(a, b, rtol, scale):
    return abs(a - b) <= rtol * max(abs(a), abs(b), scale * 1e-12)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("reference"); ap.add_argument("new")
    ap.add_argument("--rtol", type=float, default=1e-6)
    ap.add_argument("--all-keys", action="store_true")
    a = ap.parse_args(argv)
    p1, r1 = parse(a.reference, a.all_keys)
    p2, r2 = parse(a.new, a.all_keys)
    bad = 0
    if len(p1) != len(p2):
        print("different number of summary quantities: %d vs %d" % (len(p1), len(p2))); bad += 1
    scale = max([abs(v) for _, v in p1 + p2] + [0.0])
    for (k1, v1), (k2, v2) in zip(p1, p2):
        if k1 != k2 or not close(v1, v2, a.rtol, scale):
            print("DIFF %-5s reference=%.8e new=%.8e (%s)" % (k1, v1, v2, k2)); bad += 1
    if len(r1) != len(r2):
        print("different number of table rows: %d vs %d" % (len(r1), len(r2))); bad += 1
    tscale = max([abs(x) for r in r1 + r2 for x in r[1:]] + [0.0])
    for i, (u, v) in enumerate(zip(r1, r2)):
        if len(u) != len(v):
            print("row %d: different number of columns" % i); bad += 1; continue
        for j, (x, y) in enumerate(zip(u, v)):
            tol_scale = abs(x) if j == 0 else tscale
            if not close(x, y, a.rtol, tol_scale):
                print("DIFF row %d column %d: reference=%.8e new=%.8e" % (i, j, x, y)); bad += 1
    print("%d summary quantities, %d table rows compared: %s" %
          (min(len(p1), len(p2)), min(len(r1), len(r2)), "OK" if not bad else "%d difference(s)" % bad))
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())

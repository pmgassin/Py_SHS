"""Console entry points (``pyshs-hrs``, ``pyshs-shs`` ...).

Each command keeps the historical protocol of the scripts:
    pyshs-<module> <keyword> <beta_file> [<orientation_file>] <output_file>
with the numerical parameters read from stdin.
"""
import importlib
import sys


def _run(module_name):
    mod = importlib.import_module(f"pyshs.{module_name}")
    sys.exit(mod.main(sys.argv))


def hrs():
    _run("HRS")


def shs():
    _run("SHS")


def sphere():
    _run("sphere_SHS")


def cylinder():
    _run("cylinder_SHS")


def ellipsoid():
    _run("ellipsoid_SHS")


def cube():
    _run("cube_SHS")


def mie():
    """pyshs-mie <keyword> <beta_file> <output_file>  (parameters on stdin)."""
    import os
    from pyshs import shs_mie
    if len(sys.argv) < 4:
        print("Usage: pyshs-mie <keyword> <beta_file> <output_file>   (answers on stdin)")
        sys.exit(1)
    kw, beta_f, out_f = sys.argv[1:4]
    answers = [line.strip() for line in sys.stdin if line.strip()]
    work = os.path.dirname(os.path.abspath(out_f)) or "."
    shs_mie.run_mie_cli(kw, beta_f, out_f, answers, work_dir=work)

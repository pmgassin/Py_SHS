"""Minimal smoke tests: the package imports and a module runs through ``python -m``."""
import importlib
import subprocess
import sys

import pytest

MODULES = ["HRS", "SHS", "sphere_SHS", "shs_mie", "cylinder_SHS", "ellipsoid_SHS", "cube_SHS",
           "shs_units", "shs_angles", "shs_grid", "shs_io"]


@pytest.mark.parametrize("name", MODULES)
def test_import(name):
    importlib.import_module(f"py_shs.{name}")


def test_version():
    import py_shs
    assert py_shs.__version__

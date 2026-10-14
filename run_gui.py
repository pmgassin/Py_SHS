#!/usr/bin/env python3
"""Launch the PySHS graphical interface."""
import runpy
from pathlib import Path
runpy.run_path(str(Path(__file__).resolve().parent / "Src" / "pyshs_gui.py"), run_name="__main__")

#!/usr/bin/env python3
"""
PySHS V4.0 – Tkinter GUI (100 % Python / NumPy)
Tabs: HRS | SHS | Sphere | Cylinder | Ellipsoid | Cube
"""

import tkinter as tk
from tkinter import ttk, filedialog, messagebox, scrolledtext
import subprocess
import os
import sys
import shutil
import threading
import tempfile
from pathlib import Path
import matplotlib
matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
import numpy as np

# Package layout (installed with pip):
#   py_shs/               pure-Python modules + GUI + generation scripts  (SRC_DIR)
#   ~/PySHS_Work/        input/output files, temp jobs, out_plot          (WORK_DIR)
# The working directory is NOT inside the package (site-packages may be read-only).
# It can be changed with the environment variable PYSHS_WORK.
SRC_DIR = Path(__file__).resolve().parent
BASE_DIR = SRC_DIR
WORK_DIR = Path(os.environ.get("PYSHS_WORK", Path.home() / "PySHS_Work")).expanduser()
WORK_DIR.mkdir(parents=True, exist_ok=True)

# Pure-Python modules (no C compilation)
SRC_HRS = SRC_DIR / "HRS.py"
SRC_SHS = SRC_DIR / "SHS.py"
SRC_SPHERE = SRC_DIR / "sphere_SHS.py"
SRC_MIE = SRC_DIR / "shs_mie.py"
SRC_CYLINDER = SRC_DIR / "cylinder_SHS.py"
SRC_ELLIPSOID = SRC_DIR / "ellipsoid_SHS.py"
SRC_CUBE = SRC_DIR / "cube_SHS.py"

# Known generation scripts: display name → (file, list of (label, key, default))
WORK_KEEP_DAYS = 7      # run directories (Work/pyshs_*) older than this are removed


def _cleanup_old_runs(max_age_days=WORK_KEEP_DAYS):
    """Remove the temporary run directories created by the GUI (Work/pyshs_*) that are older
    than `max_age_days`; they used to accumulate for ever.  Only directories created by the
    GUI itself (prefix 'pyshs_') are touched."""
    import time
    limit = time.time() - max_age_days * 86400.0
    try:
        for d in WORK_DIR.iterdir():
            if d.is_dir() and d.name.startswith("pyshs_"):
                try:
                    if d.stat().st_mtime < limit:
                        shutil.rmtree(d, ignore_errors=True)
                except OSError:
                    pass
    except OSError:
        pass


UNITS_CHOICES = ("module default", "rad : r\u00b2|E|\u00b2 (with radiative prefactor)",
                 "raw : coherent sum (no prefactor)")
UNITS_ENV = {UNITS_CHOICES[0]: "", UNITS_CHOICES[1]: "rad", UNITS_CHOICES[2]: "raw"}

GEN_SCRIPTS = {
    "Sphere (script_sphere_random.py)": (
        "script_sphere_random.py",
        [
            ("Radius (nm)", "radius", "50"),
            ("Number of dipoles", "n_dipoles", "100"),
            ("σ disorder (rad)", "sigma", "0.3"),
        ],
    ),
    "Spheroid (script_spheroid_random.py)": (
        "script_spheroid_random.py",
        [
            ("Radius a (nm)", "radius_a", "50"),
            ("Radius c (nm)", "radius_c", "30"),
            ("Number of dipoles", "n_dipoles", "100"),
            ("σ disorder (rad)", "sigma", "0.3"),
        ],
    ),
    "Cylinder (script_cylinder_random.py)": (
        "script_cylinder_random.py",
        [
            ("Radius (nm)", "radius", "20"),
            ("Length (nm)", "length", "80"),
            ("Dipoles around circumference", "n_radius", "12"),
            ("Dipoles along length", "n_length", "10"),
            ("σ disorder (rad)", "sigma", "0.3"),
        ],
    ),
    "Custom script…": (None, []),
}

EXAMPLE_BETA = """# Example hyperpolarizability tensor (a.u.)
xxx 0.0
xxy 0.0
xxz 0.0
xyx 0.0
xyy 0.0
xyz 0.0
xzx 0.0
xzy 0.0
xzz 0.0
yxx 0.0
yxy 0.0
yxz 0.0
yyx 0.0
yyy 0.0
yyz 0.0
yzx 0.0
yzy 0.0
yzz 0.0
zxx 0.0
zxy 0.0
zxz 0.0
zyx 0.0
zyy 0.0
zyz 0.0
zzx 0.0
zzy 0.0
zzz 1.0
"""


def check_python_module(src: Path, log_callback=None) -> bool:
    """Verify that a pure-Python calculation module is present."""
    if not src.exists():
        if log_callback:
            log_callback(f"ERROR: Python module not found: {src}")
        return False
    if log_callback:
        log_callback(f"OK: {src.name}")
    return True


class PySHSApp(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("PySHS V4.0 – Graphical Interface (100 % Python)")
        self.geometry("1280x900")
        self.minsize(1000, 750)

        self.job_running = False
        self._current_proc = None
        self._last_intensity_fig = None
        self._last_assembly_fig = None

        # HRS
        self.hrs_keyword = tk.StringVar(value="polarplot_single")
        self.hrs_wavelength = tk.StringVar(value="800")
        self.hrs_n_real = tk.StringVar(value="1.33")
        self.hrs_n_imag = tk.StringVar(value="0.0")
        self.hrs_theta = tk.StringVar(value="90")
        self.hrs_theta_start = tk.StringVar(value="0")
        self.hrs_theta_end = tk.StringVar(value="180")
        self.hrs_n_points = tk.StringVar(value="20")

        # SHS
        self.shs_keyword = tk.StringVar(value="polarplot_single")
        self.shs_n_dipoles = tk.StringVar(value="100")
        self.shs_wavelength = tk.StringVar(value="800")
        self.shs_n_real = tk.StringVar(value="1.33")
        self.shs_n_imag = tk.StringVar(value="0.0")
        self.shs_theta = tk.StringVar(value="90")
        self.shs_theta_start = tk.StringVar(value="0")
        self.shs_theta_end = tk.StringVar(value="180")
        self.shs_n_points = tk.StringVar(value="20")

        # Sphere colloidal (NLWKB: n_p particle, n_m medium)
        self.sph_keyword = tk.StringVar(value="polarplot_single")
        self.sph_wavelength = tk.StringVar(value="800")
        self.sph_np_real = tk.StringVar(value="1.5")
        self.sph_np_imag = tk.StringVar(value="0.0")
        self.sph_nm_real = tk.StringVar(value="1.33")
        self.sph_nm_imag = tk.StringVar(value="0.0")
        self.sph_radius = tk.StringVar(value="50")
        # Intensity units / scale: one setting shared by ALL calculation modules
        # (the "Intensity scale factor" row of the Sphere panel edits the same variable)
        self.out_units = tk.StringVar(value=UNITS_CHOICES[0])
        self.out_i_scale = tk.StringVar(value="1")
        self.sph_i_scale = self.out_i_scale
        self.sph_theta = tk.StringVar(value="90")
        self.sph_theta_start = tk.StringVar(value="0")
        self.sph_theta_end = tk.StringVar(value="180")
        self.sph_n_points = tk.StringVar(value="20")

        # Cylinder colloidal (NLWKB: n_p particle, n_m medium)
        self.cyl_keyword = tk.StringVar(value="polarplot_single")
        self.cyl_wavelength = tk.StringVar(value="800")
        self.cyl_np_real = tk.StringVar(value="1.5")
        self.cyl_np_imag = tk.StringVar(value="0.0")
        self.cyl_nm_real = tk.StringVar(value="1.33")
        self.cyl_nm_imag = tk.StringVar(value="0.0")
        self.cyl_radius = tk.StringVar(value="20")
        self.cyl_length = tk.StringVar(value="80")
        self.cyl_theta = tk.StringVar(value="90")
        self.cyl_theta_start = tk.StringVar(value="0")
        self.cyl_theta_end = tk.StringVar(value="180")
        self.cyl_n_points = tk.StringVar(value="20")

        # Ellipsoid of revolution (NLWKB)
        self.ell_keyword = tk.StringVar(value="polarplot_single")
        self.ell_wavelength = tk.StringVar(value="800")
        self.ell_np_real = tk.StringVar(value="1.5")
        self.ell_np_imag = tk.StringVar(value="0.0")
        self.ell_nm_real = tk.StringVar(value="1.33")
        self.ell_nm_imag = tk.StringVar(value="0.0")
        self.ell_radius = tk.StringVar(value="20")   # equatorial radius a
        self.ell_length = tk.StringVar(value="80")   # length along axis = 2c
        self.ell_theta = tk.StringVar(value="90")
        self.ell_theta_start = tk.StringVar(value="0")
        self.ell_theta_end = tk.StringVar(value="180")
        self.ell_n_points = tk.StringVar(value="20")

        # Cube / rectangular box (NLWKB)
        self.cube_keyword = tk.StringVar(value="polarplot_single")
        self.cube_wavelength = tk.StringVar(value="800")
        self.cube_np_real = tk.StringVar(value="1.5")
        self.cube_np_imag = tk.StringVar(value="0.0")
        self.cube_nm_real = tk.StringVar(value="1.33")
        self.cube_nm_imag = tk.StringVar(value="0.0")
        self.cube_a = tk.StringVar(value="20")   # edge along x
        self.cube_b = tk.StringVar(value="20")   # edge along y
        self.cube_c = tk.StringVar(value="20")   # edge along z
        self.cube_theta = tk.StringVar(value="90")
        self.cube_theta_start = tk.StringVar(value="0")
        self.cube_theta_end = tk.StringVar(value="180")
        self.cube_n_points = tk.StringVar(value="20")

        # Generation
        self.gen_script_choice = tk.StringVar(value=list(GEN_SCRIPTS.keys())[0])
        self.gen_vars = {}  # key → StringVar
        self.custom_script_path = tk.StringVar(value="")
        self.custom_stdin = tk.StringVar(value="")

        self._build_ui()
        self._update_hrs_params()
        self._update_shs_params()
        self._update_sph_params()
        self._update_cyl_params()
        self._update_ell_params()
        self._update_cube_params()
        self._on_gen_script_changed()
        self.after(200, self._initial_compile)

    # ==================================================================
    def _build_ui(self):
        header = tk.Frame(self, bg="#1e40af", height=52)
        header.pack(fill=tk.X)
        header.pack_propagate(False)
        tk.Label(header, text="PySHS V4.0 – Second Harmonic Scattering (pure Python)",
                 font=("Segoe UI", 15, "bold"), fg="white", bg="#1e40af").pack(side=tk.LEFT, padx=16, pady=10)

        lf_units = ttk.LabelFrame(self, text="Intensity units (all modules)", padding=(8, 2))
        lf_units.pack(fill=tk.X, padx=6, pady=(6, 0))
        ttk.Label(lf_units, text="Units:").pack(side=tk.LEFT)
        ttk.Combobox(lf_units, textvariable=self.out_units, values=UNITS_CHOICES,
                     state="readonly", width=34).pack(side=tk.LEFT, padx=(4, 12))
        ttk.Label(lf_units, text="Intensity scale factor:").pack(side=tk.LEFT)
        ttk.Entry(lf_units, textvariable=self.out_i_scale, width=10).pack(side=tk.LEFT, padx=(4, 12))
        ttk.Label(lf_units, foreground="#64748b",
                  text="rad = r\u00b2|E|\u00b2 with (|K|\u00b2/4\u03c0)\u00b2 ; raw = coherent sum, no prefactor. "
                       "Module default: HRS/SHS raw, others rad.").pack(side=tk.LEFT)

        paned = tk.PanedWindow(self, orient=tk.VERTICAL, sashrelief=tk.RAISED, sashwidth=5)
        paned.pack(fill=tk.BOTH, expand=True, padx=6, pady=6)

        top = ttk.Frame(paned)
        # Larger top pane (HRS/SHS tabs) so all SHS controls are visible
        paned.add(top, height=480)

        self.main_nb = ttk.Notebook(top)
        self.main_nb.pack(fill=tk.BOTH, expand=True)
        self.tab_hrs = ttk.Frame(self.main_nb, padding=6)
        self.tab_shs = ttk.Frame(self.main_nb, padding=6)
        self.tab_sphere = ttk.Frame(self.main_nb, padding=6)
        self.tab_cylinder = ttk.Frame(self.main_nb, padding=6)
        self.tab_ellipsoid = ttk.Frame(self.main_nb, padding=6)
        self.tab_cube = ttk.Frame(self.main_nb, padding=6)
        self.main_nb.add(self.tab_hrs, text="  HRS Calculation  ")
        self.main_nb.add(self.tab_shs, text="  SHS Calculation  ")
        self.main_nb.add(self.tab_sphere, text="  Sphere Object  ")
        self.main_nb.add(self.tab_cylinder, text="  Cylinder Object  ")
        self.main_nb.add(self.tab_ellipsoid, text="  Ellipsoid Object  ")
        self.main_nb.add(self.tab_cube, text="  Cube Object  ")
        self._build_hrs_tab(self.tab_hrs)
        self._build_shs_tab(self.tab_shs)
        self._build_sphere_tab(self.tab_sphere)
        self._build_cylinder_tab(self.tab_cylinder)
        self._build_ellipsoid_tab(self.tab_ellipsoid)
        self._build_cube_tab(self.tab_cube)

        bottom = ttk.Frame(paned)
        paned.add(bottom)
        self._build_results(bottom)

    # ------------------------------------------------------------------
    # HRS
    # ------------------------------------------------------------------
    def _build_hrs_tab(self, parent):
        left = ttk.Frame(parent)
        left.pack(side=tk.LEFT, fill=tk.Y, padx=(0, 8))

        lf_kw = ttk.LabelFrame(left, text="Calculation type", padding=6)
        lf_kw.pack(fill=tk.X, pady=(0, 6))
        cb = ttk.Combobox(lf_kw, textvariable=self.hrs_keyword,
                          values=["polarplot_single", "polarplot_integrate", "angle_scattering"],
                          state="readonly", width=26)
        cb.pack(fill=tk.X)
        cb.bind("<<ComboboxSelected>>", lambda e: self._update_hrs_params())

        self.hrs_params_frame = ttk.LabelFrame(left, text="Parameters", padding=6)
        self.hrs_params_frame.pack(fill=tk.X, pady=(0, 6))
        self.hrs_param_widgets = {}
        for i, (key, label, var) in enumerate([
            ("wavelength", "Wavelength (nm)", self.hrs_wavelength),
            ("n_real", "n – real part", self.hrs_n_real),
            ("n_imag", "n – imaginary part", self.hrs_n_imag),
            ("theta", "Scattering angle Θ (°)", self.hrs_theta),
            ("theta_start", "Start angle (°)", self.hrs_theta_start),
            ("theta_end", "End angle (°)", self.hrs_theta_end),
            ("n_points", "Number of points", self.hrs_n_points),
        ]):
            lbl = ttk.Label(self.hrs_params_frame, text=label)
            lbl.grid(row=i, column=0, sticky=tk.W, pady=2)
            ent = ttk.Entry(self.hrs_params_frame, textvariable=var, width=12)
            ent.grid(row=i, column=1, sticky=tk.W, padx=6, pady=2)
            self.hrs_param_widgets[key] = (lbl, ent)

        self.btn_run_hrs = ttk.Button(left, text="▶  Run HRS calculation",
                                      command=lambda: self._run_calculation("HRS"))
        self.btn_run_hrs.pack(fill=tk.X, pady=(10, 3), ipady=3)

        right = ttk.LabelFrame(parent, text="Hyperpolarizability tensor (β)", padding=4)
        right.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        self.txt_beta_hrs = scrolledtext.ScrolledText(right, height=12, font=("Consolas", 9))
        self.txt_beta_hrs.pack(fill=tk.BOTH, expand=True, pady=(0, 4))
        self.txt_beta_hrs.insert(tk.END, EXAMPLE_BETA)
        bf = ttk.Frame(right)
        bf.pack(fill=tk.X)
        ttk.Button(bf, text="Load example",
                   command=lambda: self._load_example_beta(self.txt_beta_hrs)).pack(side=tk.LEFT, padx=2)
        ttk.Button(bf, text="Load file…",
                   command=lambda: self._load_file_into(self.txt_beta_hrs)).pack(side=tk.LEFT, padx=2)
        ttk.Button(bf, text="Load Gaussian file…",
                   command=lambda: self._load_gaussian_beta(self.txt_beta_hrs)).pack(side=tk.LEFT, padx=2)

    def _update_hrs_params(self):
        kw = self.hrs_keyword.get()
        show = {
            "wavelength": True,
            "n_real": True,
            "n_imag": True,
            "theta": kw == "polarplot_single",
            "theta_start": kw == "polarplot_integrate",
            "theta_end": kw == "polarplot_integrate",
            "n_points": kw in ("polarplot_integrate", "angle_scattering"),
        }
        for key, (lbl, ent) in self.hrs_param_widgets.items():
            if show.get(key, False):
                lbl.grid()
                ent.grid()
            else:
                lbl.grid_remove()
                ent.grid_remove()


    def _build_shs_tab(self, parent):
        # 3 zones: parameters (scrollable) | files | buttons pinned at bottom of left column
        outer = ttk.Frame(parent)
        outer.pack(fill=tk.BOTH, expand=True)

        # Left column: scrollable parameters + fixed buttons at bottom
        left_col = ttk.Frame(outer, width=300)
        left_col.pack(side=tk.LEFT, fill=tk.Y, padx=(0, 6))
        left_col.pack_propagate(False)

        # Scrollable canvas for parameters
        canvas = tk.Canvas(left_col, highlightthickness=0, width=290)
        scroll = ttk.Scrollbar(left_col, orient=tk.VERTICAL, command=canvas.yview)
        canvas.configure(yscrollcommand=scroll.set)
        scroll.pack(side=tk.RIGHT, fill=tk.Y)
        canvas.pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        params_inner = ttk.Frame(canvas, padding=4)
        win = canvas.create_window((0, 0), window=params_inner, anchor="nw")

        def _cfg(e):
            canvas.configure(scrollregion=canvas.bbox("all"))
            canvas.itemconfig(win, width=e.width)
        params_inner.bind("<Configure>", lambda e: canvas.configure(scrollregion=canvas.bbox("all")))
        canvas.bind("<Configure>", _cfg)

        # Parameter content
        lf_kw = ttk.LabelFrame(params_inner, text="Calculation type", padding=6)
        lf_kw.pack(fill=tk.X, pady=(0, 6))
        cb = ttk.Combobox(lf_kw, textvariable=self.shs_keyword,
                          values=["polarplot_single", "polarplot_integrate", "angle_scattering"],
                          state="readonly", width=24)
        cb.pack(fill=tk.X)
        cb.bind("<<ComboboxSelected>>", lambda e: self._update_shs_params())

        self.shs_params_frame = ttk.LabelFrame(params_inner, text="Parameters", padding=6)
        self.shs_params_frame.pack(fill=tk.X, pady=(0, 6))
        self.shs_param_widgets = {}
        for i, (key, label, var) in enumerate([
            ("n_dipoles", "Number of dipoles", self.shs_n_dipoles),
            ("wavelength", "Wavelength (nm)", self.shs_wavelength),
            ("n_real", "Refractive index – real part", self.shs_n_real),
            ("n_imag", "Refractive index – imaginary part", self.shs_n_imag),
            ("theta", "Scattering angle Θ (°)", self.shs_theta),
            ("theta_start", "Start angle (°)", self.shs_theta_start),
            ("theta_end", "End angle (°)", self.shs_theta_end),
            ("n_points", "Number of points", self.shs_n_points),
        ]):
            lbl = ttk.Label(self.shs_params_frame, text=label)
            lbl.grid(row=i, column=0, sticky=tk.W, pady=2)
            ent = ttk.Entry(self.shs_params_frame, textvariable=var, width=11)
            ent.grid(row=i, column=1, sticky=tk.W, padx=4, pady=2)
            self.shs_param_widgets[key] = (lbl, ent)

        # Buttons ALWAYS visible at bottom of left column
        btn_frame = ttk.Frame(left_col, padding=4)
        btn_frame.pack(side=tk.BOTTOM, fill=tk.X)
        self.btn_run_shs = ttk.Button(btn_frame, text="▶  Run SHS calculation",
                                      command=lambda: self._run_calculation("SHS"))
        self.btn_run_shs.pack(fill=tk.X, pady=(2, 3), ipady=3)

        # Right column: files
        right = ttk.Frame(outer)
        right.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

        nb = ttk.Notebook(right)
        nb.pack(fill=tk.BOTH, expand=True)

        # --- Beta ---
        tab_beta = ttk.Frame(nb, padding=4)
        nb.add(tab_beta, text="  Hyperpolarizability (β)  ")
        self.txt_beta_shs = scrolledtext.ScrolledText(tab_beta, height=10, font=("Consolas", 9))
        self.txt_beta_shs.pack(fill=tk.BOTH, expand=True, pady=(0, 4))
        self.txt_beta_shs.insert(tk.END, EXAMPLE_BETA)
        bf = ttk.Frame(tab_beta)
        bf.pack(fill=tk.X)
        ttk.Button(bf, text="Load example",
                   command=lambda: self._load_example_beta(self.txt_beta_shs)).pack(side=tk.LEFT, padx=2)
        ttk.Button(bf, text="Load file…",
                   command=lambda: self._load_file_into(self.txt_beta_shs)).pack(side=tk.LEFT, padx=2)

        # --- Orientation ---
        tab_orient = ttk.Frame(nb, padding=4)
        nb.add(tab_orient, text="  Positions / Orientations  ")

        # 1) Generation block FIRST (side=BOTTOM) so it stays visible
        gen = ttk.LabelFrame(tab_orient, text="Generate positions / orientations file", padding=8)
        gen.pack(side=tk.BOTTOM, fill=tk.X, pady=(6, 0))

        ttk.Label(gen, text="Choose a script:").pack(anchor=tk.W)
        self.gen_combo = ttk.Combobox(
            gen, textvariable=self.gen_script_choice,
            values=list(GEN_SCRIPTS.keys()), state="readonly", width=40
        )
        self.gen_combo.pack(fill=tk.X, pady=(2, 4))
        self.gen_combo.bind("<<ComboboxSelected>>", lambda e: self._on_gen_script_changed())

        self.gen_params_frame = ttk.Frame(gen)
        self.gen_params_frame.pack(fill=tk.X, pady=2)

        # Custom script area (hidden by default)
        self.custom_frame = ttk.Frame(gen)
        row1 = ttk.Frame(self.custom_frame)
        row1.pack(fill=tk.X)
        ttk.Label(row1, text="Script file:").pack(side=tk.LEFT)
        ttk.Entry(row1, textvariable=self.custom_script_path, width=28).pack(side=tk.LEFT, padx=3)
        ttk.Button(row1, text="Browse…",
                   command=self._browse_custom_script).pack(side=tk.LEFT)
        ttk.Label(self.custom_frame, text="stdin answers (one per line):").pack(anchor=tk.W, pady=(4, 0))
        self.txt_custom_stdin = scrolledtext.ScrolledText(self.custom_frame, height=2, width=40, font=("Consolas", 9))
        self.txt_custom_stdin.pack(fill=tk.X, pady=2)

        self.btn_generate = ttk.Button(
            gen,
            text="▶  Generate positions / orientations file",
            command=self._run_generation_script
        )
        self.btn_generate.pack(fill=tk.X, pady=(6, 0), ipady=5)

        # 2) Load / save buttons (just above generation block)
        btn_o = ttk.Frame(tab_orient)
        btn_o.pack(side=tk.BOTTOM, fill=tk.X, pady=(0, 2))
        ttk.Button(btn_o, text="Load file…", command=self._load_orient_file).pack(side=tk.LEFT, padx=2)
        ttk.Button(btn_o, text="Save as…", command=self._save_orient_file).pack(side=tk.LEFT, padx=2)

        # 3) Text area (remaining space at top)
        self.txt_orient = scrolledtext.ScrolledText(tab_orient, height=6, font=("Consolas", 9))
        self.txt_orient.pack(side=tk.TOP, fill=tk.BOTH, expand=True, pady=(0, 4))

    def _update_shs_params(self):
        kw = self.shs_keyword.get()
        always = {"n_dipoles", "wavelength", "n_real", "n_imag"}
        show = set(always)
        if kw == "polarplot_single":
            show.add("theta")
        elif kw == "polarplot_integrate":
            show.update({"theta_start", "theta_end", "n_points"})
        elif kw == "angle_scattering":
            show.add("n_points")
        for key, (lbl, ent) in self.shs_param_widgets.items():
            (lbl.grid() if key in show else lbl.grid_remove())
            (ent.grid() if key in show else ent.grid_remove())

    def _on_gen_script_changed(self):
        # Clear previous parameter widgets
        for w in self.gen_params_frame.winfo_children():
            w.destroy()
        self.gen_vars.clear()
        self.custom_frame.pack_forget()

        choice = self.gen_script_choice.get()
        script_file, fields = GEN_SCRIPTS.get(choice, (None, []))

        if choice == "Custom script…":
            # Show before the Generate button
            self.custom_frame.pack(fill=tk.X, pady=4, before=self.btn_generate)
            return

        for i, (label, key, default) in enumerate(fields):
            ttk.Label(self.gen_params_frame, text=label).grid(
                row=i // 2, column=(i % 2) * 2, sticky=tk.W, padx=2, pady=2
            )
            var = tk.StringVar(value=default)
            self.gen_vars[key] = var
            ttk.Entry(self.gen_params_frame, textvariable=var, width=10).grid(
                row=i // 2, column=(i % 2) * 2 + 1, sticky=tk.W, padx=2, pady=2
            )

    def _browse_custom_script(self):
        path = filedialog.askopenfilename(
            title="Choose a generation Python script",
            filetypes=[("Python", "*.py"), ("All", "*.*")]
        )
        if path:
            self.custom_script_path.set(path)

    # ------------------------------------------------------------------
    # Sphere (colloidal continuum model)
    # ------------------------------------------------------------------
    def _build_sphere_tab(self, parent):
        left = ttk.Frame(parent)
        left.pack(side=tk.LEFT, fill=tk.Y, padx=(0, 8))
        right = ttk.Frame(parent)
        right.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

        lf_kw = ttk.LabelFrame(left, text="Calculation type", padding=6)
        lf_kw.pack(fill=tk.X, pady=4)
        ttk.Combobox(
            lf_kw, textvariable=self.sph_keyword,
            values=["polarplot_single", "polarplot_integrate", "angle_scattering"],
            state="readonly", width=22
        ).pack(fill=tk.X)
        self.sph_keyword.trace_add("write", lambda *_: self._update_sph_params())

        self.sph_params_frame = ttk.LabelFrame(left, text="Parameters", padding=6)
        self.sph_params_frame.pack(fill=tk.X, pady=4)

        btn = ttk.Frame(left)
        btn.pack(fill=tk.X, pady=8)
        self.btn_run_sphere = ttk.Button(
            btn, text="▶  Run sphere SHS calculation (NLRGD theory)",
            command=lambda: self._run_calculation("Sphere")
        )
        self.btn_run_sphere.pack(fill=tk.X, pady=2)
        self.btn_run_sphere_mie = ttk.Button(
            btn, text="▶  Run sphere SHS calculation (Mie theory)",
            command=lambda: self._run_calculation("SphereMie")
        )
        self.btn_run_sphere_mie.pack(fill=tk.X, pady=2)

        lf_b = ttk.LabelFrame(right, text="Surface chi(2) / hyperpolarizability tensor", padding=4)
        lf_b.pack(fill=tk.BOTH, expand=True)

        ttk.Label(
            lf_b,
            text="Choose a sub-tab below: full 27-component tensor, or Mie isotropic (rrr, rtt, trt).",
            foreground="#334155",
        ).pack(anchor=tk.W, pady=(0, 4))

        self.sph_beta_nb = ttk.Notebook(lf_b)
        self.sph_beta_nb.pack(fill=tk.BOTH, expand=True)

        # ----- Tab 1: full tensor -----
        tab_full = ttk.Frame(self.sph_beta_nb, padding=4)
        self.sph_beta_nb.add(tab_full, text="Full tensor")
        self.txt_beta_sph = scrolledtext.ScrolledText(tab_full, height=14, font=("Consolas", 9))
        self.txt_beta_sph.pack(fill=tk.BOTH, expand=True)
        self.txt_beta_sph.insert(tk.END, EXAMPLE_BETA)
        bf = ttk.Frame(tab_full)
        bf.pack(fill=tk.X, pady=4)
        ttk.Button(bf, text="Load example",
                   command=lambda: self._load_example_beta(self.txt_beta_sph)).pack(side=tk.LEFT, padx=2)
        ttk.Button(bf, text="Load file...",
                   command=lambda: self._load_file_into(self.txt_beta_sph)).pack(side=tk.LEFT, padx=2)

        # ----- Tab 2: Mie isotropic -----
        tab_mie = ttk.Frame(self.sph_beta_nb, padding=8)
        self.sph_beta_nb.add(tab_mie, text="Mie isotropic (rrr, rtt, trt)")
        ttk.Label(
            tab_mie,
            text="Isotropic achiral surface. Editing these fields updates the Full tensor tab:",
            wraplength=380, justify=tk.LEFT,
        ).pack(anchor=tk.W, pady=(0, 4))
        ttk.Label(
            tab_mie,
            text="zzz = chi_rrr ;  zxx = zyy = chi_rtt ;  xzx = xxz = yzy = yyz = chi_trt",
            font=("Consolas", 8), foreground="#475569",
        ).pack(anchor=tk.W, pady=(0, 10))

        self.sph_chi_rrr = tk.StringVar(value="1.0")
        self.sph_chi_rtt = tk.StringVar(value="0.0")
        self.sph_chi_trt = tk.StringVar(value="0.0")
        form = ttk.Frame(tab_mie)
        form.pack(fill=tk.X)
        for i, (lab, var, hint) in enumerate((
            ("chi_rrr  (perp-perp-perp)", self.sph_chi_rrr, "beta_zzz"),
            ("chi_rtt  (perp-par-par)", self.sph_chi_rtt, "beta_zxx = beta_zyy = chi_rtt"),
            ("chi_trt  (par-perp-par)", self.sph_chi_trt, "beta_xzx = xxz = yzy = yyz = chi_trt"),
        )):
            ttk.Label(form, text=lab).grid(row=i, column=0, sticky=tk.W, pady=4)
            ttk.Entry(form, textvariable=var, width=14).grid(row=i, column=1, sticky=tk.W, padx=8)
            ttk.Label(form, text=hint, foreground="#64748b").grid(row=i, column=2, sticky=tk.W)

        bf2 = ttk.Frame(tab_mie)
        bf2.pack(fill=tk.X, pady=(14, 4))
        ttk.Button(bf2, text="Apply -> Full tensor",
                   command=self._apply_mie_chi_to_beta).pack(side=tk.LEFT, padx=2)
        ttk.Button(bf2, text="Read <- Full tensor",
                   command=self._read_mie_chi_from_beta).pack(side=tk.LEFT, padx=2)

        for v in (self.sph_chi_rrr, self.sph_chi_rtt, self.sph_chi_trt):
            v.trace_add("write", lambda *_: self._apply_mie_chi_to_beta())

    def _update_sph_params(self):
        for w in self.sph_params_frame.winfo_children():
            w.destroy()
        rows = [
            ("Wavelength (nm)", self.sph_wavelength),
            ("n_p (particle) – real part", self.sph_np_real),
            ("n_p (particle) – imaginary part", self.sph_np_imag),
            ("n_m (medium) – real part", self.sph_nm_real),
            ("n_m (medium) – imaginary part", self.sph_nm_imag),
            ("Radius (nm)", self.sph_radius),
            ("Intensity scale factor", self.sph_i_scale),
        ]
        kw = self.sph_keyword.get()
        if kw == "polarplot_single":
            rows.append(("Scattering angle Θ (°)", self.sph_theta))
        elif kw == "polarplot_integrate":
            rows += [
                ("Start angle (°)", self.sph_theta_start),
                ("End angle (°)", self.sph_theta_end),
                ("Number of points", self.sph_n_points),
            ]
        else:
            rows.append(("Number of points", self.sph_n_points))
        for i, (lab, var) in enumerate(rows):
            ttk.Label(self.sph_params_frame, text=lab).grid(row=i, column=0, sticky=tk.W, pady=2)
            ttk.Entry(self.sph_params_frame, textvariable=var, width=12).grid(row=i, column=1, sticky=tk.W, padx=4)

    # ------------------------------------------------------------------
    # Cylinder (colloidal continuum model)
    # ------------------------------------------------------------------
    def _build_cylinder_tab(self, parent):
        left = ttk.Frame(parent)
        left.pack(side=tk.LEFT, fill=tk.Y, padx=(0, 8))
        right = ttk.Frame(parent)
        right.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

        lf_kw = ttk.LabelFrame(left, text="Calculation type", padding=6)
        lf_kw.pack(fill=tk.X, pady=4)
        ttk.Combobox(
            lf_kw, textvariable=self.cyl_keyword,
            values=["polarplot_single", "polarplot_integrate", "angle_scattering"],
            state="readonly", width=22
        ).pack(fill=tk.X)
        self.cyl_keyword.trace_add("write", lambda *_: self._update_cyl_params())

        self.cyl_params_frame = ttk.LabelFrame(left, text="Parameters", padding=6)
        self.cyl_params_frame.pack(fill=tk.X, pady=4)

        btn = ttk.Frame(left)
        btn.pack(fill=tk.X, pady=8)
        self.btn_run_cylinder = ttk.Button(
            btn, text="▶  Run Cylinder calculation",
            command=lambda: self._run_calculation("Cylinder")
        )
        self.btn_run_cylinder.pack(fill=tk.X, pady=2)

        lf_b = ttk.LabelFrame(right, text="Surface χ⁽²⁾ / hyperpolarizability tensor", padding=4)
        lf_b.pack(fill=tk.BOTH, expand=True)
        self.txt_beta_cyl = scrolledtext.ScrolledText(lf_b, height=16, font=("Consolas", 9))
        self.txt_beta_cyl.pack(fill=tk.BOTH, expand=True)
        self.txt_beta_cyl.insert(tk.END, EXAMPLE_BETA)
        bf = ttk.Frame(lf_b)
        bf.pack(fill=tk.X, pady=4)
        ttk.Button(bf, text="Load example",
                   command=lambda: self._load_example_beta(self.txt_beta_cyl)).pack(side=tk.LEFT, padx=2)
        ttk.Button(bf, text="Load file…",
                   command=lambda: self._load_file_into(self.txt_beta_cyl)).pack(side=tk.LEFT, padx=2)

    def _update_cyl_params(self):
        for w in self.cyl_params_frame.winfo_children():
            w.destroy()
        rows = [
            ("Wavelength (nm)", self.cyl_wavelength),
            ("n_p (particle) – real part", self.cyl_np_real),
            ("n_p (particle) – imaginary part", self.cyl_np_imag),
            ("n_m (medium) – real part", self.cyl_nm_real),
            ("n_m (medium) – imaginary part", self.cyl_nm_imag),
            ("Radius (nm)", self.cyl_radius),
            ("Length (nm)", self.cyl_length),
        ]
        kw = self.cyl_keyword.get()
        if kw == "polarplot_single":
            rows.append(("Scattering angle Θ (°)", self.cyl_theta))
        elif kw == "polarplot_integrate":
            rows += [
                ("Start angle (°)", self.cyl_theta_start),
                ("End angle (°)", self.cyl_theta_end),
                ("Number of points", self.cyl_n_points),
            ]
        else:
            rows.append(("Number of points", self.cyl_n_points))
        for i, (lab, var) in enumerate(rows):
            ttk.Label(self.cyl_params_frame, text=lab).grid(row=i, column=0, sticky=tk.W, pady=2)
            ttk.Entry(self.cyl_params_frame, textvariable=var, width=12).grid(row=i, column=1, sticky=tk.W, padx=4)

    # ------------------------------------------------------------------
    # Results + plots (intensity + assembly)
    # ------------------------------------------------------------------

    # ------------------------------------------------------------------
    # Ellipsoid of revolution (colloidal continuum model)
    # ------------------------------------------------------------------
    def _build_ellipsoid_tab(self, parent):
        left = ttk.Frame(parent)
        left.pack(side=tk.LEFT, fill=tk.Y, padx=(0, 8))
        right = ttk.Frame(parent)
        right.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

        lf_kw = ttk.LabelFrame(left, text="Calculation type", padding=6)
        lf_kw.pack(fill=tk.X, pady=4)
        ttk.Combobox(
            lf_kw, textvariable=self.ell_keyword,
            values=["polarplot_single", "polarplot_integrate", "angle_scattering"],
            state="readonly", width=22
        ).pack(fill=tk.X)
        self.ell_keyword.trace_add("write", lambda *_: self._update_ell_params())

        self.ell_params_frame = ttk.LabelFrame(left, text="Parameters", padding=6)
        self.ell_params_frame.pack(fill=tk.X, pady=4)

        btn = ttk.Frame(left)
        btn.pack(fill=tk.X, pady=8)
        self.btn_run_ellipsoid = ttk.Button(
            btn, text="▶  Run Ellipsoid calculation",
            command=lambda: self._run_calculation("Ellipsoid")
        )
        self.btn_run_ellipsoid.pack(fill=tk.X, pady=2)

        lf_b = ttk.LabelFrame(right, text="Surface chi(2) / hyperpolarizability tensor", padding=4)
        lf_b.pack(fill=tk.BOTH, expand=True)
        self.txt_beta_ell = scrolledtext.ScrolledText(lf_b, height=16, font=("Consolas", 9))
        self.txt_beta_ell.pack(fill=tk.BOTH, expand=True)
        self.txt_beta_ell.insert(tk.END, EXAMPLE_BETA)
        bf = ttk.Frame(lf_b)
        bf.pack(fill=tk.X, pady=4)
        ttk.Button(bf, text="Load example",
                   command=lambda: self._load_example_beta(self.txt_beta_ell)).pack(side=tk.LEFT, padx=2)
        ttk.Button(bf, text="Load file...",
                   command=lambda: self._load_file_into(self.txt_beta_ell)).pack(side=tk.LEFT, padx=2)

        self._update_ell_params()

    def _update_ell_params(self):
        for w in self.ell_params_frame.winfo_children():
            w.destroy()
        rows = [
            ("Wavelength (nm)", self.ell_wavelength),
            ("n_p (particle) – real part", self.ell_np_real),
            ("n_p (particle) – imaginary part", self.ell_np_imag),
            ("n_m (medium) – real part", self.ell_nm_real),
            ("n_m (medium) – imaginary part", self.ell_nm_imag),
            ("Equatorial radius a (nm)", self.ell_radius),
            ("Length along axis 2c (nm)", self.ell_length),
        ]
        kw = self.ell_keyword.get()
        if kw == "polarplot_single":
            rows.append(("Scattering angle Theta (deg)", self.ell_theta))
        elif kw == "polarplot_integrate":
            rows += [
                ("Theta start (deg)", self.ell_theta_start),
                ("Theta end (deg)", self.ell_theta_end),
                ("Number of points", self.ell_n_points),
            ]
        else:
            rows.append(("Number of angle points", self.ell_n_points))
        for i, (lab, var) in enumerate(rows):
            ttk.Label(self.ell_params_frame, text=lab).grid(row=i, column=0, sticky=tk.W, pady=2)
            ttk.Entry(self.ell_params_frame, textvariable=var, width=12).grid(row=i, column=1, sticky=tk.W, padx=4)


    # ------------------------------------------------------------------
    # Cube / rectangular box (colloidal continuum model)
    # ------------------------------------------------------------------
    def _build_cube_tab(self, parent):
        left = ttk.Frame(parent)
        left.pack(side=tk.LEFT, fill=tk.Y, padx=(0, 8))
        right = ttk.Frame(parent)
        right.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

        lf_kw = ttk.LabelFrame(left, text="Calculation type", padding=6)
        lf_kw.pack(fill=tk.X, pady=4)
        ttk.Combobox(
            lf_kw, textvariable=self.cube_keyword,
            values=["polarplot_single", "polarplot_integrate", "angle_scattering"],
            state="readonly", width=22
        ).pack(fill=tk.X)
        self.cube_keyword.trace_add("write", lambda *_: self._update_cube_params())

        self.cube_params_frame = ttk.LabelFrame(left, text="Parameters", padding=6)
        self.cube_params_frame.pack(fill=tk.X, pady=4)

        btn = ttk.Frame(left)
        btn.pack(fill=tk.X, pady=8)
        self.btn_run_cube = ttk.Button(
            btn, text="▶  Run Cube calculation",
            command=lambda: self._run_calculation("Cube")
        )
        self.btn_run_cube.pack(fill=tk.X, pady=2)

        lf_b = ttk.LabelFrame(right, text="Surface chi(2) / hyperpolarizability tensor", padding=4)
        lf_b.pack(fill=tk.BOTH, expand=True)
        self.txt_beta_cube = scrolledtext.ScrolledText(lf_b, height=16, font=("Consolas", 9))
        self.txt_beta_cube.pack(fill=tk.BOTH, expand=True)
        self.txt_beta_cube.insert(tk.END, EXAMPLE_BETA)
        bf = ttk.Frame(lf_b)
        bf.pack(fill=tk.X, pady=4)
        ttk.Button(bf, text="Load example",
                   command=lambda: self._load_example_beta(self.txt_beta_cube)).pack(side=tk.LEFT, padx=2)
        ttk.Button(bf, text="Load file...",
                   command=lambda: self._load_file_into(self.txt_beta_cube)).pack(side=tk.LEFT, padx=2)

        self._update_cube_params()

    def _update_cube_params(self):
        for w in self.cube_params_frame.winfo_children():
            w.destroy()
        rows = [
            ("Wavelength (nm)", self.cube_wavelength),
            ("n_p (particle) – real part", self.cube_np_real),
            ("n_p (particle) – imaginary part", self.cube_np_imag),
            ("n_m (medium) – real part", self.cube_nm_real),
            ("n_m (medium) – imaginary part", self.cube_nm_imag),
            ("Edge a along x (nm)", self.cube_a),
            ("Edge b along y (nm)", self.cube_b),
            ("Edge c along z (nm)", self.cube_c),
        ]
        kw = self.cube_keyword.get()
        if kw == "polarplot_single":
            rows.append(("Scattering angle Theta (deg)", self.cube_theta))
        elif kw == "polarplot_integrate":
            rows += [
                ("Theta start (deg)", self.cube_theta_start),
                ("Theta end (deg)", self.cube_theta_end),
                ("Number of points", self.cube_n_points),
            ]
        else:
            rows.append(("Number of angle points", self.cube_n_points))
        for i, (lab, var) in enumerate(rows):
            ttk.Label(self.cube_params_frame, text=lab).grid(row=i, column=0, sticky=tk.W, pady=2)
            ttk.Entry(self.cube_params_frame, textvariable=var, width=12).grid(row=i, column=1, sticky=tk.W, padx=4)

    def _build_results(self, parent):
        paned = tk.PanedWindow(parent, orient=tk.HORIZONTAL, sashrelief=tk.RAISED)
        paned.pack(fill=tk.BOTH, expand=True)

        left = ttk.LabelFrame(paned, text="Calculation results", padding=4)
        paned.add(left, width=480)
        self.txt_result = scrolledtext.ScrolledText(
            left, height=10, font=("Consolas", 9),
            bg="#0f172a", fg="#e2e8f0", insertbackground="white"
        )
        self.txt_result.pack(fill=tk.BOTH, expand=True)
        btn_r = ttk.Frame(left)
        btn_r.pack(fill=tk.X, pady=3)
        ttk.Button(btn_r, text="Save…", command=self._save_results).pack(side=tk.LEFT, padx=2)
        ttk.Button(btn_r, text="Clear",
                   command=lambda: self.txt_result.delete(1.0, tk.END)).pack(side=tk.LEFT, padx=2)
        self.btn_cancel_calc = ttk.Button(btn_r, text="⏹ Stop calculation",
                                          command=self._cancel_calculation, state=tk.DISABLED)
        self.btn_cancel_calc.pack(side=tk.RIGHT, padx=2)

        right = ttk.LabelFrame(paned, text="Plots", padding=4)
        paned.add(right)

        # Export buttons FIRST at bottom so they stay visible under the graphs
        exp = ttk.Frame(right)
        exp.pack(side=tk.BOTTOM, fill=tk.X, pady=(6, 2))
        ttk.Button(exp, text="Export intensity plot…",
                   command=lambda: self._export_fig(self._last_intensity_fig, "intensity")).pack(side=tk.LEFT, padx=3)
        ttk.Button(exp, text="Export assembly plot…",
                   command=lambda: self._export_fig(self._last_assembly_fig, "assembly")).pack(side=tk.LEFT, padx=3)

        # Notebook: Intensity | Assembly (fills remaining space above the buttons)
        self.plot_nb = ttk.Notebook(right)
        self.plot_nb.pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        self.tab_plot_intensity = ttk.Frame(self.plot_nb)
        self.tab_plot_assembly = ttk.Frame(self.plot_nb)
        self.plot_nb.add(self.tab_plot_intensity, text="  SHS Intensity  ")
        self.plot_nb.add(self.tab_plot_assembly, text="  Assembly / Geometry  ")

        self.plot_intensity_container = ttk.Frame(self.tab_plot_intensity)
        self.plot_intensity_container.pack(fill=tk.BOTH, expand=True)
        self.plot_assembly_container = ttk.Frame(self.tab_plot_assembly)
        self.plot_assembly_container.pack(fill=tk.BOTH, expand=True)

    # ==================================================================
    def _log(self, msg: str):
        self.txt_result.insert(tk.END, msg + "\n")
        self.txt_result.see(tk.END)
        self.update_idletasks()

    def _log_live(self, msg: str):
        """Immediately show a progress line (called on UI thread via after)."""
        self.txt_result.insert(tk.END, msg + "\n")
        self.txt_result.see(tk.END)
        # force refresh even during a long calculation
        try:
            self.update_idletasks()
        except tk.TclError:
            pass

    def _parse_beta_text(self, content: str) -> dict:
        """Parse 'ijk value' lines into a dict of lowercase labels."""
        beta = {a + b + c: 0.0 for a in "xyz" for b in "xyz" for c in "xyz"}
        for line in content.splitlines():
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) >= 2:
                try:
                    beta[parts[0].lower()] = float(parts[1])
                except ValueError:
                    continue
        return beta

    def _format_beta_text(self, beta: dict) -> str:
        lines = ["# Surface χ⁽²⁾ / hyperpolarizability tensor (a.u.)"]
        for a in "xyz":
            for b in "xyz":
                for c in "xyz":
                    key = a + b + c
                    lines.append(f"{key} {beta.get(key, 0.0)}")
        return "\n".join(lines) + "\n"

    def _apply_mie_chi_to_beta(self):
        """
        Map isotropic Mie surface components onto the full 27-tensor:
          β_zzz = χ_rrr
          β_zxx = β_zyy = χ_rtt
          β_xzx = β_xxz = β_yzy = β_yyz = χ_trt
        Other components are zeroed (isotropic achiral surface).
        """
        if not hasattr(self, "txt_beta_sph"):
            return
        try:
            chi_rrr = float(self.sph_chi_rrr.get())
            chi_rtt = float(self.sph_chi_rtt.get())
            chi_trt = float(self.sph_chi_trt.get())
        except ValueError:
            return
        beta = {a + b + c: 0.0 for a in "xyz" for b in "xyz" for c in "xyz"}
        beta["zzz"] = chi_rrr
        beta["zxx"] = chi_rtt
        beta["zyy"] = chi_rtt
        for key in ("xzx", "xxz", "yzy", "yyz"):
            beta[key] = chi_trt
        text_beta = self._format_beta_text(beta)
        self.txt_beta_sph.delete(1.0, tk.END)
        self.txt_beta_sph.insert(tk.END, text_beta)

    def _read_mie_chi_from_beta(self):
        """Extract χ_rrr, χ_rtt, χ_trt from the full tensor editor."""
        if not hasattr(self, "txt_beta_sph"):
            return
        beta = self._parse_beta_text(self.txt_beta_sph.get(1.0, tk.END))
        chi_rrr = beta.get("zzz", 0.0)
        # Independent isotropic components (not sums): zxx = zyy = chi_rtt
        zxx, zyy = beta.get("zxx", 0.0), beta.get("zyy", 0.0)
        chi_rtt = 0.5 * (zxx + zyy) if (zxx or zyy) else 0.0
        comps = [beta.get(k, 0.0) for k in ("xzx", "xxz", "yzy", "yyz")]
        chi_trt = sum(comps) / 4.0 if any(comps) else 0.0
        self.sph_chi_rrr.set(str(chi_rrr))
        self.sph_chi_rtt.set(str(chi_rtt))
        self.sph_chi_trt.set(str(chi_trt))

    def _load_example_beta(self, w):

        w.delete(1.0, tk.END)
        w.insert(tk.END, EXAMPLE_BETA)

    def _load_file_into(self, w):
        path = filedialog.askopenfilename(filetypes=[("Text", "*.txt"), ("All", "*.*")])
        if path:
            with open(path, encoding="utf-8", errors="replace") as f:
                w.delete(1.0, tk.END)
                w.insert(tk.END, f.read())

    def _load_gaussian_beta(self, w):
        """Load a Gaussian .log file and convert Beta(-2w;w,w) to PySHS beta format."""
        path = filedialog.askopenfilename(
            title="Select Gaussian log file",
            filetypes=[
                ("Gaussian log", "*.log"),
                ("Output files", "*.out"),
                ("All files", "*.*"),
            ],
        )
        if not path:
            return

        script = SRC_DIR / "script_gauss2PySHS.py"
        if not script.exists():
            messagebox.showerror("Error", f"Script not found:\n{script}")
            return

        out_file = WORK_DIR / f"beta_from_gauss_{os.getpid()}.txt"
        self._log(f"Converting Gaussian log → PySHS beta:\n  {path}")

        def worker():
            try:
                proc = subprocess.run(
                    [sys.executable, str(script), path, str(out_file)],
                    capture_output=True,
                    text=True,
                    timeout=60,
                    cwd=str(WORK_DIR),
                )
                if proc.returncode != 0:
                    err = (proc.stderr or proc.stdout or "unknown error").strip()
                    self.after(0, lambda: messagebox.showerror(
                        "Error", f"Gaussian conversion failed:\n{err}"
                    ))
                    self.after(0, lambda: self._log(f"Gaussian conversion failed:\n{err}"))
                    return
                if not out_file.exists():
                    self.after(0, lambda: messagebox.showerror(
                        "Error", "Conversion script did not produce an output file."
                    ))
                    return
                content = out_file.read_text(encoding="utf-8", errors="replace")
                try:
                    out_file.unlink()
                except OSError:
                    pass

                def ui():
                    w.delete(1.0, tk.END)
                    w.insert(tk.END, content)
                    self._log("Hyperpolarizability tensor loaded from Gaussian log.")
                self.after(0, ui)
            except Exception as e:
                self.after(0, lambda: messagebox.showerror("Error", str(e)))
                self.after(0, lambda: self._log(f"Exception: {e}"))

        threading.Thread(target=worker, daemon=True).start()

    def _count_orient_lines(self, content: str) -> int:
        """Count valid orientation rows (at least 6 numeric fields)."""
        n = 0
        for line in content.splitlines():
            parts = line.split()
            if len(parts) < 6:
                continue
            try:
                [float(x) for x in parts[:6]]
                n += 1
            except ValueError:
                continue
        return n

    def _sync_n_dipoles_from_orient(self):
        """Set Number of dipoles field from orientation text."""
        content = self.txt_orient.get(1.0, tk.END)
        n = self._count_orient_lines(content)
        if n > 0:
            self.shs_n_dipoles.set(str(n))
            return n
        return 0

    def _load_orient_file(self):
        path = filedialog.askopenfilename(title="Load positions/orientations",
                                          filetypes=[("Text", "*.txt"), ("All", "*.*")])
        if path:
            with open(path, encoding="utf-8", errors="replace") as f:
                self.txt_orient.delete(1.0, tk.END)
                self.txt_orient.insert(tk.END, f.read())
            self._log(f"Orientation file loaded: {path}")
            n = self._sync_n_dipoles_from_orient()
            self._log(f"Number of dipoles set to {n} (from file).")
            self._show_assembly_from_text(self.txt_orient.get(1.0, tk.END))

    def _save_orient_file(self):
        content = self.txt_orient.get(1.0, tk.END).strip()
        if not content:
            messagebox.showinfo("Info", "Text area is empty.")
            return
        path = filedialog.asksaveasfilename(defaultextension=".txt", filetypes=[("Text", "*.txt")])
        if path:
            with open(path, "w", encoding="utf-8") as f:
                f.write(content + "\n")
            self._log(f"Orientation file saved → {path}")

    def _save_results(self):
        content = self.txt_result.get(1.0, tk.END)
        if not content.strip():
            return
        path = filedialog.asksaveasfilename(defaultextension=".txt", filetypes=[("Text", "*.txt")])
        if path:
            with open(path, "w", encoding="utf-8") as f:
                f.write(content)

    def _export_fig(self, fig, default_name):
        if fig is None:
            messagebox.showinfo("Info", "No plot to export yet.")
            return
        path = filedialog.asksaveasfilename(
            defaultextension=".png",
            initialfile=f"{default_name}.png",
            filetypes=[("PNG", "*.png"), ("PDF", "*.pdf")]
        )
        if path:
            fig.savefig(path, dpi=150)
            self._log(f"Plot exported → {path}")

    # ------------------------------------------------------------------
    # Multi-script generation
    # ------------------------------------------------------------------
    def _run_generation_script(self):
        choice = self.gen_script_choice.get()
        out_file = WORK_DIR / f"orient_gen_{os.getpid()}.txt"

        if choice == "Custom script…":
            script = self.custom_script_path.get().strip()
            if not script or not Path(script).exists():
                messagebox.showerror("Error", "Please select a valid script file.")
                return
            stdin_data = self.txt_custom_stdin.get(1.0, tk.END)
            script_path = Path(script)
        else:
            script_name, fields = GEN_SCRIPTS[choice]
            script_path = SRC_DIR / script_name
            if not script_path.exists():
                messagebox.showerror("Error", f"Script not found: {script_path}")
                return
            # Construire stdin dans l'ordre des champs
            answers = []
            try:
                for _, key, _ in fields:
                    answers.append(str(self.gen_vars[key].get()))
            except KeyError:
                messagebox.showerror("Error", "Parameters manquants.")
                return
            stdin_data = "\n".join(answers) + "\n"

        self._log(f"Generating with {script_path.name}…")

        def worker():
            try:
                env = os.environ.copy()
                env["MPLBACKEND"] = "Agg"
                proc = subprocess.run(
                    [sys.executable, str(script_path), str(out_file)],
                    input=stdin_data, capture_output=True, text=True,
                    timeout=90, cwd=str(WORK_DIR), env=env
                )
                if proc.returncode != 0:
                    self.after(0, lambda: self._log(
                        "Error script :\n" + (proc.stderr or proc.stdout or "")
                    ))
                    return
                if not out_file.exists():
                    self.after(0, lambda: self._log("Output file was not created."))
                    return
                content = out_file.read_text(encoding="utf-8", errors="replace")
                try:
                    out_file.unlink()
                except OSError:
                    pass
                for extra in (Path(str(out_file) + ".png"),):
                    if extra.exists():
                        try:
                            extra.unlink()
                        except OSError:
                            pass

                def ui():
                    self.txt_orient.delete(1.0, tk.END)
                    self.txt_orient.insert(tk.END, content)
                    n = self._sync_n_dipoles_from_orient()
                    self._log(f"Generation OK ({n} dipoles). Number of dipoles field updated.")
                    self._show_assembly_from_text(content)
                self.after(0, ui)
            except Exception as e:
                self.after(0, lambda: self._log(f"Exception: {e}"))

        threading.Thread(target=worker, daemon=True).start()

    # ------------------------------------------------------------------
    # Module check (no compilation)
    # ------------------------------------------------------------------
    def _initial_compile(self):
        self._log("Checking pure-Python modules…")
        threading.Thread(target=self._check_all_modules, daemon=True).start()

    def _check_all_modules(self):
        ok1 = check_python_module(SRC_HRS, self._log)
        ok2 = check_python_module(SRC_SHS, self._log)
        ok3 = check_python_module(SRC_SPHERE, self._log)
        ok4 = check_python_module(SRC_CYLINDER, self._log)
        ok5 = check_python_module(SRC_ELLIPSOID, self._log)
        ok6 = check_python_module(SRC_CUBE, self._log)
        ok7 = check_python_module(SRC_MIE, self._log)
        self._log(
            "All Python modules are ready (no compilation needed).\n"
            if (ok1 and ok2 and ok3 and ok4 and ok5 and ok6 and ok7)
            else "Some Python modules are missing.\n"
        )

    # ------------------------------------------------------------------
    # Calcul
    # ------------------------------------------------------------------
    def _run_calculation(self, module: str):
        if self.job_running:
            messagebox.showwarning("Warning", "A calculation is already running.")
            return

        orient_content = None

        if module == "HRS":
            # HRS: wavelength, n_real, n_imag, then angle params
            src = SRC_HRS
            kw = self.hrs_keyword.get()
            beta_content = self.txt_beta_hrs.get(1.0, tk.END).strip()
            try:
                answers = [
                    str(float(self.hrs_wavelength.get())),
                    str(float(self.hrs_n_real.get())),
                    str(float(self.hrs_n_imag.get())),
                ]
                if kw == "polarplot_single":
                    answers.append(str(float(self.hrs_theta.get())))
                elif kw == "polarplot_integrate":
                    answers += [str(float(self.hrs_theta_start.get())),
                                str(float(self.hrs_theta_end.get())),
                                str(int(self.hrs_n_points.get()))]
                else:
                    answers.append(str(int(self.hrs_n_points.get())))
            except ValueError as e:
                messagebox.showerror("Error", f"Invalid value: {e}")
                return

        elif module == "SHS":
            src = SRC_SHS
            kw = self.shs_keyword.get()
            beta_content = self.txt_beta_shs.get(1.0, tk.END).strip()
            orient_content = self.txt_orient.get(1.0, tk.END).strip()
            if not orient_content:
                messagebox.showerror("Error", "Positions/orientations are empty.")
                return
            # CRITICAL: N must match the orientation file (SHS intensity scales as N^2)
            n_from_file = self._count_orient_lines(orient_content)
            if n_from_file <= 0:
                messagebox.showerror("Error", "No valid orientation lines found.")
                return
            try:
                n_field = int(self.shs_n_dipoles.get())
            except ValueError:
                n_field = -1
            if n_field != n_from_file:
                self._log(
                    f"WARNING: Number of dipoles field ({n_field}) != "
                    f"orientation file lines ({n_from_file}). "
                    f"Using {n_from_file} from the file."
                )
                self.shs_n_dipoles.set(str(n_from_file))
            try:
                answers = [
                    str(n_from_file),
                    str(float(self.shs_wavelength.get())),
                    str(float(self.shs_n_real.get())),
                    str(float(self.shs_n_imag.get())),
                ]
                if kw == "polarplot_single":
                    answers.append(str(float(self.shs_theta.get())))
                elif kw == "polarplot_integrate":
                    answers += [str(float(self.shs_theta_start.get())),
                                str(float(self.shs_theta_end.get())),
                                str(int(self.shs_n_points.get()))]
                else:
                    answers.append(str(int(self.shs_n_points.get())))
            except ValueError as e:
                messagebox.showerror("Error", f"Invalid value: {e}")
                return

        elif module == "Sphere":
            # NLWKB sphere: use_nlwkb=1, wavelength, n_p, n_m, radius + angles
            src = SRC_SPHERE
            kw = self.sph_keyword.get()
            beta_content = self.txt_beta_sph.get(1.0, tk.END).strip()
            try:
                answers = [
                    "1",  # always use NLWKB (n_p and n_m)
                    str(float(self.sph_wavelength.get())),
                    str(float(self.sph_np_real.get())),
                    str(float(self.sph_np_imag.get())),
                    str(float(self.sph_nm_real.get())),
                    str(float(self.sph_nm_imag.get())),
                    str(float(self.sph_radius.get())),
                ]
                if kw == "polarplot_single":
                    answers.append(str(float(self.sph_theta.get())))
                elif kw == "polarplot_integrate":
                    answers += [str(float(self.sph_theta_start.get())),
                                str(float(self.sph_theta_end.get())),
                                str(int(self.sph_n_points.get()))]
                else:
                    answers.append(str(int(self.sph_n_points.get())))
            except ValueError as e:
                messagebox.showerror("Error", f"Invalid value: {e}")
                return

        elif module == "SphereMie":
            # Nonlinear Mie theory (Python, Src/shs_mie.py)
            # stdin answers: wavelength, n_p, n_m, radius + angles (no use_nlwkb flag)
            src = SRC_MIE
            kw = self.sph_keyword.get()
            beta_content = self.txt_beta_sph.get(1.0, tk.END).strip()
            try:
                answers = [
                    str(float(self.sph_wavelength.get())),
                    str(float(self.sph_np_real.get())),
                    str(float(self.sph_np_imag.get())),
                    str(float(self.sph_nm_real.get())),
                    str(float(self.sph_nm_imag.get())),
                    str(float(self.sph_radius.get())),
                ]
                if kw == "polarplot_single":
                    answers.append(str(float(self.sph_theta.get())))
                elif kw == "polarplot_integrate":
                    answers += [str(float(self.sph_theta_start.get())),
                                str(float(self.sph_theta_end.get())),
                                str(int(self.sph_n_points.get()))]
                else:
                    answers.append(str(int(self.sph_n_points.get())))
            except ValueError as e:
                messagebox.showerror("Error", f"Invalid value: {e}")
                return

        elif module == "Cylinder":
            # NLWKB cylinder: wavelength, n_p, n_m, R, L + angles
            src = SRC_CYLINDER
            kw = self.cyl_keyword.get()
            beta_content = self.txt_beta_cyl.get(1.0, tk.END).strip()
            try:
                answers = [
                    str(float(self.cyl_wavelength.get())),
                    str(float(self.cyl_np_real.get())),
                    str(float(self.cyl_np_imag.get())),
                    str(float(self.cyl_nm_real.get())),
                    str(float(self.cyl_nm_imag.get())),
                    str(float(self.cyl_radius.get())),
                    str(float(self.cyl_length.get())),
                ]
                if kw == "polarplot_single":
                    answers.append(str(float(self.cyl_theta.get())))
                elif kw == "polarplot_integrate":
                    answers += [str(float(self.cyl_theta_start.get())),
                                str(float(self.cyl_theta_end.get())),
                                str(int(self.cyl_n_points.get()))]
                else:
                    answers.append(str(int(self.cyl_n_points.get())))
            except ValueError as e:
                messagebox.showerror("Error", f"Invalid value: {e}")
                return

        elif module == "Ellipsoid":
            # ellipsoid_SHS: use_nlwkb, wavelength, n_p, n_m, a, length=2c + angles
            src = SRC_ELLIPSOID
            beta_content = self.txt_beta_ell.get(1.0, tk.END).strip()
            kw = self.ell_keyword.get()
            try:
                answers = [
                    "1",  # NLWKB (n_p and n_m)
                    str(float(self.ell_wavelength.get())),
                    str(float(self.ell_np_real.get())),
                    str(float(self.ell_np_imag.get())),
                    str(float(self.ell_nm_real.get())),
                    str(float(self.ell_nm_imag.get())),
                    str(float(self.ell_radius.get())),
                    str(float(self.ell_length.get())),
                ]
                if kw == "polarplot_single":
                    answers.append(str(float(self.ell_theta.get())))
                elif kw == "polarplot_integrate":
                    answers += [str(float(self.ell_theta_start.get())),
                                str(float(self.ell_theta_end.get())),
                                str(int(self.ell_n_points.get()))]
                else:
                    answers.append(str(int(self.ell_n_points.get())))
            except ValueError as e:
                messagebox.showerror("Error", f"Invalid value: {e}")
                return

        elif module == "Cube":
            # cube_SHS: use_nlwkb, wavelength, n_p, n_m, a, b, c + angles
            src = SRC_CUBE
            beta_content = self.txt_beta_cube.get(1.0, tk.END).strip()
            kw = self.cube_keyword.get()
            try:
                answers = [
                    "1",  # NLWKB
                    str(float(self.cube_wavelength.get())),
                    str(float(self.cube_np_real.get())),
                    str(float(self.cube_np_imag.get())),
                    str(float(self.cube_nm_real.get())),
                    str(float(self.cube_nm_imag.get())),
                    str(float(self.cube_a.get())),
                    str(float(self.cube_b.get())),
                    str(float(self.cube_c.get())),
                ]
                if kw == "polarplot_single":
                    answers.append(str(float(self.cube_theta.get())))
                elif kw == "polarplot_integrate":
                    answers += [str(float(self.cube_theta_start.get())),
                                str(float(self.cube_theta_end.get())),
                                str(int(self.cube_n_points.get()))]
                else:
                    answers.append(str(int(self.cube_n_points.get())))
            except ValueError as e:
                messagebox.showerror("Error", f"Invalid value: {e}")
                return

        else:
            messagebox.showerror("Error", f"Unknown module: {module}")
            return

        if not beta_content:
            messagebox.showerror("Error", "Beta tensor is empty.")
            return
        if not src.exists():
            messagebox.showerror("Error", f"Python module not found:\n{src}")
            return

        _cleanup_old_runs()
        work = tempfile.mkdtemp(prefix="pyshs_", dir=str(WORK_DIR))
        beta_file = Path(work) / "input_beta.txt"
        orient_file = Path(work) / "input_orient.txt"
        out_file = Path(work) / "output.txt"
        with open(beta_file, "w", encoding="utf-8") as f:
            f.write(beta_content + "\n")
        if orient_content is not None:
            with open(orient_file, "w", encoding="utf-8") as f:
                f.write(orient_content + "\n")

        stdin_data = "\n".join(answers) + "\n"
        # All modules are pure Python: python3 script.py keyword [files...]
        # HRS / Sphere / Cylinder / Ellipsoid / Cube / Mie: keyword + beta + output
        # SHS: keyword + beta + orientation + output
        # Modules are run with "python -m py_shs.<module>" so that they work once installed.
        modname = f"py_shs.{src.stem}"
        if module == "SHS":
            cmd = [sys.executable, "-u", "-m", modname, kw, str(beta_file), str(orient_file), str(out_file)]
        else:
            cmd = [sys.executable, "-u", "-m", modname, kw, str(beta_file), str(out_file)]

        # Intensity units / scale, common to all modules (see shs_units.py)
        try:
            i_scale_val = float(self.out_i_scale.get())
        except (ValueError, AttributeError):
            i_scale_val = 1.0
        if not (i_scale_val > 0):
            i_scale_val = 1.0
        units_choice = self.out_units.get() if hasattr(self, "out_units") else UNITS_CHOICES[0]
        units_env = UNITS_ENV.get(units_choice, "")
        self._log(f"Intensity scale factor = {i_scale_val:g}   units = {units_choice}")

        self.job_running = True
        self._current_proc = None
        self._set_run_buttons_state(tk.DISABLED)
        if hasattr(self, "btn_cancel_calc"):
            self.btn_cancel_calc.config(state=tk.NORMAL)
        self._log(f"\n{'='*55}\nStarting {module}\nCommand: {' '.join(cmd)}\n"
                  f"stdin: {answers}\n"
                  f"(live progress – no time limit)\n"
                  f"{'='*55}\n")

        def worker():
            stdout_lines = []
            stderr_lines = []
            try:
                run_env = {**os.environ, "LC_ALL": "C"}
                run_env["PYSHS_I_SCALE"] = str(i_scale_val)
                if units_env:
                    run_env["PYSHS_UNITS"] = units_env
                else:
                    run_env.pop("PYSHS_UNITS", None)
                proc = subprocess.Popen(
                    cmd,
                    stdin=subprocess.PIPE,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    text=True,
                    bufsize=1,          # line-buffered on Python side
                    cwd=str(work),
                    env=run_env,
                )
                self._current_proc = proc

                try:
                    proc.stdin.write(stdin_data)
                    proc.stdin.flush()
                    proc.stdin.close()
                except Exception as e:
                    self.after(0, lambda err=str(e): self._log(f"stdin write error: {err}"))

                def reader(pipe, collector, is_stderr=False):
                    try:
                        while True:
                            line = pipe.readline()
                            if line == "":
                                break
                            collector.append(line)
                            msg = line.rstrip("\n\r")
                            if msg:
                                prefix = "[stderr] " if is_stderr else ""
                                # default-arg to capture current msg value
                                self.after(0, lambda m=prefix + msg: self._log_live(m))
                    finally:
                        try:
                            pipe.close()
                        except Exception:
                            pass

                t_out = threading.Thread(target=reader, args=(proc.stdout, stdout_lines, False), daemon=True)
                t_err = threading.Thread(target=reader, args=(proc.stderr, stderr_lines, True), daemon=True)
                t_out.start()
                t_err.start()

                returncode = proc.wait()
                t_out.join(timeout=10)
                t_err.join(timeout=10)

                stdout = "".join(stdout_lines)
                stderr = "".join(stderr_lines)
                self._current_proc = None
                self.after(0, lambda: self._on_calc_finished(
                    returncode, stdout, stderr, work, out_file, kw, module, orient_content
                ))
            except Exception as e:
                if self._current_proc is not None:
                    try:
                        self._current_proc.kill()
                    except Exception:
                        pass
                    self._current_proc = None
                self.after(0, lambda err=str(e): self._log(f"Exception: {err}"))
                self.after(0, self._reset_run_buttons)

        threading.Thread(target=worker, daemon=True).start()

    def _cancel_calculation(self):
        """Stop the running calculation process (if any)."""
        proc = self._current_proc
        if proc is None:
            self._log("No calculation to stop.")
            return
        try:
            proc.kill()
            self._log(">>> Calculation stopped by user.")
        except Exception as e:
            self._log(f"Could not stop process: {e}")

    def _set_run_buttons_state(self, state):
        for attr in ("btn_run_hrs", "btn_run_shs", "btn_run_sphere", "btn_run_sphere_mie", "btn_run_cylinder", "btn_run_ellipsoid", "btn_run_cube"):
            if hasattr(self, attr):
                getattr(self, attr).config(state=state)

    def _reset_run_buttons(self):
        self.job_running = False
        self._current_proc = None
        self._set_run_buttons_state(tk.NORMAL)
        if hasattr(self, "btn_cancel_calc"):
            self.btn_cancel_calc.config(state=tk.DISABLED)

    def _on_calc_finished(self, returncode, stdout, stderr, work, out_file, keyword, module, orient_content):
        self._reset_run_buttons()
        # stdout was already shown live line by line
        self._log(f"\n--- Process finished (code {returncode}) ---")
        if returncode != 0:
            self._log(f"*** Error code {returncode} ***")
            if stderr and stderr.strip():
                self._log("--- stderr (rappel) ---\n" + stderr)
            return
        if out_file.exists():
            self._log("\n--- Fichier de sortie ---\n" + out_file.read_text(encoding="utf-8", errors="replace"))

        out_plot = Path(work) / "out_plot"
        if out_plot.exists():
            self._show_intensity_plot(out_plot, keyword)
        else:
            self._log("No out_plot file.")

        # Dipole assembly (SHS) or particle geometry (continuum models)
        if module == "SHS" and orient_content:
            self._show_assembly_from_text(orient_content)
        elif module in ("Sphere", "SphereMie", "Cylinder", "Ellipsoid", "Cube"):
            self._show_particle_geometry(module)

        self._log(f"\n*** Calculation {module} finished ***\n")

    # ------------------------------------------------------------------
    # Plots
    # ------------------------------------------------------------------
    def _show_intensity_plot(self, out_plot: Path, keyword: str):
        for w in self.plot_intensity_container.winfo_children():
            w.destroy()
        try:
            data = self._parse_out_plot(out_plot)
        except Exception as e:
            self._log(f"Failed to parse out_plot: {e}")
            return

        fig = Figure(figsize=(8.5, 5.2), dpi=110)
        ax = fig.add_subplot(111)
        if keyword in ("polarplot_single", "polarplot_integrate") and "IV" in data:
            ax.plot(data["gamma"], data["IV"], label="IV", lw=2)
            ax.plot(data["gamma"], data["IH"], label="IH", lw=2)
            ax.set_xlabel("γ (°)")
            ax.set_ylabel("SHS Intensity")
            ax.set_title("Polarization-resolved SHS")
            ax.legend()
            ax.grid(True, alpha=0.3)
        elif "IH_0" in data:
            g = data["gamma"]
            g_full = g + [-x for x in reversed(g[:-1])]
            for key, lab in [("IV_0", "IV(0°)"), ("IV_90", "IV(90°)"),
                             ("IH_0", "IH(0°)"), ("IH_90", "IH(90°)")]:
                vals = data[key]
                ax.plot(g_full, vals + list(reversed(vals[:-1])), label=lab)
            ax.set_xlabel("Angle (°)")
            ax.set_ylabel("SHS Intensity")
            ax.set_title("Angle-resolved SHS")
            ax.legend()
            ax.grid(True, alpha=0.3)
        else:
            ax.text(0.5, 0.5, "Unrecognized data", ha="center", transform=ax.transAxes)

        fig.subplots_adjust(left=0.12, right=0.96, top=0.90, bottom=0.14)
        self._last_intensity_fig = fig
        canvas = FigureCanvasTkAgg(fig, master=self.plot_intensity_container)
        canvas.draw()
        canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True, padx=4, pady=4)
        self.plot_nb.select(self.tab_plot_intensity)

    def _show_assembly_from_text(self, content: str):
        for w in self.plot_assembly_container.winfo_children():
            w.destroy()

        x, y, z, phi, theta = [], [], [], [], []
        for line in content.splitlines():
            parts = line.split()
            if len(parts) < 6:
                continue
            try:
                vals = [float(p) for p in parts[:6]]
            except ValueError:
                continue
            if abs(vals[0]) <= 3.2 and abs(vals[1]) <= 3.2:
                phi.append(vals[0]); theta.append(vals[1])
                x.append(vals[3]); y.append(vals[4]); z.append(vals[5])
            else:
                x.append(vals[0]); y.append(vals[1]); z.append(vals[2])
                phi.append(vals[3]); theta.append(vals[4])

        if not x:
            self._log("No assembly data to plot.")
            return

        x, y, z = np.array(x), np.array(y), np.array(z)
        phi, theta = np.array(phi), np.array(theta)
        # Arrow length proportional to scale
        scale = max(np.ptp(x), np.ptp(y), np.ptp(z), 1.0) * 0.15

        fig = Figure(figsize=(7.5, 6.2), dpi=110)
        ax = fig.add_subplot(111, projection="3d")
        ax.quiver(x, y, z,
                  scale * np.cos(phi) * np.sin(theta),
                  scale * np.sin(phi) * np.sin(theta),
                  scale * np.cos(theta),
                  color="crimson", arrow_length_ratio=0.3)
        ax.set_xlabel("X")
        ax.set_ylabel("Y")
        ax.set_zlabel("Z")
        ax.set_title(f"Dipole assembly (N={len(x)})", pad=12)
        fig.subplots_adjust(left=0.02, right=0.98, top=0.90, bottom=0.05)
        self._last_assembly_fig = fig
        canvas = FigureCanvasTkAgg(fig, master=self.plot_assembly_container)
        canvas.draw()
        canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True, padx=4, pady=4)


    def _show_particle_geometry(self, module: str):
        """Draw the continuum particle shape in the Assembly / Geometry tab."""
        for w in self.plot_assembly_container.winfo_children():
            w.destroy()

        try:
            if module in ("Sphere", "SphereMie"):
                R = float(self.sph_radius.get())
                title = f"Sphere  R = {R:g} nm"
                shape = "sphere"
                dims = (R,)
            elif module == "Cylinder":
                R = float(self.cyl_radius.get())
                L = float(self.cyl_length.get())
                title = f"Cylinder  R = {R:g} nm,  L = {L:g} nm"
                shape = "cylinder"
                dims = (R, L)
            elif module == "Ellipsoid":
                a = float(self.ell_radius.get())
                length = float(self.ell_length.get())
                c = length / 2.0
                title = f"Ellipsoid  a = {a:g} nm,  2c = {length:g} nm"
                shape = "ellipsoid"
                dims = (a, c)
            elif module == "Cube":
                a = float(self.cube_a.get())
                b = float(self.cube_b.get())
                c = float(self.cube_c.get())
                title = f"Box  a={a:g}, b={b:g}, c={c:g} nm"
                shape = "cube"
                dims = (a, b, c)
            else:
                return
        except (ValueError, AttributeError) as e:
            self._log(f"Cannot draw geometry: {e}")
            return

        fig = Figure(figsize=(7.5, 6.2), dpi=110)
        ax = fig.add_subplot(111, projection="3d")

        if shape == "sphere":
            R = dims[0]
            u = np.linspace(0, 2 * np.pi, 40)
            v = np.linspace(0, np.pi, 24)
            x = R * np.outer(np.cos(u), np.sin(v))
            y = R * np.outer(np.sin(u), np.sin(v))
            z = R * np.outer(np.ones_like(u), np.cos(v))
            ax.plot_surface(x, y, z, color="#3b82f6", alpha=0.55, linewidth=0, antialiased=True)
            # equator + meridians for readability
            t = np.linspace(0, 2 * np.pi, 80)
            ax.plot(R * np.cos(t), R * np.sin(t), 0 * t, color="#1e3a8a", lw=1.2)
            ax.plot(R * np.cos(t), 0 * t, R * np.sin(t), color="#1e3a8a", lw=0.8, alpha=0.7)
            ax.plot(0 * t, R * np.cos(t), R * np.sin(t), color="#1e3a8a", lw=0.8, alpha=0.7)
            lim = R * 1.15

        elif shape == "cylinder":
            R, L = dims
            z = np.linspace(-L / 2, L / 2, 30)
            th = np.linspace(0, 2 * np.pi, 40)
            Th, Z = np.meshgrid(th, z)
            X = R * np.cos(Th)
            Y = R * np.sin(Th)
            ax.plot_surface(X, Y, Z, color="#10b981", alpha=0.55, linewidth=0, antialiased=True)
            # top / bottom disks
            for z0 in (-L / 2, L / 2):
                ax.plot(R * np.cos(th), R * np.sin(th), np.full_like(th, z0), color="#065f46", lw=1.2)
            lim = max(R, L / 2) * 1.2

        elif shape == "ellipsoid":
            a, c = dims
            u = np.linspace(0, 2 * np.pi, 40)
            v = np.linspace(0, np.pi, 24)
            x = a * np.outer(np.cos(u), np.sin(v))
            y = a * np.outer(np.sin(u), np.sin(v))
            z = c * np.outer(np.ones_like(u), np.cos(v))
            ax.plot_surface(x, y, z, color="#a855f7", alpha=0.55, linewidth=0, antialiased=True)
            t = np.linspace(0, 2 * np.pi, 80)
            ax.plot(a * np.cos(t), a * np.sin(t), 0 * t, color="#6b21a8", lw=1.2)
            ax.plot(a * np.cos(t), 0 * t, c * np.sin(t), color="#6b21a8", lw=0.8, alpha=0.7)
            lim = max(a, c) * 1.15

        else:  # cube / box
            a, b, c = dims
            # 8 vertices of the box centred at origin
            verts = np.array([
                [-a/2, -b/2, -c/2], [a/2, -b/2, -c/2], [a/2, b/2, -c/2], [-a/2, b/2, -c/2],
                [-a/2, -b/2,  c/2], [a/2, -b/2,  c/2], [a/2, b/2,  c/2], [-a/2, b/2,  c/2],
            ])
            edges = [
                (0, 1), (1, 2), (2, 3), (3, 0),
                (4, 5), (5, 6), (6, 7), (7, 4),
                (0, 4), (1, 5), (2, 6), (3, 7),
            ]
            for i, j in edges:
                ax.plot(*zip(verts[i], verts[j]), color="#ea580c", lw=2)
            # translucent faces
            from mpl_toolkits.mplot3d.art3d import Poly3DCollection
            faces = [
                [verts[i] for i in (0, 1, 2, 3)],
                [verts[i] for i in (4, 5, 6, 7)],
                [verts[i] for i in (0, 1, 5, 4)],
                [verts[i] for i in (2, 3, 7, 6)],
                [verts[i] for i in (1, 2, 6, 5)],
                [verts[i] for i in (0, 3, 7, 4)],
            ]
            ax.add_collection3d(Poly3DCollection(
                faces, facecolors="#fb923c", alpha=0.35, edgecolors="#c2410c", linewidths=0.5
            ))
            lim = max(a, b, c) / 2 * 1.25

        ax.set_xlim(-lim, lim)
        ax.set_ylim(-lim, lim)
        ax.set_zlim(-lim, lim)
        try:
            ax.set_box_aspect((1, 1, 1))
        except Exception:
            pass
        ax.set_xlabel("X (nm)")
        ax.set_ylabel("Y (nm)")
        ax.set_zlabel("Z (nm)")
        ax.set_title(title, pad=14)
        fig.subplots_adjust(left=0.02, right=0.98, top=0.88, bottom=0.05)
        self._last_assembly_fig = fig
        canvas = FigureCanvasTkAgg(fig, master=self.plot_assembly_container)
        canvas.draw()
        canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True, padx=4, pady=4)
        # switch to geometry tab so the user sees it
        try:
            self.plot_nb.select(self.tab_plot_assembly)
        except Exception:
            pass

    def _parse_out_plot(self, path: Path) -> dict:

        lines = path.read_text(encoding="utf-8", errors="replace").strip().splitlines()
        start = 0
        if lines:
            try:
                float(lines[0].split()[0])
            except (ValueError, IndexError):
                start = 1
        gamma, IV, IH = [], [], []
        IH_0, IH_90, IV_0, IV_90 = [], [], [], []
        for line in lines[start:]:
            parts = line.split()
            if len(parts) < 3:
                continue
            try:
                vals = [float(p) for p in parts]
            except ValueError:
                continue
            if len(vals) == 3:
                gamma.append(vals[0])          # unit decided for the whole file below
                IV.append(vals[1]); IH.append(vals[2])
            elif len(vals) >= 5:
                gamma.append(vals[0])
                IH_0.append(vals[1]); IH_90.append(vals[2])
                IV_0.append(vals[3]); IV_90.append(vals[4])
        if IV and gamma and len(gamma) == len(IV) and max(gamma) <= 2 * np.pi + 0.2:
            # polar plots are written with gamma in radians (HRS, SHS, sphere, cylinder,
            # ellipsoid, cube, Mie); decided once for the whole file (it used to be decided
            # row by row, which mixed units in a file)
            gamma = [g * 180.0 / np.pi for g in gamma]
        result = {"gamma": gamma}
        if IV:
            result["IV"] = IV
            result["IH"] = IH
        if IH_0:
            result["IH_0"] = IH_0
            result["IH_90"] = IH_90
            result["IV_0"] = IV_0
            result["IV_90"] = IV_90
        return result


def main():
    """Entry point of the ``pyshs-optics-gui`` command."""
    app = PySHSApp()
    app.mainloop()


if __name__ == "__main__":
    main()

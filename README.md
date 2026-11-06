# Py_SHS – A Python Second Harmonic Scattering package.

[![PyPI version](https://img.shields.io/pypi/v/pyshs-optics.svg)](https://pypi.org/project/pyshs-optics/)
[![Python versions](https://img.shields.io/pypi/pyversions/pyshs-optics.svg)](https://pypi.org/project/pyshs-optics/)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](https://github.com/pmgassin/Py_SHS/blob/Py_SHS_V4/LICENSE)
[![User guide](https://img.shields.io/badge/docs-User%20Manual%20(PDF)-orange.svg)](https://github.com/pmgassin/Py_SHS/blob/Py_SHS_V4/Documentation/PySHS_User_Manual.pdf)

**Py_SHS** is a collection of programs that compute the **second-harmonic intensity** scattered by
a molecule, or by a colloidal / supramolecular structure (sphere, cylinder, ellipsoid, box,
assembly of oriented dipoles). It covers incoherent Hyper-Rayleigh Scattering (HRS) as well as
coherent Second Harmonic Scattering (SHS), and comes with a graphical interface.

Written in Python on top of NumPy: **no C compiler is needed**, and it installs with a single
`pip` command.

> 📖 **User guide:** [PySHS_User_Manual.pdf](https://github.com/pmgassin/Py_SHS/blob/Py_SHS_V4/Documentation/PySHS_User_Manual.pdf)
> 🐛 **Issues and questions:** [github.com/pmgassin/Py_SHS/issues](https://github.com/pmgassin/Py_SHS/issues)

---

## ✨ What it computes

| Module | Description |
|--------|-------------|
| **HRS** | Incoherent Hyper-Rayleigh Scattering (single molecule) |
| **SHS** | Coherent SHS of an assembly of N oriented dipoles |
| **Sphere** | Colloidal sphere (NLRGD / NLWKB continuum models) |
| **Sphere Mie** | Nonlinear Mie theory for spheres |
| **Cylinder** | Continuum model for colloidal cylinders |
| **Ellipsoid** | Spheroid / ellipsoid of revolution |
| **Cube / Box** | Rectangular box (or cube) |

Keywords shared by most modules:

- `polarplot_single` – polar plot at one scattering angle
- `polarplot_integrate` – average over a range of angles
- `angle_scattering` – angular distribution: you enter N, and the code computes N+1 angles from 0° to 180° (step 180°/N), identically in all modules

---

## 📦 Installation

### 1. Requirements

- **Python 3.8 or newer** ([python.org/downloads](https://www.python.org/downloads/))
- **tkinter**, only for the graphical interface. It is included in the official Python installers for
  Windows and macOS, but is a separate system package on some systems (see the note below).

NumPy, SciPy and matplotlib are installed automatically by `pip`.

### 2. Install with pip

```bash
pip install pyshs-optics
```

That's it. Check that the installation worked:

```bash
python -c "import py_shs; print(py_shs.__version__)"
```

> 💡 **Recommended:** install inside a virtual environment, so that Py_SHS and its dependencies stay
> separate from the rest of your Python installation.
>
> ```bash
> python3 -m venv pyshs-env
> source pyshs-env/bin/activate        # Windows: pyshs-env\Scripts\activate
> pip install pyshs-optics
> ```

### 3. Update or uninstall

```bash
pip install --upgrade pyshs-optics
pip uninstall pyshs-optics
```

### tkinter not found?

If launching the GUI fails with `No module named 'tkinter'`, install it for your system:

| System | Command |
|--------|---------|
| Debian / Ubuntu | `sudo apt install python3-tk` |
| Fedora | `sudo dnf install python3-tkinter` |
| macOS (Homebrew Python) | `brew install python-tk` |
| Windows | re-run the Python installer and tick *tcl/tk and IDLE* |

The command-line programs do not need tkinter.

---

## 🚀 Launching the GUI

Once the package is installed, open a terminal and type:

```bash
pyshs-optics-gui
```

or, equivalently:

```bash
python -m py_shs
```

The graphical interface opens. Choose a module, enter your parameters and input files, and run the
calculation: progress is streamed live in the log panel and results are plotted in the window.
An *Intensity units* selector and a scale factor apply to every module.

**Where are my files?** The GUI writes its input / output files and temporary calculation folders in

```text
~/PySHS_Work/
```

(that is `C:\Users\<you>\PySHS_Work` on Windows). You can choose another folder with the environment
variable `PYSHS_WORK`:

```bash
export PYSHS_WORK=/path/to/my/calculations      # macOS / Linux
pyshs-optics-gui
```

If the `pyshs-optics-gui` command is not found, use `python -m py_shs` instead, which always works.

📖 The full description of the interface and of the physical models is in the
**[User Manual (PDF)](https://github.com/pmgassin/Py_SHS/blob/Py_SHS_V4/Documentation/PySHS_User_Manual.pdf)**.

---

## 💻 Command-line use (without the GUI)

Each module can also be run from a terminal. The protocol is: a keyword, the beta (or ki2) file, an
optional orientation file and the output file; the numerical parameters are read from standard input.

| Command | Module |
|---------|--------|
| `pyshs-optics-hrs` | HRS |
| `pyshs-optics-shs` | SHS (dipole assembly) |
| `pyshs-optics-sphere` | Sphere (NLRGD / NLWKB) |
| `pyshs-optics-mie` | Nonlinear Mie theory |
| `pyshs-optics-cylinder` | Cylinder |
| `pyshs-optics-ellipsoid` | Ellipsoid |
| `pyshs-optics-cube` | Cube / box |

Examples (macOS / Linux shell; on Windows use Git Bash or WSL), with a beta file `beta.txt`
in the current folder:

```bash
# HRS polar plot at 90°
printf "800\n1.33\n0.0\n90\n" | pyshs-optics-hrs polarplot_single beta.txt out_hrs.txt

# SHS (N dipoles + orientation file)
printf "100\n800\n1.33\n0.0\n90\n" | pyshs-optics-shs polarplot_single beta.txt orient.txt out_shs.txt

# Sphere NLWKB
printf "1\n800\n1.5\n0.0\n1.33\n0.0\n50\n90\n" | pyshs-optics-sphere polarplot_single beta.txt out_sph.txt

# Cylinder
printf "800\n1.5\n0.0\n1.33\n0.0\n20\n80\n90\n" | pyshs-optics-cylinder polarplot_single beta.txt out_cyl.txt
```

Every module can also be started with `python -m py_shs.<module>` (for example `python -m py_shs.HRS`).

---

## 📄 Input file formats

**Beta / ki(2) file** – one header line, then 27 lines `name value` (order `xxx` … `zzz`):

```text
# surface hyperpolarizability
xxx 0.0
xxy 0.0
...
zzz 1.0
```

**Orientation file** (SHS only) – one line per dipole:

```text
phi  theta  psi  x  y  z
```

(angles in radians, positions in nm).

---

## ⚙️ Environment variables

| Variable | Effect |
|----------|--------|
| `PYSHS_WORK=dir` | Working folder of the GUI (default: `~/PySHS_Work`) |
| `PYSHS_THREADS=k` | Thread count for HRS / SHS (default: CPU cores) |
| `PYSHS_LEGACY=1` | Reproduce two old C quirks of HRS / SHS (affects only bv1, bv2, bh1, bh2, which are zero after the isotropic average) |
| `PYSHS_I_SCALE=x` | Intensity scale factor, applied by **every** module (HRS, SHS, sphere, cylinder, ellipsoid, cube, Mie) |
| `PYSHS_UNITS=rad\|raw` | Intensity units common to all modules (see *Intensity units* below) |
| `PYSHS_ANGLE_RANGE=half\|negative\|full` | Angle range of `angle_scattering` (see *Scattering-angle grid* below) |
| `PYSHS_GRID=fixed` | Use the historical fixed quadrature grids (see *Automatic quadrature grids*) |
| `PYSHS_GRID_SCALE=x` | Multiplies the size parameter used by the grid rules (2 = safer and slower, 0.5 = faster) |
| `PYSHS_GRID_MAX=n` | Upper bound on the grid size (memory / time guard) |
| `PYSHS_THETA_RULE=midpoint` | Former mid-point rule in theta |

---

## 🔬 Technical notes

The notes below describe the numerical conventions of this version. Click a title to expand it.

<details>
<summary><b>Package layout</b></summary>

```text
py_shs/                      (installed in site-packages)
├── pyshs_gui.py            # Tkinter GUI  -> command `pyshs-optics-gui`
├── HRS.py, SHS.py
├── sphere_SHS.py, cylinder_SHS.py, ellipsoid_SHS.py, cube_SHS.py
├── shs_mie.py              # nonlinear Mie theory
├── graph_*.py, script_*.py # plots and dipole generators (sphere, spheroid, cylinder)
├── shs_grid.py, shs_units.py, shs_angles.py, shs_io.py   # shared helpers
└── cli.py                  # console commands pyshs-optics-hrs, pyshs-optics-shs, ...
```

All calculation modules are pure Python and are invoked with the same interpreter
(`python -m py_shs.<module>`).

</details>

<details>
<summary><b>What changed in V4 (pure Python)</b></summary>

- All former `.c` sources (`HRS.c`, `SHS.c`, `sphere_SHS.c`, `cylinder_SHS.c`, `ellipsoid_SHS.c`,
  `cube_SHS.c`) are replaced by equivalent NumPy implementations with the **same physics and the same
  stdin / file protocol**.
- The GUI no longer compiles anything; it launches `python -u -m py_shs.<module> …`.
- Dependencies: Python 3 + NumPy (+ SciPy for `shs_mie.py`; tkinter / matplotlib for the GUI).

Equivalence with the previous C programs: the C sources/binaries are **not** part of this package, so
this package does not itself demonstrate it. What *is* verified (see the `tests/` folder of the GitHub
repository): analytic values, an independent brute-force evaluation of the polarisation formulas, mutual
consistency of all the models (absolute scale included) and a frozen snapshot of the current outputs. To
check against the C programs, run both on the same input and compare the output files.

`PYSHS_LEGACY=1` (HRS and SHS only) reproduces two defects of the old C code (normalisation of bv1, bv2,
bh1, bh2, and the pairing of the sin*cos cross terms of bh1, bh2). These four coefficients are identically
zero after the isotropic orientational average (shown analytically and tested numerically), so the option
changes only rounding noise in them and none of av, bv, cv, ah, bh, ch, I2 and I4.

</details>

<details>
<summary><b>Hyperpolarizability tensor (symmetrisation)</b></summary>

All modules assume beta_ijk = beta_ikj (intrinsic permutation symmetry in the last two indices). The
beta/surface-susceptibility file reader of HRS, SHS, sphere_SHS, cylinder_SHS, ellipsoid_SHS and cube_SHS
symmetrises the input (beta_ijk <- (beta_ijk + beta_ikj)/2) and prints a WARNING if the file was not
symmetric. `shs_mie.py` already averages the equivalent components (xzx/xxz, yzy/yyz).

</details>

<details>
<summary><b>Orientation quadrature (normalisation)</b></summary>

The orientation average over the polar angle theta uses a Gauss-Legendre rule in cos(theta) in HRS, SHS,
cylinder_SHS and ellipsoid_SHS (cube_SHS already did). The weights then sum exactly to 8*pi^2 (SHS/HRS) so
the 1/(8*pi^2) normalisation is exact, and all modules normalise consistently. The former midpoint rule had
a weight sum of 1.0041 * 8*pi^2, i.e. a +0.4 % bias on rotation-invariant quantities. Set
`PYSHS_THETA_RULE=midpoint` to restore the previous grid (e.g. to reproduce earlier results). Grid sizes
(N_PHI, N_THETA, ...) are unchanged.

</details>

<details>
<summary><b>Automatic quadrature grids (convergence)</b></summary>

All intensities contain the phase factor exp(i q.r) with |q| <= q_max = 2 Re(k_w) + Re(K_2w). A grid that is
fine for a small object silently becomes under-resolved for a large one, so the grids are sized from the
size parameter `a = q_max x (largest half-extent)` by `py_shs/shs_grid.py` (rules calibrated numerically
against very fine reference grids; typical error <= 1e-4, usually ~1e-7):

| Module | What is sized | Rule |
|---|---|---|
| `sphere_SHS` | surface grid (Gauss-Legendre in cos theta_c x trapezoid in phi_c) | ny >= a/2 + 20, nx = 2 ny |
| `cylinder_SHS` | azimuth points, orientation grid; axial integral is analytic (sinc) | n_beta >= qR + 14 ; N >= q sqrt(R^2+L^2/4) + 12 |
| `ellipsoid_SHS` | surface (alpha, beta) and orientation grids | n_alpha >= a/2+20, n_beta >= a+16, N >= a+12 |
| `cube_SHS` | orientation grid (faces are analytic) | N >= q d + 12, d = half diagonal |
| `SHS` | orientation grid (phi, theta, psi) | N_phi = N_psi >= 2 a + 8, N_theta = N/2, a = q max&#124;r&#124; |
| `HRS`, `shs_mie` | nothing to size (HRS: single molecule; Mie: Wiscombe cutoff already adaptive) | |

The grid used is printed in the run log (`Grid: ...` / `Orientation grid: ...`). The historical values (the
module constants) are used as minimum sizes: the grids are never coarser than before (for sphere, dipole
assemblies and box, small objects use exactly the historical sizes).

* `PYSHS_GRID=fixed` : use the historical fixed grids (reproduces earlier runs; may be inaccurate for large objects).
* `PYSHS_GRID_SCALE=<float>` : multiplies the size parameter used by the rules (2 = safer and slower, 0.5 = faster).
* `PYSHS_GRID_MAX=<int>` : upper bound on the grid size (memory / time guard). When the bound prevents the rule from being satisfied, a `WARNING` is printed.
* `PYSHS_THETA_RULE=midpoint` : former mid-point rule in theta.

Other changes: `cylinder_SHS` integrates the axial direction analytically (the former `N_Z` no longer
exists, so no aliasing for long rods); `ellipsoid_SHS`, `cube_SHS` and `cylinder_SHS` are vectorised
(5-10x faster).

</details>

<details>
<summary><b>Intensity units (comparing modules)</b></summary>

Two conventions exist, related by the radiative prefactor `P = (|K_2w|^2 / 4 pi)^2`
(`K = 4 pi n / lambda`, in nm^-1):

| Convention | Meaning | Default for |
|---|---|---|
| `rad` | `I = r^2 \|E_2w\|^2 / \|E_w\|^4`, prefactor `P` **included** | sphere, cylinder, ellipsoid, cube, `shs_mie` |
| `raw` | orientation average of `\|sum beta exp(i q.r)\|^2`, **no** prefactor (beta in atomic units; e.g. `<beta_ZZZ^2> = 1/7` for `beta_zzz = 1`) | HRS, SHS |

`I_rad = P * I_raw`. `PYSHS_UNITS=rad` or `raw` forces one convention on every module, and `PYSHS_I_SCALE`
multiplies every module's result, so any two modules can be put on the same scale. The convention used is
printed in the run log and in the header of the output file. Default outputs are unchanged. (HRS/SHS switch
to scientific notation when the multiplier is not 1, since `rad` values are ~1e-9.) The GUI has an
*Intensity units* selector and a scale factor that apply to all modules.

Absolute consistency (checked numerically): a continuum model with surface susceptibility `chi` is
equivalent to a dipole assembly of surface density `sigma` [dipoles/nm^2] and molecular `beta` with
`chi = sigma * beta` (so chi is in units of beta per nm^2). With both on the same scale,
`I_SHS(assembly) = sigma^2 * I_continuum` (agreement to 1e-5 for a sphere, cylinder and cube), and the
nonlinear Mie theory equals the sphere NLRGD model at matched index.

</details>

<details>
<summary><b>Scattering-angle grid (<code>angle_scattering</code>)</b></summary>

All modules use the same grid: the number N you enter is the number of intervals between 0 and 180
degrees, so N+1 angles `Theta = 0, 180/N, ..., 180`. Previously HRS, SHS and the cylinder covered
[-180, 0] with N+1 points, while the sphere, ellipsoid, cube and Mie covered [-180, +180] with 2N+1 points,
so the same N did not mean the same thing (and half of the sphere/ellipsoid/cube/Mie work was redundant).

A half range is sufficient because the four reported intensities satisfy `I(Theta) = I(-Theta)`; this was
checked numerically in every module with a random chiral tensor.

`PYSHS_ANGLE_RANGE` restores the former ranges: `negative` (-180...0, N+1 points; former HRS/SHS/cylinder
output) or `full` (-180...+180, 2N+1 points; former sphere/ellipsoid/cube/Mie output). The values at a given
angle are identical in all modes. `graph_angle.py` accepts any of the three ranges and plots the full
-180...+180 curve.

</details>

<details>
<summary><b>Models: notes and limitations</b></summary>

**Output files.** `out_plot` and `out_arrows_gnuplot` are written in the directory of the output file given
on the command line. Running twice in the same directory overwrites them.

**Models.**

- `cylinder_SHS` models the lateral surface only (no end caps): use it for objects long compared with their
  radius. The surface frame uses the exact polar angle pi/2 (the old C code used the truncated 1.570796;
  `PYSHS_THETA_MESO_LEGACY=1` restores it, effect ~1e-7).
- The classic Rayleigh-Gans sphere model is the NLWKB model with `n_p = n_m`: in the GUI enter the same
  index twice.
- `shs_mie.py` has three independent surface constants (isotropic achiral surface: `zzz`, `zxx`/`zyy`,
  `xzx`/`xxz`/`yzy`/`yyz`); other tensor components are ignored and a warning lists them. The
  command-line/GUI interface uses the same index at omega and 2 omega and drops the imaginary part of the
  medium index (`shs_sphere()` accepts distinct indices).
- With an absorbing medium (complex index) the phase factor is `exp(i q.r)` with complex `q`, so the result
  depends on the origin of the positions (the attenuation is `exp(-Im q . r)`): put the origin where the
  incident beam enters the object.
- `polarplot_integrate` is a plain average over the scattering angles, not weighted by `sin(Theta)`.
- The sphere and Mie outputs do not contain `b1`, `b2` (they are zero anyway, see the legacy note).

**Generators and converters (`py_shs/`).**

- `script_sphere.py`, `script_sphere_random.py`, `script_spheroid_random.py`, `script_cylinder_random.py`:
  write `phi theta psi x y z` files for `SHS.py` (dipoles along the outward normal; disorder = Gaussian
  noise on `phi`, `theta`; a noise on `phi` is a real tilt of `sigma*sin(theta)`).
  `script_spheroid_random.py` (a, c, N, sigma) places the points with a uniform *surface* density and
  reproduces `ellipsoid_SHS.py` through `chi = sigma * beta`.
- `script_gauss2PySHS.py job.log beta_out [--marker TEXT] [--occurrence K] [--offset N]` reads the first
  `Beta(-2w;w,w)` block (see its docstring; the layout was not checked against every Gaussian version) and
  writes a tensor symmetric in the last two indices.
- `graph_assembly.py orientation_file image.png [--positions-first]` expects `phi theta psi x y z`.

**GUI.** The *Intensity units* selector and the scale factor apply to every module; run directories
`~/PySHS_Work/pyshs_*` older than 7 days are removed automatically; the unit of polar plots (radians) is
decided once per file.

</details>

<details>
<summary><b>Tests (GitHub repository)</b></summary>

A regression and validation suite (standard library only, no pytest needed) is kept in the `tests/` folder
of the [GitHub repository](https://github.com/pmgassin/Py_SHS):

```bash
python tests/run_tests.py        # everything, about 2 minutes
python tests/run_tests.py -f     # fast subset (skips the slow convergence tests)
```

| File | What it checks |
|---|---|
| `test_core.py` | Euler matrices, exact orientation quadrature, analytic HRS values, HRS = SHS(N=1), the polarisation expansion against a brute-force evaluation, legacy mode, Theta -> -Theta symmetry |
| `test_models.py` | absolute scale of sphere, ellipsoid, cylinder, cube and Mie against dipole assemblies (chi = sigma beta), symmetrisation of beta in every module |
| `test_grids.py` | automatic grids against much finer references for large objects |
| `test_units_angles.py` | `PYSHS_UNITS`, `PYSHS_I_SCALE`, the angle grid and its compatibility ranges, `graph_angle.py` |
| `test_golden.py` | frozen snapshot of command-line outputs of every module (it detects unintended changes; it is not a C reference) |
| `test_tools.py` | `compare_outputs.py` |

The Mie tests are skipped when SciPy is not installed.

</details>

---

## 📚 How to cite

If you use Py_SHS in your work, please cite:

- J. Chem. Inf. Model. **2020**, 60, 12, 5912–5917
- J. Chem. Inf. Model. **2021**, 61, 11, 5719

## 👥 Authors and license

Dr Pierre-Marie GASSIN and Dr Gaelle GASSIN – ICGM / ENSCM.

Open-source research software, released under the MIT License (see the
[LICENSE](https://github.com/pmgassin/Py_SHS/blob/Py_SHS_V4/LICENSE) file).

## 📘 Documentation note

The PDF [User Manual] is available here: (https://github.com/pmgassin/Py_SHS/blob/Py_SHS_V4/Documentation/PySHS_User_Manual.pdf)

SHS_simulation is a collection of different programs which can compute the Second Harmonic Intensity Scattered by a supramolecular structure. 
Don't forget to cite this work:
J. Chem. Inf. Model. 2020, 60, 12, 5912–5917 J. Chem. Inf. Model. 2021, 61, 11, 5719


# PySHS Graphical Interface


Desktop GUI for the **PySHS** Second Harmonic Scattering codes:

| Module | Description |
|--------|-------------|
| **HRS** | Incoherent Hyper-Rayleigh Scattering |
| **SHS** | Coherent SHS of supramolecular aggregates |
| **Sphere Object** | Continuum model for colloidal spheres |
| **Cylinder Object** | Continuum model for colloidal cylinders |

---

## Directory layout

```text
PySHS_GUI/
├── run_gui.py              # Launch the GUI (recommended)
├── README.md
├── Src/                    # All source code
│   ├── pyshs_gui.py
│   ├── HRS.c, SHS.c, sphere_SHS.c, cylinder_SHS.c
│   ├── graph_*.py, script_*.py
│   └── script_gauss2PySHS.py
├── bin/                    # Compiled C executables (created automatically)
│   ├── HRS, SHS, sphere_SHS, cylinder_SHS   # Linux / macOS
│   └── *.exe                                # Windows
├── Work/                   # Input / output of calculations only
│   └── pyshs_*/ , out_plot, temporary files
└── Documentation/          # User manuals (pdf, docx, md, tex)
```

---

## Requirements (all platforms)

| Component | Purpose |
|-----------|---------|
| **Python 3.8+** | GUI (`tkinter`, `matplotlib`, `numpy`) |
| **C compiler** (`gcc` or `clang`) | Compile the C modules |

**Automatic compilation:** when you start the GUI (`run_gui.py` or `Src/pyshs_gui.py`), it compiles (or updates) the C sources into `bin/` if a compiler is available. You do **not** need a separate compile script on Linux or macOS.

---

## Installation on Linux

### 1. System packages

**Debian / Ubuntu:**

```bash
sudo apt update
sudo apt install python3 python3-tk python3-matplotlib python3-numpy build-essential
```

**Fedora:**

```bash
sudo dnf install python3 python3-tkinter python3-matplotlib python3-numpy gcc
```

### 2. Get the code

```bash
git clone https://github.com/pmgassin/PySHS_GUI.git
cd PySHS_GUI
```

(Or extract a zip archive of the project.)

### 3. Run the GUI

```bash
python3 run_gui.py
# or
python3 Src/pyshs_gui.py
```

The C modules are compiled into `bin/` automatically at startup (or the first time a calculation is launched).

**Manual compile (optional, e.g. for command-line use without GUI):**

```bash
gcc -o bin/HRS Src/HRS.c -lm
gcc -o bin/SHS Src/SHS.c -lm
gcc -o bin/sphere_SHS Src/sphere_SHS.c -lm
gcc -o bin/cylinder_SHS Src/cylinder_SHS.c -lm
```

---

## Installation on macOS

### 1. Tools

Install **Xcode Command Line Tools** (provides `clang`):

```bash
xcode-select --install
```

Install **Python 3** if needed (from [python.org](https://www.python.org/) or Homebrew):

```bash
brew install python
```

Install Python packages:

```bash
pip3 install matplotlib numpy
# tkinter is usually included with the official Python installer;
# with Homebrew you may need: brew install python-tk
```

### 2. Get the code

```bash
git clone https://github.com/pmgassin/SHS_simulation.git
cd PySHS_GUI
```

### 3. Run the GUI

```bash
python3 run_gui.py
# or
python3 Src/pyshs_gui.py
```

Compilation into `bin/` is handled by the GUI.

**Manual compile (optional):**

```bash
clang -o bin/HRS Src/HRS.c -lm
clang -o bin/SHS Src/SHS.c -lm
clang -o bin/sphere_SHS Src/sphere_SHS.c -lm
clang -o bin/cylinder_SHS Src/cylinder_SHS.c -lm
```


---

## Installation on Windows

### 1. Install Python

1. Download Python 3 from [https://www.python.org/downloads/](https://www.python.org/downloads/)
2. During setup, enable **“Add python.exe to PATH”**
3. Open a new Command Prompt and check:

```bat
python --version
```

Install scientific packages:

```bat
python -m pip install matplotlib numpy
```

(`tkinter` is included with the official Windows installer.)

### 2. Install a C compiler

You need **gcc** on the PATH. Recommended options:

- **MSYS2**: [https://www.msys2.org/](https://www.msys2.org/)  
  then: `pacman -S mingw-w64-ucrt-x86_64-gcc` and add `C:\msys64\ucrt64\bin` to PATH  
- **WinLibs** (standalone MinGW): [https://winlibs.com/](https://winlibs.com/)

Check:

```bat
gcc --version
```

If you see errors like `clock_gettime64` / `cc1.exe`, your MinGW PATH is mixed (see Troubleshooting).

### 3. Get the code

Clone with Git, or download the ZIP from GitHub and extract it.

### 4. Run the GUI

```bat
cd path\to\PySHS_GUI
python run_gui.py
```

or:

```bat
python Src\pyshs_gui.py
```

The GUI compiles the C sources into `bin\` when needed.


**Manual compile (optional):**

```bat
gcc -o bin\HRS.exe Src\HRS.c -lm
gcc -o bin\SHS.exe Src\SHS.c -lm
gcc -o bin\sphere_SHS.exe Src\sphere_SHS.c -lm
gcc -o bin\cylinder_SHS.exe Src\cylinder_SHS.c -lm
```


---

## Quick reference

| Platform | Launch GUI |
|----------|------------|
| **Linux / macOS** | `python3 run_gui.py` |
| **Windows** | `python run_gui.py` |

| Platform | Optional helper |
|----------|-----------------|
| **Windows** | `python run_gui.py` |

Binaries are written to **`bin/`**. Calculation inputs/outputs go to **`Work/`**.

---

## Command-line usage (without GUI)

After the binaries exist in `bin/`:

**Linux / macOS:**

```bash
./bin/HRS  <keyword> <input_beta> <output>
./bin/SHS  <keyword> <input_beta> <input_orient> <output>
./bin/sphere_SHS   <keyword> <input_ki2> <output>
./bin/cylinder_SHS <keyword> <input_ki2> <output>
```

**Windows:**

```bat
bin\HRS.exe  <keyword> <input_beta> <output>
bin\SHS.exe  <keyword> <input_beta> <input_orient> <output>
bin\sphere_SHS.exe   <keyword> <input_ki2> <output>
bin\cylinder_SHS.exe <keyword> <input_ki2> <output>
```

**Keywords:** `polarplot_single` | `polarplot_integrate` | `angle_scattering`

Runtime files such as `out_plot` are written in the **current working directory**.  
When using the GUI, that directory is a temporary folder under `Work/`.

---

## Documentation

User manuals are in the `Documentation/` folder:

- `PySHS_GUI_User_Manual.pdf`
- `PySHS_GUI_User_Manual.docx`
- and markdown / LaTeX sources

---

## Troubleshooting

| Problem | What to try |
|---------|-------------|
| `No module named tkinter` | Install `python3-tk` (Linux) or reinstall Python from python.org (Windows/macOS) |
| `No C compiler found` | Install `build-essential` / Xcode CLT / MinGW and reopen the terminal |
| Windows: `clock_gettime64` / `cc1.exe` | MSYS2 PATH mix-up — use only `ucrt64\bin`, or WinLibs; update with `pacman -Syu` |
| Calculation fails / wrong results | Ensure **Number of dipoles** matches the orientation file (SHS); intensities scale as N² |
| GUI opens but no progress text | Normal for some modes; `angle_scattering` is usually more verbose |

---

## Citation

Please cite:

> J. Chem. Inf. Model. 2020, 60, 12, 5912–5917

when publishing results obtained with PySHS (command line or GUI).

PySHS was developed at Institut Charles Gerhardt Montpellier – ENSCM  
(Support: ANR-15-CE21-0002, CAMOMILS).

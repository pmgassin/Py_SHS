#!/usr/bin/env python3
"""Renomme le module Python 'pyshs' en 'py_shs' dans tout le depot.

A lancer UNE FOIS, depuis la racine du depot (le dossier qui contient pyproject.toml) :
    python3 rename_to_py_shs.py
Les noms de commandes (pyshs-gui, pyshs-hrs ...), le fichier pyshs_gui.py et les
variables PYSHS_WORK / ~/PySHS_Work ne changent pas.
"""
import os
import re
import subprocess
import sys

OLD, NEW = "src/pyshs", "src/py_shs"

if not os.path.isfile("pyproject.toml"):
    sys.exit("Erreur : lancez ce script depuis la racine du depot (pyproject.toml introuvable).")

# 1. Renommer le dossier (git mv si possible)
if os.path.isdir(OLD):
    if os.path.isdir(".git"):
        subprocess.check_call(["git", "mv", OLD, NEW])
    else:
        os.rename(OLD, NEW)
    print(f"Dossier renomme : {OLD} -> {NEW}")
elif os.path.isdir(NEW):
    print(f"{NEW} existe deja, on continue avec les remplacements.")
else:
    sys.exit(f"Erreur : {OLD} introuvable.")

# 2. Remplacements dans les fichiers texte
REPL = [
    ("from pyshs import", "from py_shs import"),
    ("from pyshs.", "from py_shs."),
    ("import pyshs", "import py_shs"),
    ('f"pyshs.{', 'f"py_shs.{'),
    ('= "pyshs.', '= "py_shs.'),                 # points d'entree dans pyproject.toml
    ("-m pyshs", "-m py_shs"),
    ("assert pyshs.", "assert py_shs."),
    ("pip install pyshs", "pip install py-shs"),
    ('name = "pyshs"', 'name = "py_shs"'),
    ("python -m pyshs", "python -m py_shs"),
    ("\npyshs/", "\npy_shs/"),                  # arborescence dans le README
    ("`pyshs/", "`py_shs/"),
    ("#   pyshs/", "#   py_shs/"),
]
VERSION_FROM, VERSION_TO = "4.0.0", "4.0.1"

changed = []
for root, dirs, files in os.walk("."):
    dirs[:] = [d for d in dirs if d not in (".git", "dist", "build", "__pycache__") and not d.endswith(".egg-info")]
    for f in files:
        if not f.endswith((".py", ".toml", ".md", ".yml", ".yaml")) or f == "rename_to_py_shs.py":
            continue
        p = os.path.join(root, f)
        s = open(p, encoding="utf-8").read()
        t = s
        for a, b in REPL:
            t = t.replace(a, b)
        if f == "pyproject.toml":
            t = t.replace(f'version = "{VERSION_FROM}"', f'version = "{VERSION_TO}"')
        if f == "__init__.py":
            t = t.replace(f'__version__ = "{VERSION_FROM}"', f'__version__ = "{VERSION_TO}"')
        if t != s:
            open(p, "w", encoding="utf-8").write(t)
            changed.append(p)

print(f"{len(changed)} fichier(s) modifie(s).")

# 3. Controle : occurrences restantes de l'ancien nom de module
pat = re.compile(r"(?<![\w-])pyshs(?![\w-]|_gui\.py)")
left = []
for root, dirs, files in os.walk("."):
    dirs[:] = [d for d in dirs if d not in (".git", "dist", "build", "__pycache__") and not d.endswith(".egg-info")]
    for f in files:
        if f.endswith((".py", ".toml", ".md", ".yml")) and f != "rename_to_py_shs.py":
            p = os.path.join(root, f)
            for i, line in enumerate(open(p, encoding="utf-8"), 1):
                if pat.search(line):
                    left.append(f"{p}:{i}: {line.strip()[:110]}")
print("\nOccurrences restantes de 'pyshs' (a relire, certaines sont normales) :")
print("\n".join(left) if left else "  aucune")

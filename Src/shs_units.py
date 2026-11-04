"""Intensity units shared by all the PySHS modules.

Two conventions exist, related by the radiative prefactor  P = (|K_2w|^2 / 4 pi)^2 :

  rad : I = r^2 |E_2w|^2 / |E_w|^4.  Radiative prefactor included.  This is the absolute
        scale of the continuum models (sphere, cylinder, ellipsoid, cube) and of the
        nonlinear Mie theory.  The surface susceptibility chi(2) is expressed in units of
        "beta per nm^2" (a dipole assembly of surface density sigma [dipoles/nm^2] and
        molecular hyperpolarizability beta is equivalent to chi = sigma * beta).
  raw : I = <|sum beta exp(i q.r)|^2> (the square of the coherent sum, orientation
        averaged), WITHOUT the prefactor.  This is the natural scale of HRS and SHS, where
        beta is in atomic units  (e.g. <beta_ZZZ^2> = 1/7 for a single beta_zzz = 1).

  I_rad = P * I_raw        (K = 4 pi n / lambda, in nm^-1)

Defaults (historical): HRS and SHS -> raw ; sphere/cylinder/ellipsoid/cube/Mie -> rad.

Environment variables
---------------------
PYSHS_UNITS     rad | raw   force one convention for every module (so that all the
                modules can be compared on the same scale).  Default: historical.
PYSHS_I_SCALE   positive float (default 1): extra multiplier, applied by every module.
"""
import math
import os


class Settings:
    def __init__(self, mode, pref, i_scale, forced):
        self.mode = mode              # "rad" | "raw"
        self.pref = pref              # |K|^2 / (4 pi)
        self.i_scale = i_scale
        self.forced = forced          # True if chosen through PYSHS_UNITS
        self.multiplier = i_scale * (pref * pref if mode == "rad" else 1.0)

    def lines(self):
        """Human readable description (list of strings)."""
        P = self.pref * self.pref
        if self.mode == "rad":
            first = ("Intensity units: r^2|E_2w|^2/|E_w|^4, radiative prefactor "
                     "(|K|^2/4pi)^2 = %.6e included" % P)
        else:
            first = ("Intensity units: raw coherent-sum average |sum beta exp(iq.r)|^2, "
                     "radiative prefactor NOT included (multiply by (|K|^2/4pi)^2 = %.6e "
                     "for r^2|E|^2)" % P)
        src = "PYSHS_UNITS" if self.forced else "module default"
        return [first,
                "Units convention: %s (%s) ; PYSHS_I_SCALE = %.6e ; total multiplier = %.6e"
                % (self.mode, src, self.i_scale, self.multiplier)]


def _i_scale():
    s = os.environ.get("PYSHS_I_SCALE", "").strip()
    if s:
        try:
            v = float(s)
            if v > 0.0 and math.isfinite(v):
                return v
        except ValueError:
            pass
    return 1.0


def settings(grandk, default):
    """grandk: wave vector K_2w (complex or real, nm^-1); default: "rad" or "raw"."""
    k2 = abs(complex(grandk)) ** 2
    pref = k2 / (4.0 * math.pi)
    env = os.environ.get("PYSHS_UNITS", "").strip().lower()
    if env in ("rad", "raw"):
        return Settings(env, pref, _i_scale(), True)
    return Settings(default, pref, _i_scale(), False)


def write_header_lines(fout, unit):
    for line in unit.lines():
        fout.write("  %s\n" % line)

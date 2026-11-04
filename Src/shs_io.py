"""Location of the side files (out_plot, out_arrows_gnuplot) written next to the output file.

They used to be written in the *current directory*, so running a module from another
directory silently created or overwrote files there.  They are now written in the directory
of the output file given on the command line (as shs_mie.py always did); the GUI runs every
job in its own work directory that also holds the output file, so nothing changes for it.
"""
import os

_dir = "."


def set_output(output_file):
    """Call once, as soon as the output file name is known."""
    global _dir
    _dir = os.path.dirname(os.path.abspath(output_file)) or "."


def path(name):
    return os.path.join(_dir, name)

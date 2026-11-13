"""Which file of a run directory holds a quantity.

PHARE writes one file per quantity, named after the diagnostic (``EM_B.h5``,
``ions_pop_protons_density.vtkhdf``, ...). A quantity is asked for by

* master's short key: ``"B"``, ``"E"``, ``"Vi"``, ``"CF"``, ``"Ni"``, ``"rho"``, ``"nppc"``,
  ``"mhd_rho"`` ... (``pharesee/run/run.py:27-40``, the keys of ``Run.all_times()``),
* the file stem itself: ``"EM_B"``, ``"ions_charge_density"``,
  ``"ions_pop_protons_density"`` (per-population files have no short key on master),
* ``pop_stem(pop, what)`` for per-population quantities: ``pop_stem("protons", "density")``.

Nothing is opened here: only the directory listing is read.
"""

from __future__ import annotations

import os

#: master's short keys (run.py:27-40, quantities_per_file) and the "tags" file
SHORT_TO_STEM = {
    "B": "EM_B",
    "E": "EM_E",
    "Vi": "ions_bulkVelocity",
    "CF": "ions_charge_flux",
    "Ni": "ions_charge_density",
    "rho": "ions_mass_density",
    "nppc": "particle_count",
    "mhd_rho": "mhd_rho",
    "mhd_V": "mhd_V",
    "mhd_P": "mhd_P",
    "mhd_rhoV": "mhd_rhoV",
    "mhd_Etot": "mhd_Etot",
    "tags": "tags",
}
STEM_TO_SHORT = {stem: short for short, stem in SHORT_TO_STEM.items()}

#: file extension of each format name
EXTENSIONS = {"h5": ".h5", "vtkhdf": ".vtkhdf"}
_FORMAT_ALIASES = {
    "h5": "h5", "phareh5": "h5", ".h5": "h5",
    "vtkhdf": "vtkhdf", "pharevtkhdf": "vtkhdf", "vtk": "vtkhdf", ".vtkhdf": "vtkhdf",
}

#: per-population file kinds (pharein diagnostics, fromh5.pop_attribute_from_h5)
POP_KINDS = ("density", "charge_density", "flux", "momentum_tensor", "domain", "levelGhost")


def normalize_format(fmt):
    """``"h5"`` / ``"vtkhdf"`` for the accepted spellings, None for None.

    Accepted: "h5", "phareh5", "vtkhdf", "pharevtkhdf", "vtk". Raises ValueError otherwise.
    """
    if fmt is None:
        return None
    key = str(fmt).strip().lower()
    if key not in _FORMAT_ALIASES:
        raise ValueError(f"format={fmt!r}: use 'h5' (phareh5 files) or 'vtkhdf'")
    return _FORMAT_ALIASES[key]


def pop_stem(pop, what):
    """File stem of a per-population quantity, e.g. ``pop_stem("protons", "flux")``."""
    return f"ions_pop_{pop}_{what}"


def quantity_stem(key):
    """File stem for a quantity key (short key or stem). Raises TypeError for a non-string."""
    if not isinstance(key, str) or not key:
        raise TypeError(f"a quantity is a non-empty string such as 'B' or 'EM_B', got {key!r}")
    return SHORT_TO_STEM.get(key, key)


def short_name(stem):
    """Master's key of a file stem (``"EM_B"`` -> ``"B"``), else the stem itself."""
    return STEM_TO_SHORT.get(stem, stem)


def scan(directory):
    """``{stem: {format: path}}`` of every ``.h5`` / ``.vtkhdf`` file in ``directory``.

    Reads the directory listing only. Raises FileNotFoundError when ``directory`` does not
    exist or is not a directory.
    """
    directory = os.fspath(directory)
    if not os.path.isdir(directory):
        raise FileNotFoundError(f"{directory}: no such directory (give the diagnostics folder)")
    found = {}
    for name in sorted(os.listdir(directory)):
        for fmt, ext in EXTENSIONS.items():
            if name.endswith(ext) and len(name) > len(ext):
                stem = name[: -len(ext)]
                found.setdefault(stem, {})[fmt] = os.path.join(directory, name)
    return found


def describe(found):
    """Short text listing the quantities of a ``scan`` result, for error messages."""
    if not found:
        return "no .h5 or .vtkhdf file"
    return ", ".join(f"{stem} ({'/'.join(sorted(fmts))})" for stem, fmts in sorted(found.items()))


def find_quantity_file(directory, key, format=None):
    """Path of the file holding quantity ``key`` in ``directory``.

    key    : master's short key ("B", "Ni", ...) or a file stem ("EM_B",
             "ions_pop_protons_density")
    format : None, "h5" or "vtkhdf" (spellings of ``normalize_format``)
    Raises
      ValueError        both ``<stem>.h5`` and ``<stem>.vtkhdf`` exist and ``format`` is None
                        (the message asks for ``format=``)
      FileNotFoundError no file for that quantity (in that format): the message lists what
                        the directory holds
    The file is not opened: a vtkhdf file written without data (master writes empty
    ``*momentum_tensor.vtkhdf``) is returned, and the reader's FormatError comes when it is
    opened.
    """
    fmt = normalize_format(format)
    stem = quantity_stem(key)
    found = scan(directory)
    have = found.get(stem, {})
    if fmt is not None:
        if fmt in have:
            return have[fmt]
        other = "" if not have else f" (there is {stem}{EXTENSIONS[next(iter(have))]})"
        raise FileNotFoundError(
            f"{directory}: no {stem}{EXTENSIONS[fmt]} for quantity {key!r}{other}; the "
            f"directory has: {describe(found)}"
        )
    if len(have) > 1:
        raise ValueError(
            f"{directory}: both {stem}.h5 and {stem}.vtkhdf exist for quantity {key!r}; say "
            "which one to read with format='h5' or format='vtkhdf' (e.g. Run(path, "
            "format='vtkhdf'))"
        )
    if not have:
        raise FileNotFoundError(
            f"{directory}: no file for quantity {key!r} ({stem}.h5 or {stem}.vtkhdf); the "
            f"directory has: {describe(found)}"
        )
    return next(iter(have.values()))


__all__ = [
    "EXTENSIONS",
    "POP_KINDS",
    "SHORT_TO_STEM",
    "STEM_TO_SHORT",
    "describe",
    "find_quantity_file",
    "normalize_format",
    "pop_stem",
    "quantity_stem",
    "scan",
    "short_name",
]

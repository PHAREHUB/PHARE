"""Readers of PHARE output files.

A backend turns one file (one quantity) into small index tables and fills flat
buffers on request. Nothing here knows about fields, geometry or numpy
operations; h5py is imported only inside this sub-package.
"""

from __future__ import annotations

from .meta import Backend, FileMeta, LevelIndex

__all__ = ["Backend", "FileMeta", "LevelIndex", "open_backend"]


def open_backend(path, **kwargs):
    """Return the backend for `path`, chosen from the file extension."""
    path = str(path)
    if path.endswith(".vtkhdf"):
        from .vtkhdf import VtkhdfBackend

        return VtkhdfBackend(path, **kwargs)
    if path.endswith(".h5"):
        from .phareh5 import H5Backend

        return H5Backend(path, **kwargs)
    raise ValueError(f"{path}: not a .vtkhdf or .h5 file")

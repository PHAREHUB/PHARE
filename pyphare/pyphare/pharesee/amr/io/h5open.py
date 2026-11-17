"""Open an HDF5 file read-only, the same way for the vtkhdf and phareh5 readers.

Files are opened without HDF5 file locking: a running simulation may hold the file
open to write the next dump, and a lock taken by a reader must not get in its way.

HDF5 refuses to open a file a second time in one Python process with other lock
settings than a handle that is still open. That happens when another reader of the
same session (Ulysse's ``vtkamr``, master's pharesee, plain ``h5py.File(path)``)
keeps the file open with h5py's default settings. The open is then retried with
those default settings, which HDF5 accepts.
"""

from __future__ import annotations

import os

from ..errors import PhareseeIOError

#: text of the HDF5 error raised when an open's lock settings differ from those of a
#: handle of this process that keeps the file open
_LOCK_CLASH = "file locking"


def open_h5(path):
    """``h5py.File(path, "r")`` without file locking; with h5py's default locking when
    another handle of this process already holds the file that way.

    Raises PhareseeIOError naming the file and the reason (missing, empty, truncated
    or not HDF5, or held elsewhere in this session with other lock settings).
    """
    import h5py

    try:
        return h5py.File(path, "r", locking=False)
    except TypeError:  # older h5py without the locking argument
        try:
            return h5py.File(path, "r")
        except OSError as e:
            raise _open_error(path, e) from e
    except OSError as e:
        if _LOCK_CLASH not in str(e):
            raise _open_error(path, e) from e
        first = e
    try:
        return h5py.File(path, "r")
    except OSError as e:
        raise PhareseeIOError(
            f"{path}: cannot open as HDF5: the file is already open elsewhere in this "
            f"Python session with other HDF5 lock settings ({first}), and opening it with "
            f"h5py's default settings failed too ({e}). Close the other handle and retry"
        ) from e


def _open_error(path, e):
    try:
        size = os.path.getsize(path)
    except OSError:
        size = None
    what = ("empty file (0 bytes)" if size == 0 else
            "file missing" if size is None else "truncated or not an HDF5 file")
    return PhareseeIOError(
        f"{path}: cannot open as HDF5: {what} ({e}). If a simulation is still writing "
        "it, wait for the dump to finish; otherwise the file is damaged"
    )

"""What a backend returns and promises: file metadata, per-level index tables, and
the Backend protocol (kept apart from ``io/__init__`` so that the readers can import
it without importing the package that loads them)."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Protocol

import numpy as np


@dataclass(frozen=True)
class FileMeta:
    """What a file says about itself, read once when it is opened.

    path             absolute path of the file
    format           "vtkhdf" or "phareh5"
    quantity         name of the quantity as in the file name (e.g. "EM_B")
    dim              1, 2 or 3
    dl0              (dim,) float64 mesh size of level 0, rounded to the
                     shortest decimal that round-trips the file's dtype
    domain_cells0    (dim,) int64 number of level-0 cells
    periodic         (dim,) bool
    interp_order     int, or None when the file does not have it (MHD)
    model            "hybrid" or "mhd"
    names            component names in file order, e.g. ("x", "y", "z")
    dtype            numpy dtype of the field values in the file
    colocated        True when all components of a patch are stored together
    writer_signature short text (chunk shape, dtypes, read strategy) for repr
    attrs            other file attributes (pop_mass, pop_charge, ...)
    """

    path: str
    format: str
    quantity: str
    dim: int
    dl0: np.ndarray
    domain_cells0: np.ndarray
    periodic: np.ndarray
    interp_order: int | None
    model: str
    names: tuple
    dtype: np.dtype
    colocated: bool
    writer_signature: str = ""
    attrs: dict = field(default_factory=dict)


@dataclass(frozen=True)
class LevelIndex:
    """Patches of one level at one time, in FILE order.

    lower, upper : (npatch, dim) int64 inclusive cell boxes
    payload      : backend-private tables (row starts, byte offsets, ...)
    """

    level: int
    lower: np.ndarray
    upper: np.ndarray
    payload: Any = None


class Backend(Protocol):
    """What every file reader provides (see DESIGN section 5)."""

    def meta(self) -> FileMeta:
        """What the file says about itself (FileMeta), read once at opening."""

    def times(self) -> np.ndarray:
        """Sorted, unique dump times as float64."""

    def index(self, time: float) -> tuple:
        """Tuple of LevelIndex for the levels that have patches at `time`."""

    def read_into(self, time, level, file_ids, out, dst_start) -> None:
        """Fill `out` with the patches `file_ids` of `level` at `time`.

        file_ids  sorted int64 file-order patch indices
        out       {component name: flat buffer} for a non-empty subset of
                  meta().names: only those components are read (all of them in
                  one pass when `out` holds every name)
        dst_start int64 (len(file_ids),): node offset of each patch in `out`
        Each patch is written x fastest (Fortran order), nodes lower..upper+1,
        in the dtype of `out`. Values do not depend on which other components
        are read with them. Raises ValueError/KeyError for an empty `out` or an
        unknown name, PhareseeIOError on read failure; on failure the caller
        must consider nothing loaded.
        """

    def nbytes(self, time, level, file_ids, names=None) -> int:
        """Bytes read_into will bring into memory for components `names`
        (None: all of meta().names). No I/O."""

    def spec(self) -> tuple:
        """Picklable description from which open_backend can reopen the file."""

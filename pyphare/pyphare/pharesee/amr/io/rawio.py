"""Raw, threaded reads of byte ranges, and file identity checks.

Format-agnostic: used by the vtkhdf reader (chunk spans found by walking the
HDF5 chunk B-tree) and by the phareh5 reader (contiguous datasets at known
byte offsets). Nothing here knows about HDF5.

Why: on network file systems (glusterfs, Lustre) one reader thread issuing
large reads gets ~200 MB/s cold; 8 threads issuing 512 KiB ``pread`` calls get
500-670 MB/s (measured, see the design notes). Larger blocks are slower.

Main entry points

``read_ranges(fd, offsets, lengths, out)``
    read many ``[offset, offset + length)`` byte ranges of an open file into
    memory, splitting large ranges into blocks spread over threads.
``file_identity(path)`` / ``check_identity(recorded)`` / ``modified_since``
    detect that a file was replaced or truncated since it was indexed (a
    file that only grew, e.g. a running simulation appending, is accepted),
    or modified in place (the backends then re-validate what they indexed).
``stat_signature(fd)``
    (size, mtime) of an open file, to detect a file modified during a read.
``env_int`` / ``env_float`` / ``env_flag``
    strict parsing of the ``PHARESEE_*`` environment variables (a bad value
    raises ValueError naming the variable and the accepted values).
``time_tolerance`` / ``nearest_time``
    the time-matching rule shared by both backends.

No thread pool and no file handle exist at import time: a pool is created
per call (or passed in by the caller), so the module is fork-safe.
"""

from __future__ import annotations

import os
from concurrent.futures import ThreadPoolExecutor
from typing import NamedTuple

import numpy as np

from ..errors import PhareseeFileChanged, PhareseeIOError

#: default size of one read call; 128-512 KiB is the measured optimum
BLOCK_BYTES = 512 << 10
#: default number of reader threads, overridden by PHARESEE_IO_THREADS
DEFAULT_THREADS = 8
#: smallest time-matching tolerance (the phareh5 time keys have 10 decimals)
MIN_TIME_TOLERANCE = 5e-11

_TRUE = ("1", "true", "yes", "on")
_FALSE = ("0", "false", "no", "off")


# ---------------------------------------------------------------- environment


def env_int(name, default, minimum=None):
    """Integer value of environment variable `name` (`default` when unset or empty).

    Values below `minimum` are raised to it. Raises ValueError naming the
    variable when the value is not a whole number.
    """
    value = os.environ.get(name, "").strip()
    if not value:
        return default
    try:
        n = int(value)
    except ValueError:
        raise ValueError(
            f"{name}={value!r} is not a whole number; unset it or give e.g. {name}={default}"
        ) from None
    return n if minimum is None else max(minimum, n)


def env_float(name, default):
    """Float value of environment variable `name` (`default` when unset or empty).

    Raises ValueError naming the variable when the value is not a number.
    """
    value = os.environ.get(name, "").strip()
    if not value:
        return default
    try:
        return float(value)
    except ValueError:
        raise ValueError(
            f"{name}={value!r} is not a number; unset it or give e.g. {name}={default}"
        ) from None


def env_flag(name, default):
    """Boolean value of environment variable `name` (`default` when unset or empty).

    Accepted (any case): 1, true, yes, on / 0, false, no, off. Raises
    ValueError naming the variable and the accepted values otherwise.
    """
    value = os.environ.get(name, "").strip().lower()
    if not value:
        return bool(default)
    if value in _TRUE:
        return True
    if value in _FALSE:
        return False
    raise ValueError(
        f"{name}={os.environ.get(name)!r} is not understood; use one of "
        f"{', '.join(_TRUE)} (on) or {', '.join(_FALSE)} (off), or unset it"
    )


def io_threads():
    """Number of reader threads: ``PHARESEE_IO_THREADS`` (default 8, at least 1).

    Raises ValueError when the variable is set to something else than a whole number.
    """
    return env_int("PHARESEE_IO_THREADS", DEFAULT_THREADS, minimum=1)


def raw_io_enabled():
    """False when ``PHARESEE_RAW_IO`` is off (0/false/no/off); ValueError on other values."""
    return env_flag("PHARESEE_RAW_IO", True)


# ---------------------------------------------------------------- times


def time_tolerance(t, file_dtype=np.float64, given=None):
    """Matching tolerance between a requested time and a stored time near `t`.

    max(2 ulp of `file_dtype` at t, 5e-11), and also at least 2 ulp at t of the
    dtype of `given` when `given` (the value the caller passed) is a numpy
    floating scalar or 0-d array narrower than float64: ``np.float32(0.2)`` is
    0.20000000298..., which must still find a dump stored as 0.2 in float64.
    """
    t = abs(float(t))
    tol = max(2.0 * float(np.spacing(np.asarray(t, dtype=np.dtype(file_dtype)))),
              MIN_TIME_TOLERANCE)
    dtype = getattr(given, "dtype", None)
    if dtype is not None and dtype.kind == "f" and dtype.itemsize < 8:
        tol = max(tol, 2.0 * float(np.spacing(np.asarray(t, dtype=dtype))))
    return tol


def nearest_time(times, time, file_dtype=np.float64):
    """Index of the entry of the sorted float64 array `times` matching `time`, or None.

    Matching rule: ``time_tolerance`` (the dtype of `time` counts when it is a
    numpy float narrower than float64). Raises TypeError if `time` is not a number.
    """
    try:
        t = float(time)
    except (TypeError, ValueError):
        raise TypeError(f"a time must be a number, got {time!r}") from None
    if len(times) == 0 or not np.isfinite(t):
        return None
    i = int(np.searchsorted(times, t))
    cand = [j for j in (i - 1, i) if 0 <= j < len(times)]
    best = min(cand, key=lambda j: abs(times[j] - t))
    if abs(times[best] - t) <= time_tolerance(t, file_dtype, time):
        return best
    return None


# ---------------------------------------------------------------- identity


class FileIdentity(NamedTuple):
    """What identifies a file on disk at the time it was indexed.

    realpath     absolute path with symlinks resolved
    st_ino       inode number (a replaced file gets a new inode)
    st_size      size in bytes (a running simulation makes it grow)
    st_mtime_ns  modification time in ns
    """

    realpath: str
    st_ino: int
    st_size: int
    st_mtime_ns: int


def file_identity(path):
    """Return the FileIdentity of `path` (one ``os.stat``).

    Raises PhareseeIOError if the file cannot be stat'ed.
    """
    real = os.path.realpath(path)
    try:
        st = os.stat(real)
    except OSError as e:
        raise PhareseeIOError(f"{path}: cannot access the file ({e})") from e
    return FileIdentity(real, int(st.st_ino), int(st.st_size), int(st.st_mtime_ns))


def check_identity(recorded, path=None):
    """Check that the file recorded in `recorded` is still the same file.

    Accepted: same inode and size >= recorded size (the file grew: a running
    simulation appends). Whether what was indexed is still valid is the
    caller's check (``modified_since``, then a re-validation of its indexes).
    Raises PhareseeFileChanged when the file is gone, has a different inode
    (replaced) or is smaller (truncated). Returns the current FileIdentity.
    `path` defaults to ``recorded.realpath``.
    """
    path = recorded.realpath if path is None else path
    try:
        now = file_identity(path)
    except PhareseeIOError as e:
        raise PhareseeFileChanged(
            f"{path}: the file disappeared since it was opened; "
            "create a new Run to read the new file"
        ) from e
    if now.st_ino != recorded.st_ino:
        raise PhareseeFileChanged(
            f"{path}: the file was replaced since it was opened (inode "
            f"{recorded.st_ino} -> {now.st_ino}); create a new Run to read the new file"
        )
    if now.st_size < recorded.st_size:
        raise PhareseeFileChanged(
            f"{path}: the file shrank since it was opened ({recorded.st_size} -> "
            f"{now.st_size} bytes): it was truncated or rewritten; create a new Run to read "
            "the new file"
        )
    return now


def modified_since(recorded, now):
    """True when size or modification time differ between two FileIdentity of one inode.

    A modified file is either a running simulation that appended a dump (what
    was indexed is unchanged) or a file rewritten in place (e.g. a simulation
    re-run into the same folder): the backends then re-validate every cached
    index before using it (see ``rewritten``).
    """
    return now.st_size != recorded.st_size or now.st_mtime_ns != recorded.st_mtime_ns


def rewritten(path, what):
    """PhareseeFileChanged for a file rewritten in place; `what` says what differs."""
    return PhareseeFileChanged(
        f"{path}: the file was rewritten since it was opened ({what}); the data indexed "
        "earlier is gone. Create a new Run (or open the file again) to read the new content"
    )


def stat_signature(fd):
    """(size, mtime_ns) of an open file descriptor; compare before/after a read."""
    st = os.fstat(fd)
    return int(st.st_size), int(st.st_mtime_ns)


def open_readonly(path):
    """Open `path` read-only (O_CLOEXEC) and return the file descriptor.

    The caller closes it with ``os.close``. Raises PhareseeIOError.
    """
    flags = os.O_RDONLY | getattr(os, "O_CLOEXEC", 0)
    try:
        return os.open(path, flags)
    except OSError as e:
        raise PhareseeIOError(f"{path}: cannot open for reading ({e})") from e


# ---------------------------------------------------------------- reads


def _pread_into(fd, view, offset):
    """Fill the writable memoryview `view` with bytes from `offset`.

    Loops on short reads; raises PhareseeIOError at end of file.
    """
    want = view.nbytes
    done = 0
    while done < want:
        part = view[done:]
        if _HAVE_PREADV:
            n = os.preadv(fd, [part], offset + done)
        else:  # pragma: no cover - non-Linux fallback
            data = os.pread(fd, part.nbytes, offset + done)
            n = len(data)
            part[:n] = data
        if n <= 0:
            raise PhareseeIOError(
                f"short read: {done} of {want} bytes at offset {offset} (end of file reached)"
            )
        done += n
    return want


_HAVE_PREADV = hasattr(os, "preadv")


def _as_byte_view(buf):
    """Writable 1-D uint8 memoryview over `buf` (bytearray, ndarray, memoryview)."""
    if isinstance(buf, np.ndarray):
        if not buf.flags.c_contiguous:
            raise ValueError("read_ranges: output arrays must be C-contiguous")
        if not buf.flags.writeable:
            raise ValueError("read_ranges: output array is read-only")
        return memoryview(buf.reshape(-1).view(np.uint8))
    view = memoryview(buf)
    if view.readonly:
        raise ValueError("read_ranges: output buffer is read-only")
    return view.cast("B")


def _make_jobs(offsets, lengths, views, block):
    """Split each range into (file offset, destination view) jobs of <= block bytes."""
    jobs = []
    for off, ln, view in zip(offsets.tolist(), lengths.tolist(), views):
        for p in range(0, ln, block):
            k = min(block, ln - p)
            jobs.append((off + p, view[p : p + k]))
    return jobs


def _make_batches(jobs, block):
    """Group consecutive jobs into batches of about `block` bytes (at least one job each).

    Batches are handed to the threads in file order, so at any time the threads read
    neighbouring blocks: one moving read front. Measured on a network file system
    (gluster, cold): 2.3x faster than giving each thread its own distant part of the
    range. Grouping only matters for many small ranges (it bounds the per-call overhead).
    """
    batches, current, size = [], [], 0
    for job in jobs:
        current.append(job)
        size += job[1].nbytes
        if size >= block:
            batches.append(current)
            current, size = [], 0
    if current:
        batches.append(current)
    return batches


def _run_batch(fd, batch):
    total = 0
    for off, view in batch:
        total += _pread_into(fd, view, off)
    return total


def read_ranges(
    fd, offsets, lengths, out=None, *, nthreads=None, block=BLOCK_BYTES, executor=None
):
    """Read byte ranges ``[offsets[i], offsets[i] + lengths[i])`` of an open file.

    fd        file descriptor opened for reading (see ``open_readonly``)
    offsets   int64 array-like, byte offsets in the file
    lengths   int64 array-like, byte counts (same length as offsets)
    out       None: a new uint8 ndarray receives the ranges one after the other;
              one writable buffer (bytearray, C-contiguous ndarray, memoryview)
              of at least ``sum(lengths)`` bytes: same, written in place;
              a list of writable buffers: range i goes to ``out[i]`` (each
              at least ``lengths[i]`` bytes).
    nthreads  reader threads (default ``io_threads()``); 1 reads inline
    block     maximum bytes per system call (default 512 KiB)
    executor  optional ``concurrent.futures`` executor to run the reads on
              (lets a caller reuse one pool for several calls); when None a
              pool is created for this call and shut down at the end

    Returns `out` (the new array when `out` is None). Ranges may be given in
    any order; reads are issued in that order. Raises PhareseeIOError on a
    short read (range past the end of file) or an OS error.
    """
    offsets = np.asarray(offsets, dtype=np.int64).reshape(-1)
    lengths = np.asarray(lengths, dtype=np.int64).reshape(-1)
    if offsets.shape != lengths.shape:
        raise ValueError("read_ranges: offsets and lengths differ in length")
    if (lengths < 0).any() or (offsets < 0).any():
        raise ValueError("read_ranges: negative offset or length")
    block = max(1, int(block))
    total = int(lengths.sum())

    if out is None:
        out = np.empty(total, dtype=np.uint8)
    if isinstance(out, (list, tuple)):
        if len(out) != len(offsets):
            raise ValueError("read_ranges: one output buffer per range is required")
        views = [_as_byte_view(b) for b in out]
        for v, ln in zip(views, lengths.tolist()):
            if v.nbytes < ln:
                raise ValueError("read_ranges: an output buffer is smaller than its range")
    else:
        whole = _as_byte_view(out)
        if whole.nbytes < total:
            raise ValueError(
                f"read_ranges: output buffer has {whole.nbytes} bytes, {total} needed"
            )
        starts = np.concatenate([[0], np.cumsum(lengths)[:-1]]).tolist() if len(lengths) else []
        views = [whole[s : s + ln] for s, ln in zip(starts, lengths.tolist())]

    jobs = _make_jobs(offsets, lengths, views, block)
    if not jobs:
        return out
    nthreads = io_threads() if nthreads is None else max(1, int(nthreads))
    try:
        if (nthreads == 1 and executor is None) or len(jobs) == 1:
            _run_batch(fd, jobs)
            return out
        batches = _make_batches(jobs, block)
        if executor is not None:
            for fut in [executor.submit(_run_batch, fd, b) for b in batches]:
                fut.result()
        else:
            with ThreadPoolExecutor(max_workers=min(nthreads, len(batches))) as ex:
                for fut in [ex.submit(_run_batch, fd, b) for b in batches]:
                    fut.result()
    except PhareseeIOError:
        raise
    except OSError as e:
        raise PhareseeIOError(f"read failed: {e}") from e
    return out

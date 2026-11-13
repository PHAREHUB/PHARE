"""On-disk cache of the phareh5 index of one time of one file.

Walking the patch groups of one time of a production phareh5 file costs
0.4 ms per patch even when the file is in the page cache (seconds per time on
10^4 patches, much more cold). The tables that walk produces are small, so they
are kept on disk, one ``.npz`` file per (file, time), and reused by any later
process.

Location: ``$PHARESEE_CACHE_DIR`` if set, else ``$XDG_CACHE_HOME/pharesee``,
else ``~/.cache/pharesee``. Nothing is ever written next to the data.

Key: (real path, inode, time key, address of the ``t/<time>`` object header in
the file, format version). The entry is valid only if a hash of the first 4 KiB
of the file at that address is the one stored with it: a restarted run that
rewrites the same file (same inode) gets new headers, hence a miss.

Any problem while reading an entry (missing, truncated, corrupted, written by
another version, other key) is a miss: the caller walks the file again. Entries
are written atomically (temporary file + rename) so that concurrent processes
never see half an entry. The total size is capped (LRU on modification time,
refreshed at each hit) at ``PHARESEE_CACHE_MB`` MiB (default 1024).
``PHARESEE_CACHE=0`` disables reading and writing.
"""

from __future__ import annotations

import contextlib
import hashlib
import os
import tempfile

import numpy as np

from .rawio import env_flag, env_float

FORMAT_VERSION = "phareh5-index-v1"
HEADER_BYTES = 4096
_KEY_FIELD = "__cache_key__"


def enabled() -> bool:
    """False when the environment variable ``PHARESEE_CACHE`` is off (0, false, no, off).

    Raises ValueError for a value that is neither on nor off.
    """
    return env_flag("PHARESEE_CACHE", True)


def cache_dir() -> str:
    """Directory of the cache entries (not created here)."""
    explicit = os.environ.get("PHARESEE_CACHE_DIR")
    if explicit:
        return explicit
    base = os.environ.get("XDG_CACHE_HOME") or os.path.join(os.path.expanduser("~"), ".cache")
    return os.path.join(base, "pharesee")


def cap_bytes() -> int:
    """Size cap of the cache directory in bytes, from ``PHARESEE_CACHE_MB`` (default 1024).

    Raises ValueError when the variable is not a number.
    """
    mb = env_float("PHARESEE_CACHE_MB", 1024.0)
    return max(0, int(mb * (1 << 20)))


def header_digest(path: str, address: int, nbytes: int = HEADER_BYTES) -> str:
    """sha256 hex digest of the `nbytes` bytes of `path` starting at byte `address`.

    Reads at most `nbytes` bytes (fewer at the end of the file). Raises OSError.
    """
    fd = os.open(path, os.O_RDONLY)
    try:
        data = os.pread(fd, nbytes, int(address))
    finally:
        os.close(fd)
    return hashlib.sha256(data).hexdigest()


def _key_text(realpath, inode, time_key, address) -> str:
    return f"{FORMAT_VERSION}\0{realpath}\0{int(inode)}\0{time_key}\0{int(address)}"


def entry_path(realpath, inode, time_key, address) -> str:
    """Path of the cache file for this key (whether it exists or not)."""
    name = hashlib.sha256(_key_text(realpath, inode, time_key, address).encode()).hexdigest()
    return os.path.join(cache_dir(), name[:40] + ".npz")


def load(realpath, inode, time_key, address, digest) -> dict | None:
    """Arrays stored for this key, or None (disabled, missing, invalid, any error).

    `digest` is `header_digest(realpath, address)` computed by the caller now;
    the entry is used only if it was stored with the same digest.
    """
    if not enabled():
        return None
    path = entry_path(realpath, inode, time_key, address)
    try:
        with open(path, "rb") as fh, np.load(fh, allow_pickle=False) as npz:
            stored = npz[_KEY_FIELD]
            expected = np.array([_key_text(realpath, inode, time_key, address), digest])
            if stored.shape != expected.shape or not np.array_equal(stored, expected):
                return None
            arrays = {name: npz[name] for name in npz.files if name != _KEY_FIELD}
    except Exception:  # noqa: BLE001 - any unreadable entry is a miss by design
        return None
    with contextlib.suppress(OSError):
        os.utime(path)  # LRU marker
    return arrays


def store(realpath, inode, time_key, address, digest, arrays: dict) -> bool:
    """Write the entry for this key atomically, then trim the cache to its cap.

    `arrays` maps names (not starting with "__") to numpy arrays (no object
    dtype). Returns True if the entry was written. Never raises for I/O
    problems: a cache that cannot be written is simply not used.
    """
    if not enabled():
        return False
    directory = cache_dir()
    final = entry_path(realpath, inode, time_key, address)
    payload = dict(arrays)
    payload[_KEY_FIELD] = np.array([_key_text(realpath, inode, time_key, address), digest])
    tmp = None
    try:
        os.makedirs(directory, mode=0o700, exist_ok=True)
        fd, tmp = tempfile.mkstemp(prefix=".tmp-", suffix=".npz", dir=directory)
        with os.fdopen(fd, "wb") as fh:
            np.savez(fh, **payload)
        os.replace(tmp, final)
        tmp = None
    except (OSError, ValueError, TypeError):
        return False
    finally:
        if tmp is not None:
            with contextlib.suppress(OSError):
                os.unlink(tmp)
    trim(keep=final)
    return True


def trim(keep: str | None = None) -> None:
    """Delete the least recently used entries until the cache fits its cap.

    `keep` (a path) is deleted last, only if it alone exceeds the cap.
    Errors (concurrent deletions, permissions) are ignored.
    """
    directory = cache_dir()
    try:
        entries = []
        with os.scandir(directory) as it:
            for e in it:
                if e.name.endswith(".npz") and not e.name.startswith(".tmp-"):
                    st = e.stat()
                    entries.append((st.st_mtime_ns, st.st_size, e.path))
    except OSError:
        return
    total = sum(size for _, size, _ in entries)
    cap = cap_bytes()
    if total <= cap:
        return
    entries.sort(key=lambda e: (e[2] == keep, e[0]))
    for _, size, path in entries:
        if total <= cap:
            break
        try:
            os.unlink(path)
            total -= size
        except OSError:
            # already deleted by a concurrent prune, or not deletable: keep pruning the rest
            pass

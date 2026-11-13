"""Memory checks before reads and allocations.

``MemoryBudget`` finds how many bytes this process may still use: the minimum over
``/proc/meminfo`` MemAvailable, the cgroup v2 limits of every ancestor cgroup, the cgroup v1
limit, and ``RLIMIT_AS`` minus the current virtual size (lazy buffers count there). The
environment variable ``PHARESEE_MEMORY_LIMIT`` (bytes, or "8G", "500M", ...) replaces all of
them. ``check`` turns a size into a verdict: below 80 % of what is available silently
passes, 80-100 % emits a ``ResourceWarning``, above raises ``PhareseeMemoryError`` unless
``force=True``.

Cost: reading the system takes ~7 file opens (tens of microseconds), too slow to repeat on
every small arithmetic operation. A reading is therefore reused for ``ttl`` seconds (0.5 s),
and ``check`` answers "ok" without reading anything for requests below 1 % of the last
reading and below 64 MiB. Requests of 64 MiB or more, and any request whose verdict would
not be "ok", are always decided on a fresh reading.

Process pools: every worker sees the whole machine and would each think it may use all of
it. Set ``PHARESEE_MEMORY_LIMIT`` in each worker (e.g. total / number of workers).
"""

from __future__ import annotations

import os
import re
import sys
import time
import warnings

import numpy as np

from .errors import PhareseeMemoryError

ENV_LIMIT = "PHARESEE_MEMORY_LIMIT"
WARN_FRACTION = 0.8
SMALL_BYTES = 64 * 2**20  # below this and below SMALL_FRACTION of the last reading:
SMALL_FRACTION = 0.01  # "ok" without reading the system again
DEFAULT_TTL = 0.5  # seconds a system reading is reused
_UNITS = {"": 1, "K": 2**10, "M": 2**20, "G": 2**30, "T": 2**40, "P": 2**50}


def parse_size(text):
    """Bytes (int) from ``"123456"``, ``"500M"``, ``"8G"``, ``"1.5GiB"``, ``"2 TB"``.

    Units are powers of 1024 (K, M, G, T, P; an optional "B" or "iB" suffix is ignored).
    Raises ValueError with the accepted forms otherwise.
    """
    m = re.fullmatch(r"\s*([0-9]*\.?[0-9]+)\s*([KMGTP]?)(I?B)?\s*", str(text).upper())
    if m is None:
        raise ValueError(
            f"{text!r} is not a memory size: give bytes (e.g. 8000000000) or a number with "
            "K, M, G or T (e.g. 8G, 500M)"
        )
    return int(float(m.group(1)) * _UNITS[m.group(2)])


def format_bytes(n):
    """Short human-readable size, e.g. ``'1.50 GiB'``."""
    n = float(n)
    for unit in ("B", "KiB", "MiB", "GiB", "TiB"):
        if abs(n) < 1024 or unit == "TiB":
            return f"{n:.0f} {unit}" if unit == "B" else f"{n:.2f} {unit}"
        n /= 1024
    return f"{n:.2f} TiB"


def _read_text(path):
    try:
        with open(path) as f:
            return f.read()
    except OSError:
        return None


def _rlimit_as():
    try:
        import resource
    except ImportError:  # not on Linux/Unix
        return None
    soft, _ = resource.getrlimit(resource.RLIMIT_AS)
    return None if soft == resource.RLIM_INFINITY else int(soft)


def _kv(text, scale=1):
    """``{key: int}`` from 'key value [kB]' lines (meminfo, status, memory.stat)."""
    out = {}
    for line in (text or "").splitlines():
        parts = line.replace(":", " ").split()
        if len(parts) >= 2 and parts[1].lstrip("-").isdigit():
            out[parts[0]] = int(parts[1]) * scale
    return out


class MemoryBudget:
    """Bytes this process may still use, and the verdict for a planned allocation.

    Every system reader is injectable for tests:
      environ    mapping used for ``PHARESEE_MEMORY_LIMIT`` (default ``os.environ``)
      read_text  ``path -> str | None`` (None = missing file), default reads the file
      rlimit_as  ``() -> int | None`` soft RLIMIT_AS in bytes (None = unlimited)
      proc, cgroup_root  where /proc and the cgroup mount live
      ttl        seconds a reading of the system is reused by ``limit`` (default 0.5; 0 = never)
      clock      ``() -> float`` seconds, monotonic (default ``time.monotonic``)
    A change of ``PHARESEE_MEMORY_LIMIT`` invalidates the reused reading immediately.
    """

    def __init__(self, *, environ=None, read_text=None, rlimit_as=None, proc="/proc",
                 cgroup_root="/sys/fs/cgroup", ttl=DEFAULT_TTL, clock=None):
        self.environ = os.environ if environ is None else environ
        self.read_text = _read_text if read_text is None else read_text
        self.rlimit_as = _rlimit_as if rlimit_as is None else rlimit_as
        self.proc = proc
        self.cgroup_root = cgroup_root
        self.ttl = float(ttl)
        self.clock = time.monotonic if clock is None else clock
        self._reading = None  # (PHARESEE_MEMORY_LIMIT value, clock time, (bytes, source))

    # ---- individual sources: each returns (bytes, source text) or None
    def _status(self):
        return _kv(self.read_text(f"{self.proc}/self/status"), 1024)

    def _from_env(self):
        value = self.environ.get(ENV_LIMIT)
        if value is None or str(value).strip() == "":
            return None
        try:
            limit = parse_size(value)
        except ValueError as e:
            raise ValueError(f"{ENV_LIMIT}: {e}") from None
        rss = self._status().get("VmRSS", 0)
        # DESIGN-NOTE: the variable is read as a cap on this process's resident memory, so
        # what is still available is the cap minus the current resident size.
        return max(limit - rss, 0), f"{ENV_LIMIT}={value} minus resident size {format_bytes(rss)}"

    def _meminfo(self):
        info = _kv(self.read_text(f"{self.proc}/meminfo"), 1024)
        if "MemAvailable" not in info:
            return None
        return info["MemAvailable"], "MemAvailable in /proc/meminfo"

    def _cgroup_lines(self):
        lines = []
        for line in (self.read_text(f"{self.proc}/self/cgroup") or "").splitlines():
            parts = line.strip().split(":", 2)
            if len(parts) == 3:
                lines.append(parts)
        return lines

    def _cgroup_v2(self):
        rel = next((p[2] for p in self._cgroup_lines() if p[0] == "0" and p[1] == ""), None)
        if rel is None:
            return None
        best = None
        path = rel.rstrip("/")
        while True:  # the cgroup itself, then every ancestor up to the root
            d = f"{self.cgroup_root}{path}"
            limit = (self.read_text(f"{d}/memory.max") or "").strip()
            if limit and limit != "max":
                stat = _kv(self.read_text(f"{d}/memory.stat"))
                used = stat.get("anon", 0) + stat.get(
                    "kernel",
                    stat.get("kernel_stack", 0) + stat.get("slab", 0) + stat.get("pagetables", 0),
                )
                avail = max(int(limit) - used, 0)
                if best is None or avail < best[0]:
                    best = (avail, f"cgroup v2 memory.max of {path or '/'}")
            if path in ("", "/"):
                break
            path = path.rsplit("/", 1)[0]
        return best

    def _cgroup_v1(self):
        rel = next((p[2] for p in self._cgroup_lines() if "memory" in p[1].split(",")), None)
        if rel is None:
            return None
        d = f"{self.cgroup_root}/memory{rel.rstrip('/')}"
        text = (self.read_text(f"{d}/memory.limit_in_bytes") or "").strip()
        stat = _kv(self.read_text(f"{d}/memory.stat"))
        limits = [int(text)] if text.isdigit() else []
        if "hierarchical_memory_limit" in stat:
            limits.append(stat["hierarchical_memory_limit"])
        limits = [v for v in limits if v < 2**60]  # 2**63 - page = "no limit"
        if not limits:
            return None
        used = stat.get("total_rss", stat.get("rss", 0))
        return max(min(limits) - used, 0), f"cgroup v1 memory.limit_in_bytes of {rel}"

    def _rlimit(self):
        soft = self.rlimit_as()
        if soft is None:
            return None
        vm = self._status().get("VmSize", 0)
        return max(soft - vm, 0), "RLIMIT_AS (ulimit -v) minus the virtual size"

    # ---- public
    def sources(self):
        """All limits found, as a list of ``(bytes, source)``.

        With ``PHARESEE_MEMORY_LIMIT`` set, only that one.
        """
        env = self._from_env()
        if env is not None:
            return [env]
        found = [self._meminfo(), self._cgroup_v2(), self._cgroup_v1(), self._rlimit()]
        return [s for s in found if s is not None]

    def limit(self, max_age=None):
        """``(bytes, source)``: the smallest limit and where it comes from.

        max_age : seconds a previous reading may be reused (default ``ttl``; 0 = read now).
        """
        env = self.environ.get(ENV_LIMIT)
        max_age = self.ttl if max_age is None else float(max_age)
        now = self.clock()
        r = self._reading
        if r is not None and max_age > 0 and r[0] == env and now - r[1] < max_age:
            return r[2]
        found = self.sources()
        if not found:
            result = (sys.maxsize, "no memory limit found on this system")
        else:
            result = min(found, key=lambda s: s[0])
        self._reading = (env, now, result)
        return result

    def available(self):
        """Bytes this process may still allocate (int); a reading up to ``ttl`` old."""
        return self.limit()[0]

    def check(self, needed, what, force=False):
        """Verdict for allocating ``needed`` bytes for ``what`` (a short description).

        Returns "ok" (< 80 % of available), "warn" (80-100 %: ``ResourceWarning`` emitted),
        or "forced" (over, ``force=True``: ``ResourceWarning`` emitted). Over without force
        raises ``PhareseeMemoryError`` with the sizes, the source of the limit and how to
        override it.
        Readings: below ``SMALL_BYTES`` (64 MiB) and ``SMALL_FRACTION`` (1 %) of the last
        reading -> "ok" at once, nothing read (however old that reading is); from 64 MiB ->
        a fresh reading; in between -> a reading up to ``ttl`` old, read again before any
        verdict other than "ok".
        """
        needed = int(needed)
        r = self._reading
        if (r is not None and needed < SMALL_BYTES and needed < SMALL_FRACTION * r[2][0]
                and r[0] == self.environ.get(ENV_LIMIT)):
            return "ok"
        fresh = needed >= SMALL_BYTES
        avail, source = self.limit(max_age=0 if fresh else None)
        if needed < WARN_FRACTION * avail:
            return "ok"
        if not fresh:  # never warn or raise on a reused reading
            avail, source = self.limit(max_age=0)
            if needed < WARN_FRACTION * avail:
                return "ok"
        text = (
            f"{what} needs {format_bytes(needed)}; {format_bytes(avail)} available "
            f"({source})"
        )
        if needed <= avail:
            warnings.warn(f"{text}: close to the limit", ResourceWarning, stacklevel=3)
            return "warn"
        if force:
            warnings.warn(f"{text}: going ahead because force=True", ResourceWarning,
                          stacklevel=3)
            return "forced"
        raise PhareseeMemoryError(
            f"{text}. Select less data (box=..., one level, one component), free memory, or "
            f"override with force=True or {ENV_LIMIT}=<size, e.g. 16G>"
        )


_default_budget = MemoryBudget()


def default_budget():
    """The MemoryBudget used when none is given (readings reused for 0.5 s, see the module)."""
    return _default_budget


def check(needed, what, force=False):
    """``default_budget().check(needed, what, force)``."""
    return _default_budget.check(needed, what, force=force)


class ReadPlan:
    """Bytes a set of reads and results will bring into memory, before doing anything.

    Items are ``(what, bytes)``. Reads are sized with the backend (``Backend.nbytes``) for
    the patches not yet loaded; results with ``np.result_type`` of the operands' dtypes.
    No I/O, no allocation.
    """

    def __init__(self):
        self.items = []

    def add(self, nbytes, what):
        """Add a raw amount of bytes."""
        self.items.append((str(what), int(nbytes)))
        return self

    def add_read(self, level_data, slots=None, names=None):
        """Bytes to read the missing patches ``slots`` (store order; None = all) of a LevelData.

        names : None (every component) or the component name(s) that will be read.
        In-memory LevelData and already loaded (component, patch) pairs cost nothing
        (``LevelData.read_nbytes``).
        """
        ld = level_data
        if ld.backend is None:
            return self
        missing, _ = ld.missing(slots, names)
        if len(missing) == 0:
            return self
        n = ld.read_nbytes(slots, names)
        return self.add(n, f"read level {ld.geom.level} t={ld.time:g} ({len(missing)} patches)")

    def add_result(self, nnodes, ncomp, *dtypes, what="result"):
        """Bytes of a result of ``nnodes`` x ``ncomp`` values of ``np.result_type(*dtypes)``."""
        itemsize = np.result_type(*dtypes).itemsize
        return self.add(int(nnodes) * int(ncomp) * itemsize, what)

    @property
    def total(self):
        return sum(n for _, n in self.items)

    def check(self, budget=None, force=False, what="this operation"):
        """Run ``MemoryBudget.check`` on the total; returns the verdict."""
        budget = default_budget() if budget is None else budget
        return budget.check(self.total, what, force=force)

    def __repr__(self):
        lines = [f"ReadPlan: {format_bytes(self.total)}"]
        lines += [f"  {format_bytes(n):>12}  {w}" for w, n in self.items]
        return "\n".join(lines)

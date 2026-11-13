#!/usr/bin/env python3
"""Tests of pharesee.amr.plan: memory limits from injected system files, verdicts, ReadPlan."""

import unittest
import warnings

import numpy as np

from pyphare.pharesee.amr.errors import PhareseeMemoryError
from pyphare.pharesee.amr.geometry import LevelGeometry
from pyphare.pharesee.amr.plan import (
    SMALL_BYTES,
    MemoryBudget,
    ReadPlan,
    format_bytes,
    parse_size,
)
from pyphare.pharesee.amr.store import LevelData
from pyphare.pharesee.amr.testing import FakeBackend, fill, make_amr

GiB = 2**30


def fake_fs(files):
    return files.get


def budget(files, env=None, rlimit=None):
    return MemoryBudget(environ=env or {}, read_text=fake_fs(files), rlimit_as=lambda: rlimit)


MEMINFO = "MemTotal: 65000000 kB\nMemFree: 1000 kB\nMemAvailable: 8388608 kB\n"  # 8 GiB


class TestParse(unittest.TestCase):
    def test_sizes(self):
        self.assertEqual(parse_size("123456"), 123456)
        self.assertEqual(parse_size("8G"), 8 * GiB)
        self.assertEqual(parse_size("500M"), 500 * 2**20)
        self.assertEqual(parse_size("1.5GiB"), int(1.5 * GiB))
        self.assertEqual(parse_size(" 2 tb "), 2 * 2**40)
        self.assertEqual(parse_size("4k"), 4096)
        for bad in ("", "lots", "8 X", "-1G"):
            with self.assertRaises(ValueError):
                parse_size(bad)
        self.assertEqual(format_bytes(1536 * 2**20), "1.50 GiB")


class TestSources(unittest.TestCase):
    def test_meminfo_only(self):
        b = budget({"/proc/meminfo": MEMINFO})
        self.assertEqual(b.limit(), (8 * GiB, "MemAvailable in /proc/meminfo"))

    def test_nothing_found(self):
        avail, source = budget({}).limit()
        self.assertGreater(avail, 2**60)
        self.assertIn("no memory limit", source)

    def test_cgroup_v2_ancestors(self):
        files = {
            "/proc/meminfo": MEMINFO,
            "/proc/self/cgroup": "0::/job/step\n",
            "/sys/fs/cgroup/job/step/memory.max": "max\n",
            "/sys/fs/cgroup/job/memory.max": f"{2 * GiB}\n",
            "/sys/fs/cgroup/job/memory.stat": f"anon {GiB // 2}\nfile 999999999\nkernel {2**20}\n",
            "/sys/fs/cgroup/memory.max": f"{4 * GiB}\n",
            "/sys/fs/cgroup/memory.stat": "anon 0\n",
        }
        avail, source = budget(files).limit()
        self.assertEqual(avail, 2 * GiB - GiB // 2 - 2**20)  # page cache ("file") not counted
        self.assertIn("cgroup v2", source)
        self.assertIn("/job", source)
        # without "kernel", the kernel parts are summed
        files["/sys/fs/cgroup/job/memory.stat"] = (
            "anon 10\nslab 20\nkernel_stack 30\npagetables 40\n"
        )
        self.assertEqual(budget(files).available(), 2 * GiB - 100)

    def test_cgroup_v1(self):
        files = {
            "/proc/meminfo": MEMINFO,
            "/proc/self/cgroup": "5:cpu,cpuacct:/x\n4:memory:/slurm/job1\n",
            "/sys/fs/cgroup/memory/slurm/job1/memory.limit_in_bytes": f"{3 * GiB}\n",
            "/sys/fs/cgroup/memory/slurm/job1/memory.stat":
                f"total_rss {GiB}\nhierarchical_memory_limit 9223372036854771712\n",
        }
        avail, source = budget(files).limit()
        self.assertEqual(avail, 2 * GiB)
        self.assertIn("cgroup v1", source)
        files["/sys/fs/cgroup/memory/slurm/job1/memory.limit_in_bytes"] = "9223372036854771712"
        self.assertEqual(budget(files).limit()[1], "MemAvailable in /proc/meminfo")

    def test_rlimit_as(self):
        files = {"/proc/meminfo": MEMINFO, "/proc/self/status": "VmSize:\t1048576 kB\n"}
        avail, source = budget(files, rlimit=4 * GiB).limit()
        self.assertEqual(avail, 3 * GiB)
        self.assertIn("RLIMIT_AS", source)

    def test_env_overrides(self):
        files = {"/proc/meminfo": MEMINFO, "/proc/self/status": "VmRSS:\t1048576 kB\n"}
        b = budget(files, env={"PHARESEE_MEMORY_LIMIT": "100G"}, rlimit=GiB)
        avail, source = b.limit()
        self.assertEqual(avail, 99 * GiB)
        self.assertIn("PHARESEE_MEMORY_LIMIT", source)
        with self.assertRaises(ValueError):
            budget(files, env={"PHARESEE_MEMORY_LIMIT": "huge"}).limit()

    def test_real_system(self):
        avail = MemoryBudget().available()
        self.assertGreater(avail, 0)


class TestVerdicts(unittest.TestCase):
    def setUp(self):
        self.b = budget({"/proc/self/status": "VmRSS: 0 kB\n"},
                        env={"PHARESEE_MEMORY_LIMIT": "1000"})

    def test_levels(self):
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            self.assertEqual(self.b.check(799, "x"), "ok")
        with self.assertWarns(ResourceWarning):
            self.assertEqual(self.b.check(900, "x"), "warn")
        with self.assertRaises(PhareseeMemoryError) as cm:
            self.b.check(1001, "reading B")
        msg = str(cm.exception)
        for part in ("reading B", "1001 B", "1000 B", "PHARESEE_MEMORY_LIMIT", "force=True"):
            self.assertIn(part, msg)
        self.assertIsInstance(cm.exception, MemoryError)
        with self.assertWarns(ResourceWarning):
            self.assertEqual(self.b.check(5000, "x", force=True), "forced")


class TestReadPlan(unittest.TestCase):
    def test_sizes_missing_patches_only(self):
        hier = make_amr(2, 2, seed=1)
        src = fill(hier, {"x": lambda x, y: x, "y": lambda x, y: y}, dtype=np.float32)
        fb = FakeBackend(src)
        lds = [LevelData.from_backend(fb, g, 0.0) for g in hier.levels]
        plan = ReadPlan()
        for ld in lds:
            plan.add_read(ld)
        expected = sum(g.nnodes for g in hier.levels) * 4 * 2
        self.assertEqual(plan.total, expected)
        lds[0].ensure([0, 1])
        p2 = ReadPlan().add_read(lds[0])
        self.assertEqual(p2.total, (hier.levels[0].nnodes - hier.levels[0].size[:2].sum()) * 8)
        self.assertEqual(ReadPlan().add_read(src[0]).total, 0)  # in memory
        p3 = ReadPlan().add_result(1000, 3, np.float32, np.float64)
        self.assertEqual(p3.total, 24000)
        self.assertIn("ReadPlan", repr(p3))
        b = budget({"/proc/self/status": "VmRSS: 0 kB\n"}, env={"PHARESEE_MEMORY_LIMIT": "1K"})
        with self.assertRaises(PhareseeMemoryError):
            p3.check(b)

    def test_subset_level_uses_file_ids(self):
        # level restricted to some patches: sizes must be those of the FILE's patches
        hier = make_amr(2, 1, seed=4)
        g = hier.levels[0]
        src = fill(hier, {"x": lambda x, y: x}, dtype=np.float64)
        fb = FakeBackend(src)
        lo, up = g.file_order_boxes()
        size = np.prod(up - lo + 2, axis=1)  # nodes per patch, file order
        ids = np.argsort(size)[-3:][::-1].copy()  # the three largest patches
        sub = LevelGeometry.from_file_order(0, lo[ids], up[ids], g.dl, g.domain_cells,
                                            g.periodic)
        ld = LevelData.from_backend(fb, sub, 0.0, file_ids=ids)
        expected = int(size[ids].sum()) * 8
        wrong = int(size[:3].sum()) * 8  # what positions 0..2 of the file would give
        self.assertNotEqual(expected, wrong, "layout must have patches of different sizes")
        self.assertEqual(ReadPlan().add_read(ld).total, expected)
        ld.ensure([sub.store_index[0]])
        self.assertEqual(ReadPlan().add_read(ld).total, int(size[ids[1:]].sum()) * 8)
        self.assertEqual(ReadPlan().add_read(ld, slots=[sub.store_index[2]]).total,
                         int(size[ids[2]]) * 8)


class Clock:
    def __init__(self):
        self.t = 0.0

    def __call__(self):
        return self.t


class TestReuse(unittest.TestCase):
    """Readings reused for ttl seconds; small requests answered without reading."""

    def make(self, avail_kb=8388608, env=None):
        self.reads = []
        self.files = {"/proc/meminfo": f"MemAvailable: {avail_kb} kB\n",
                      "/proc/self/status": "VmRSS: 0 kB\n"}

        def read_text(path):
            self.reads.append(path)
            return self.files.get(path)

        self.clock = Clock()
        self.env = dict(env or {})
        return MemoryBudget(environ=self.env, read_text=read_text, rlimit_as=lambda: None,
                            clock=self.clock)

    def test_ttl(self):
        b = self.make()
        self.assertEqual(b.available(), 8 * GiB)
        n = len(self.reads)
        self.assertGreater(n, 0)
        self.files["/proc/meminfo"] = "MemAvailable: 1048576 kB\n"
        self.clock.t = 0.4
        self.assertEqual(b.available(), 8 * GiB)  # reused
        self.assertEqual(len(self.reads), n)
        self.clock.t = 0.6
        self.assertEqual(b.available(), GiB)  # older than ttl: read again
        self.assertEqual(len(self.reads), 2 * n)
        self.assertEqual(b.limit(max_age=0)[0], GiB)  # forced fresh reading
        self.assertEqual(len(self.reads), 3 * n)

    def test_env_change_invalidates(self):
        b = self.make()
        b.available()
        self.env["PHARESEE_MEMORY_LIMIT"] = "1000"
        self.assertEqual(b.limit(), (1000, "PHARESEE_MEMORY_LIMIT=1000 minus resident size 0 B"))
        # the early exit does not use a reading made under another limit
        b2 = self.make()
        b2.check(1000, "x")
        self.env["PHARESEE_MEMORY_LIMIT"] = "1K"
        with self.assertRaises(PhareseeMemoryError):
            b2.check(5000, "x")

    def test_small_requests_read_nothing(self):
        b = self.make()
        self.assertEqual(b.check(10, "x"), "ok")  # first call: must read
        n = len(self.reads)
        self.clock.t = 100.0  # however old the reading
        for size in (1, 2**20, SMALL_BYTES - 1):
            self.assertEqual(b.check(size, "x"), "ok")
        self.assertEqual(len(self.reads), n)

    def test_large_requests_read_fresh(self):
        b = self.make()
        b.check(10, "x")
        n = len(self.reads)
        self.assertEqual(b.check(SMALL_BYTES, "x"), "ok")
        self.assertEqual(len(self.reads), 2 * n)  # 64 MiB: fresh even within ttl
        # memory dropped since: the large request sees it at once
        self.files["/proc/meminfo"] = "MemAvailable: 65536 kB\n"  # 64 MiB
        with self.assertRaises(PhareseeMemoryError):
            b.check(SMALL_BYTES + 1, "x")

    def test_never_warns_on_a_reused_reading(self):
        b = self.make(avail_kb=10240)  # 10 MiB: 1 MiB is not "small" (> 1 %)
        b.check(1, "x")
        n = len(self.reads)
        self.assertEqual(b.check(2**20, "x"), "ok")  # reused reading, ok
        self.assertEqual(len(self.reads), n)
        self.files["/proc/meminfo"] = "MemAvailable: 100000 kB\n"
        self.assertEqual(b.check(9 * 2**20, "x"), "ok")  # stale says warn -> fresh says ok
        self.assertEqual(len(self.reads), 2 * n)
        self.files["/proc/meminfo"] = "MemAvailable: 1024 kB\n"
        self.clock.t = 1.0
        b.available()
        self.files["/proc/meminfo"] = "MemAvailable: 10240 kB\n"
        with self.assertWarns(ResourceWarning):  # decided on the fresh 10 MiB, not 1 MiB
            self.assertEqual(b.check(9 * 2**20, "x"), "warn")


if __name__ == "__main__":
    unittest.main()

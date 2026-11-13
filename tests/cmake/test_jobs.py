import os
import re
import subprocess
import tempfile
import unittest
from pathlib import Path

JOBS = Path(__file__).resolve().parents[2] / "res/cmake/jobs.cmake"
CMAKE = os.environ.get("CMAKE_COMMAND") or "cmake"
GIB = 1024**3


def run_cmake(body):
    with tempfile.TemporaryDirectory() as tmp:
        driver = Path(tmp) / "driver.cmake"
        driver.write_text(f'cmake_minimum_required(VERSION 3.20)\ninclude("{JOBS}")\n{body}\nmessage("RESULT=${{out}}")\n')
        proc = subprocess.run(
            [CMAKE, "-P", str(driver)], capture_output=True, text=True, check=False
        )
    return proc.returncode, proc.stdout + proc.stderr


def run_driver(body):
    returncode, output = run_cmake(body)
    match = re.search(r"RESULT=(.*)", output)
    if returncode != 0 or match is None:
        raise AssertionError(output)
    return match.group(1).strip()


def cgroup_limit(content):
    with tempfile.TemporaryDirectory() as tmp:
        path = Path(tmp) / "limit"
        if content is not None:
            path.write_text(content)
        return run_driver(f'phare_cgroup_limit_mb("{path}" out)')


def cgroup_memory(proc_cgroup, files):
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp) / "cgroup"
        root.mkdir()
        proc = Path(tmp) / "proc_self_cgroup"
        if proc_cgroup is not None:
            proc.write_text(proc_cgroup)
        for name, content in files.items():
            path = root / name
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(content)
        return run_driver(f'phare_cgroup_memory_mb("{proc}" "{root}" out)')


def compile_jobs(total, cgroup, per_module):
    return run_driver(f'phare_module_compile_jobs({total} "{cgroup}" {per_module} out)')


def effective_jobs(explicit, memory, per_module):
    return run_cmake(
        f'phare_effective_module_jobs("{explicit}" {memory} "{per_module}" out)'
    )


class CgroupLimitTest(unittest.TestCase):
    def test_max_is_unlimited(self):
        self.assertEqual(cgroup_limit("max\n"), "")

    def test_bytes_are_converted_to_mib(self):
        self.assertEqual(cgroup_limit("1073741824\n"), "1024")

    def test_missing_file_is_unlimited(self):
        self.assertEqual(cgroup_limit(None), "")

    def test_garbage_is_unlimited(self):
        self.assertEqual(cgroup_limit("abc"), "")


class CgroupMemoryTest(unittest.TestCase):
    def test_v2_limit_two_levels_up(self):
        files = {
            "a/b/c/memory.max": "max\n",
            "a/b/memory.max": "max\n",
            "a/memory.max": f"{8 * GIB}\n",
        }
        self.assertEqual(cgroup_memory("0::/a/b/c\n", files), "8192")

    def test_v2_smallest_level_wins(self):
        files = {
            "a/b/c/memory.max": f"{2 * GIB}\n",
            "a/b/memory.max": "max\n",
            "a/memory.max": f"{8 * GIB}\n",
        }
        self.assertEqual(cgroup_memory("0::/a/b/c\n", files), "2048")

    def test_v2_max_everywhere(self):
        files = {
            "a/b/c/memory.max": "max\n",
            "a/b/memory.max": "max\n",
            "a/memory.max": "max\n",
        }
        self.assertEqual(cgroup_memory("0::/a/b/c\n", files), "")

    def test_v1_memory_controller_path(self):
        files = {
            "memory/job/step/memory.limit_in_bytes": "9223372036854771712\n",
            "memory/job/memory.limit_in_bytes": f"{4 * GIB}\n",
            "memory/memory.limit_in_bytes": "9223372036854771712\n",
            "cpu/job/step/memory.limit_in_bytes": f"{1 * GIB}\n",
        }
        proc = "12:cpu,cpuacct:/job/step\n4:memory:/job/step\n1:name=systemd:/job\n"
        self.assertEqual(cgroup_memory(proc, files), "4096")

    def test_root_files_are_the_fallback(self):
        self.assertEqual(cgroup_memory(None, {"memory.max": f"{1 * GIB}\n"}), "1024")
        self.assertEqual(
            cgroup_memory(None, {"memory/memory.limit_in_bytes": f"{3 * GIB}\n"}),
            "3072",
        )

    def test_namespace_root(self):
        self.assertEqual(
            cgroup_memory("0::/\n", {"memory.max": f"{6 * GIB}\n"}), "6144"
        )

    def test_no_cgroup_is_unlimited(self):
        self.assertEqual(cgroup_memory(None, {}), "")


class ModuleCompileJobsTest(unittest.TestCase):
    def test_host_memory_only(self):
        self.assertEqual(compile_jobs(54000, "", 4096), "9")

    def test_smaller_cgroup_wins(self):
        self.assertEqual(compile_jobs(54000, 16384, 2048), "6")

    def test_larger_cgroup_ignored(self):
        self.assertEqual(compile_jobs(54000, 100000, 4096), "9")

    def test_at_least_one_job(self):
        self.assertEqual(compile_jobs(2000, "", 4096), "1")


class EffectiveModuleJobsTest(unittest.TestCase):
    def assertResult(self, explicit, memory, per_module, expected):
        returncode, output = effective_jobs(explicit, memory, per_module)
        self.assertEqual(returncode, 0, output)
        self.assertEqual(re.search(r"RESULT=(.*)", output).group(1).strip(), expected)

    def assertFails(self, explicit, memory, per_module, variable):
        returncode, output = effective_jobs(explicit, memory, per_module)
        self.assertNotEqual(returncode, 0, output)
        self.assertIn(variable, output)
        self.assertNotIn("RESULT=", output)

    def test_empty_is_from_memory(self):
        self.assertResult("", 16384, 2048, "6")

    def test_explicit_is_used_as_given(self):
        self.assertResult("40", 16384, 2048, "40")

    def test_zero_disables(self):
        self.assertResult("0", 16384, 2048, "0")

    def test_invalid_jobs_stop(self):
        for value in ["abc", "-1", "2.5", "4 jobs"]:
            with self.subTest(value=value):
                self.assertFails(value, 16384, 2048, "PHARE_MODULE_COMPILE_JOBS")

    def test_invalid_memory_stops(self):
        for value in ["0", "-512", "abc", ""]:
            with self.subTest(value=value):
                self.assertFails("", 16384, value, "PHARE_MODULE_COMPILE_MEMORY_MB")


if __name__ == "__main__":
    unittest.main()

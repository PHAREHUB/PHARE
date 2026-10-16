import re
import subprocess
import tempfile
import unittest
from pathlib import Path

JOBS = Path(__file__).resolve().parents[2] / "res/cmake/jobs.cmake"


def run_driver(body):
    with tempfile.TemporaryDirectory() as tmp:
        driver = Path(tmp) / "driver.cmake"
        driver.write_text(f'include("{JOBS}")\n{body}\nmessage("RESULT=${{out}}")\n')
        proc = subprocess.run(
            ["cmake", "-P", str(driver)], capture_output=True, text=True
        )
    output = proc.stdout + proc.stderr
    match = re.search(r"RESULT=(.*)", output)
    if proc.returncode != 0 or match is None:
        raise AssertionError(output)
    return match.group(1).strip()


def cgroup_limit(content):
    with tempfile.TemporaryDirectory() as tmp:
        path = Path(tmp) / "limit"
        if content is not None:
            path.write_text(content)
        return run_driver(f'phare_cgroup_limit_mb("{path}" out)')


def compile_jobs(total, cgroup, per_module):
    return run_driver(f'phare_module_compile_jobs({total} "{cgroup}" {per_module} out)')


class CgroupLimitTest(unittest.TestCase):
    def test_max_is_unlimited(self):
        self.assertEqual(cgroup_limit("max\n"), "")

    def test_bytes_are_converted_to_mib(self):
        self.assertEqual(cgroup_limit("1073741824\n"), "1024")

    def test_missing_file_is_unlimited(self):
        self.assertEqual(cgroup_limit(None), "")

    def test_garbage_is_unlimited(self):
        self.assertEqual(cgroup_limit("abc"), "")


class ModuleCompileJobsTest(unittest.TestCase):
    def test_host_memory_only(self):
        self.assertEqual(compile_jobs(54000, "", 4096), "13")

    def test_smaller_cgroup_wins(self):
        self.assertEqual(compile_jobs(54000, 8192, 4096), "2")

    def test_larger_cgroup_ignored(self):
        self.assertEqual(compile_jobs(54000, 100000, 4096), "13")

    def test_at_least_one_job(self):
        self.assertEqual(compile_jobs(2000, "", 4096), "1")


if __name__ == "__main__":
    unittest.main()

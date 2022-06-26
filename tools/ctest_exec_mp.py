

import os
import sys
from pathlib import Path

not_root_error = "script must be run from project root directory"
assert all([os.path.exists(d) for d in ["tests", "tools", "CMakeLists.txt"]]), not_root_error
root = Path(os.getcwd())
sys.path.insert(0, ".")

from tools.python3 import run, pushd
import tools.python3.cmake as cmake
import tools.python3.git as git


def run_tests_log_to_file(data_dir, n_cores, tests):
    """ ctest seems to merge stdout and stderr so we only need one output file """
    import concurrent.futures

    with concurrent.futures.ThreadPoolExecutor(max_workers=n_cores) as executor:
        jobs = [
          executor.submit(
            run, cmake.test_cmd(test, verbose=True),
            shell=True, capture_output=False, check=True,
            stdout=open(os.path.join(str(data_dir), test), "w"))
            for i, test in enumerate(tests)
        ]
        results = []
        for future in concurrent.futures.as_completed(jobs):
            try:
                results += [future.result()]
            except Exception as exc:
                print("run_mp generated an exception: %s" % exc)
        return results

def main():
    with pushd(os.path.join(str(root), "build")):
        current_git_hash = git.hashes(1)[0]
        data_dir = Path(os.path.join(str(root), "data_out", current_git_hash, "ctest_logs"))
        os.makedirs(str(data_dir), exist_ok=True)
        N_CORES= int(os.environ["N_CORES"]) if "N_CORES" in os.environ else 1
        print(f"Launching ctests with N_CORES {N_CORES}")
        run_tests_log_to_file(data_dir, N_CORES, cmake.list_tests())

if __name__ == "__main__":
    main()

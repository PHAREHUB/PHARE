
"""

To add another test case, copy/rename the file "test_cases/uniform.py", and configure it as you wish!
new test cases are identified by the file name, minus ".py" i.e. "uniform"

example execution might be
```
export PYTHONPATH=$PWD/build:$PWD/pyphare
python3 -u tools/bench/functional/run_gen_plots.py --samrai_dir=/mkn/r/llnl/samrai/master

```

"""

import os
import sys
import atexit
from distutils.util import strtobool
from pathlib import Path

import matplotlib as mpl
mpl.use('Agg') # without GUI

assert all([os.path.exists(d) for d in ["tests", "tools", "CMakeLists.txt"]])
root = Path(os.getcwd())  # expects project root!
sys.path.insert(0, ".")

# sys.path.insert(0, os.path.join(root, "pyphare"))
this_dir = Path(os.path.dirname(os.path.abspath(__file__)))
test_cases_dir = "test_cases"

from tools.python3 import binary_exists_on_path, scan_dir
import tools.python3.perf as perf
import tools.python3.git as git
import tools.python3.cmake as cmake

def find_test_case():
    test_cases = []
    for file_name in scan_dir(os.path.join(this_dir, test_cases_dir), files_only=True):
        file_name, file_ext = os.path.splitext(file_name)
        test_cases += [file_name]
    return test_cases

available_test_cases = find_test_case()

# defaults that can be modified, or overridden by cli
default_cli_args = {
    "build_test": True,
    "test_cases": available_test_cases,  # None = scan "generated" dir and run all files to run all tests for said file
    "samrai_dir": None,  # None = ${root}/subprojects/samrai
    "cxx_flags": "-O3 -g3 -march=native -mtune=native",
    "repeat_stat": 1,
    "build_top_N_commits": 1,
    "use_ninja": binary_exists_on_path("ninja"),
    "use_ccache": binary_exists_on_path("ccache"),
    "multithreaded": False,
    "use_found_build_for_top": False,
}


def parse_cli_args():
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--build_test",
        default=default_cli_args["build_test"],
        type=strtobool,
        help="if True, builds and tests $build_top_N_commits times, useful if you have already built and just want to mess with plotting for each set of perf results",
    )
    parser.add_argument("--test_cases", default=default_cli_args["test_cases"])
    parser.add_argument("--samrai_dir", default=default_cli_args["samrai_dir"])
    parser.add_argument("--cxx_flags", default=default_cli_args["cxx_flags"])
    parser.add_argument(
        "--repeat_stat",
        type=int,
        default=default_cli_args["repeat_stat"],
        help="run perf stat this many times on each test case",
    )
    parser.add_argument(
        "--build_top_N_commits",
        type=int,
        default=default_cli_args["build_top_N_commits"],
    )
    parser.add_argument(
        "--use_ninja", default=default_cli_args["use_ninja"], type=strtobool, help="default: used if found"
    )
    parser.add_argument(
        "--use_ccache", default=default_cli_args["use_ccache"], type=strtobool, help="default: used if found"
    )
    parser.add_argument(
        "--use_found_build_for_top", default=default_cli_args["use_found_build_for_top"], type=strtobool,
        help="If true, will not delete the build directory if exists for first set of tests"
    )
    parser.add_argument(
        "--multithreaded",
        type=strtobool,
        default=default_cli_args["multithreaded"],
        help="advanced option to run $repeat_stat threads at the same time for stat operation",
    )
    args = parser.parse_args()
    if isinstance(args.test_cases, str):
        args.test_cases = args.test_cases.split(",")
    return vars(args)


# can be modified
perf_events = [
    "duration_time",
    "cycles",
    "instructions",
    "cache-references",
    "cache-misses",
    "L1-dcache-loads",
    "L1-dcache-load-misses",
]
# "perf stat" can support more events than "perf record"
stat_events = perf_events + ["bus-cycles"]

# configurable but not so recommended.
phare_exe = "build/src/phare/phare-exe"
build_dir = Path(os.path.join(str(root), "build"))
data_dir = Path(os.path.join(str(root), "data_out"))
cmake_config_user = "-DdevMode=ON -DCMAKE_BUILD_TYPE=Release"

# not to be changed!
cmake_config_extra = cmake_config_user + " -Dbench=ON"
current_git_hash = git.hashes(1)[0]

def git_branch_reset_at_exit():
    current_git_branch = git.current_branch()

    @atexit.register
    def _reset_():
        git.checkout(current_git_branch)


def perf_bin(pyfile):
    return phare_exe + f" {pyfile}"


def py_file_to_module(bin_dir, py_file):
    return os.path.relpath(os.path.join(bin_dir, py_file)).replace(os.path.sep, ".")


def git_commit_data_dir(git_hash):
    return os.path.join(str(data_dir), git_hash)


def test_file_data_dir(test_case, test_file, git_hash):
    return os.path.join(git_commit_data_dir(git_hash), test_case, test_file)


def test_case_gen_dir(test_case):
    return os.path.join(this_dir, "generated", test_case)


def run_tests(test_cases, cli_args, git_hash):
    for test_case in test_cases:
        bin_dir = test_case_gen_dir(test_case)
        for file_name in scan_dir(bin_dir, files_only=True):
            file_name, file_ext = os.path.splitext(file_name)
            assert file_ext == ".py"
            bindata_dir = test_file_data_dir(test_case, file_name, git_hash)
            Path(bindata_dir).mkdir(parents=True, exist_ok=True)
            py_file_module = py_file_to_module(bin_dir, file_name)
            exe = perf_bin(py_file_module + file_ext)
            outs = [
                os.path.join(bindata_dir, f"{repeat}.data")
                for repeat in range(cli_args["repeat_stat"])
            ]
            if cli_args["multithreaded"]:
                perf.stat_mp(exe, stat_events, outs)
            else:
                for out in outs:
                    perf.stat(exe, stat_events, out)

            perf.record(exe, perf_events, os.path.join(bindata_dir, "record.data"))


def parse_test_results(test_case, cli_args, git_hash):
    import importlib
    results = {}
    bin_dir = test_case_gen_dir(test_case)
    for file_name in scan_dir(bin_dir, files_only=True):
        file_name, file_ext = os.path.splitext(file_name)
        bindata_dir = test_file_data_dir(test_case, file_name, git_hash)
        py_file_module = py_file_to_module(bin_dir, file_name)
        module = importlib.import_module(py_file_module)

        results[file_name] = {
          "stat_data": [
              perf.parse_stat_csv(os.path.join(bindata_dir, f"{repeat}.data"))
              for repeat in range(cli_args["repeat_stat"])
          ],
          "params" : module.params
        }

    return results


def retrieve_results(test_cases, cli_args, git_hash):
    """
    returns {"uniform" : {
      "uniform_1_1_256_11" : {"stat_data" : stat_data_files_list, "params" : params_for_file_dict}}}
    """

    return {
        test_case: parse_test_results(test_case, cli_args, git_hash)
        for test_case in test_cases
    }


def plot_scatter(test_case, cli_args, file_results, getter, type, git_hash):
    # https://stackoverflow.com/a/44985520
    import numpy as np
    import matplotlib.pyplot as plt

    nresults = len(file_results)
    labels = list(range(nresults))

    width = 0.4
    fig, ax = plt.subplots(figsize=(15, 15))
    for i, k in enumerate(file_results):
        v = file_results[k]["stat_data"]
        data = np.array(
            [getter(v[repeat]) for repeat in range(cli_args["repeat_stat"])]
        )

        x = np.ones(data.shape[0]) * i + (
            np.random.rand(data.shape[0]) * width - width / 2.0
        )
        ax.scatter(x, data, s=25)
        mean = data.mean()
        ax.plot([i - width / 2.0, i + width / 2.0], [mean, mean], color="k")

    ax.set_xticks(labels)
    ax.set_xticklabels(labels)
    lgd = ax.legend(
        [k for k, v in file_results.items()], loc="center left", bbox_to_anchor=(1, 0.5)
    )
    filename = f"{test_case}_{type}_scatter.png"
    out = f"{os.path.join(git_commit_data_dir(git_hash), filename)}"
    fig.savefig(out, bbox_extra_artists=(lgd,), bbox_inches="tight")


def plot_gen(test_cases, cli_args, git_hash):
    case_results = retrieve_results(test_cases, cli_args, git_hash)
    for test_case, file_results in case_results.items():
        plot_scatter(
            test_case,
            cli_args,
            file_results,
            lambda v0: float(v0["duration_time"]) * 1e-9,
            type="time",
            git_hash=git_hash,
        )
        plot_scatter(
            test_case,
            cli_args,
            file_results,
            lambda v0: float(v0["L1-dcache-load-misses"]),
            type="L1_miss",
            git_hash=git_hash,
        )


def cmake_clean_config_build(cli_args):
    import shutil

    os.chdir(str(root))
    if os.path.exists(str(build_dir)):
        shutil.rmtree(str(build_dir))
    build_dir.mkdir(exist_ok=False)
    os.chdir(str(build_dir))

    cmake.config(
        path="..",
        samrai_dir=cli_args["samrai_dir"],
        cxx_flags=cli_args["cxx_flags"],
        use_ninja=cli_args["use_ninja"],
        use_ccache=cli_args["use_ccache"],
        extra=cmake_config_extra,
    )
    cmake.build(use_ninja=cli_args["use_ninja"])
    os.chdir(str(root))


def build_test(test_cases, cli_args, git_hash):
    git.checkout(git_hash)
    should_build = not cli_args["use_found_build_for_top"] or (
      current_git_hash == git_hash and not os.path.exists(str(build_dir))
    )
    if should_build:
        cmake_clean_config_build(cli_args)
    run_tests(test_cases, cli_args, git_hash)


def verify_cli_args(cli_args):

    for key in list(default_cli_args.keys()):
        assert key in cli_args

    assert len(cli_args["test_cases"]) > 0
    assert cli_args["build_top_N_commits"] > 0
    return cli_args


def generate_test_cases(test_cases):
    import importlib
    def generate(test_case):
        module = importlib.import_module(f"tools.bench.functional.{test_cases_dir}.{test_case}")
        module.generate_all()

    for test_case in test_cases:
        if test_case not in available_test_cases:
            raise RuntimeError(f"test_case ({test_case}) does not exist in subdirectory {test_cases_dir}")
        generate(test_case)


def main():
    cli_args = verify_cli_args(parse_cli_args())
    perf.check()

    data_dir.mkdir(exist_ok=True)
    git_branch_reset_at_exit()

    test_cases = cli_args["test_cases"]
    generate_test_cases(test_cases)
    build_top_N_commits = cli_args["build_top_N_commits"]

    with open(os.path.join(str(data_dir), "glog"), "w") as out:
        out.write(git.log(build_top_N_commits + 10, use_short=True))

    hashes = git.hashes(build_top_N_commits)

    for git_hash in hashes:
        if cli_args["build_test"]:
            build_test(test_cases, cli_args, git_hash)
        plot_gen(test_cases, cli_args, git_hash)


if __name__ == "__main__":
    main()

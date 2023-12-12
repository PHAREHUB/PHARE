from tools.python3 import decode_bytes, run


def version():
    pass

# set to false if you do not want this
_USE_NINJA=True

# None == subprojects, needs building per commit if unset, so not great
_SAMRAI_DIR="/mkn/r/llnl/samrai/master"

def make_config_str(
    path, samrai_dir=None, cxx_flags=None, use_ccache=False, extra=""
):
    """
    FULL UNPORTABLE OPTIMIZATIONS = cxx_flags="-O3 -march=native -mtune=native"
    """
    samrai_dir = _SAMRAI_DIR if samrai_dir is None else samrai_dir
    samrai_dir = "" if samrai_dir is None else f"-DSAMRAI_ROOT={samrai_dir}"
    cxx_flags = "" if cxx_flags is None else f'-DCMAKE_CXX_FLAGS="{cxx_flags}"'
    ccache = "" if use_ccache is False else "-DCMAKE_CXX_COMPILER_LAUNCHER=ccache"
    ninja = "" if not _USE_NINJA else "-G Ninja"

    return f"cmake {path} {samrai_dir} {cxx_flags} {ninja} {ccache} {extra}"


def config(
    path, samrai_dir=None, cxx_flags=None, use_ccache=False, extra=""
):
    cmd = make_config_str(path, samrai_dir, cxx_flags, extra)
    run(cmd, capture_output=False)


def build(threads=1):
    run("ninja" if _USE_NINJA else f"make -j{threads}", capture_output=False)


def list_tests():
    proc = run("ctest -N", capture_output=True)
    out  = decode_bytes(proc.stdout)
    return [
        line.split(" ")[-1] for line in out.splitlines()[1:-2]
    ]


def test_cmd(test, verbose=False):
    cmd = f"ctest -R {test}"
    if verbose:
        cmd = f"{cmd} -V"
    return cmd


def run_test(test, verbose=False, capture_output=False):
    run(test_cmd(cmd, verbose=verbose), capture_output=capture_output)


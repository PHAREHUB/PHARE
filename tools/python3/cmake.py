from tools.python3 import run


def version():
    pass


def make_config_str(
    path, samrai_dir=None, cxx_flags=None, use_ninja=False, use_ccache=False, extra=""
):
    """
    FULL UNPORTABLE OPTIMIZATIONS = cxx_flags="-O3 -march=native -mtune=native"
    """

    samrai_dir = "" if samrai_dir is None else f"-DSAMRAI_ROOT={samrai_dir}"
    cxx_flags = "" if cxx_flags is None else f'-DCMAKE_CXX_FLAGS="{cxx_flags}"'
    ccache = "" if use_ccache is False else "-DCMAKE_CXX_COMPILER_LAUNCHER=ccache"
    ninja = "" if not use_ninja else "-G Ninja"

    return f"cmake {path} {samrai_dir} {cxx_flags} {ninja} {ccache} {extra}"


def config(
    path, samrai_dir=None, cxx_flags=None, use_ninja=False, use_ccache=False, extra=""
):
    cmd = make_config_str(path, samrai_dir, cxx_flags, use_ninja, extra)
    run(cmd, capture_output=False)


def build(use_ninja=False, threads=1):
    run("ninja" if use_ninja else f"make -j{threads}", capture_output=False)

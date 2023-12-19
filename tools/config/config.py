import json
import subprocess
from pathlib import Path
from dataclasses import dataclass, field
import dataclasses

FILE_DIR = Path(__file__).resolve().parent
BUILD_DIR = FILE_DIR / "build"
ROOT_DIR = FILE_DIR.parent.parent
FILES = [f for f in BUILD_DIR.glob("PHARE_*")]
DEF_DIR = ROOT_DIR / "src" / "core" / "def"
GENERATED_CONFIGS = dict(
    mpi=DEF_DIR / "_gen_mpi.hpp",
    system=DEF_DIR / "_gen_sys.hpp",
)


@dataclass
class SystemSettings:
    cmake_binary: str  # field(default_factory=lambda: "Couldn't parse cmake binary path")
    cmake_version: str  # field(default_factory=lambda: "Couldn't parse cmake version")
    hdf5_version: str  # field(default_factory=lambda: "Couldn't parse cmake hdf5 version")
    mpi_version: str  # field(default_factory=lambda: "Couldn't parse cmake mpi version")
    uname: str  # field(default_factory=lambda: "Couldn't parse uname")


system_cpp_ = """
#include <string_view>

namespace PHARE {{
struct SystemConfig {{

    constexpr static std::string_view UNAME = R"({})";
    constexpr static std::string_view MPI_VERSION = R"({})";
    constexpr static std::string_view HDF5_VERSION = R"({})";
    constexpr static std::string_view CMAKE_VERSION = R"({})";
    constexpr static std::string_view CMAKE_BINARY = R"({})";

}};

}}

"""


def exec(cmd):
    proc = subprocess.run(cmd.split(" "), check=False, capture_output=True)
    return (proc.stdout + proc.stderr).decode().strip()


def file_string_or(filename, fail=None):
    if not fail:
        fail = f"couldn't open file {filename}"
    filepath = BUILD_DIR / filename
    if not filepath.exists():
        return fail
    with open(filepath) as f:
        return f.readline().strip()


def cmake_cache_line_to_key_value_pair(line):
    bits = line.split("=")
    return bits[0], "=".join(bits[1:])


def parse_cmake_cache_file():
    cmaka_cache_file = BUILD_DIR / "CMakeCache.txt"
    dic = {}
    with open(cmaka_cache_file) as f:
        for line in f.readlines():
            line = line.strip()
            if any([line.startswith("#"), "=" not in line]):
                continue
            key, val = cmake_cache_line_to_key_value_pair(line)
            dic[key] = val
    return dic


def cmake_version_from(cmake_cache: dict) -> str:
    if all(
        [
            "CMAKE_CACHE_MAJOR_VERSION:INTERNAL" in cmake_cache,
            "CMAKE_CACHE_MINOR_VERSION:INTERNAL" in cmake_cache,
            "CMAKE_CACHE_PATCH_VERSION:INTERNAL" in cmake_cache,
        ]
    ):
        major = cmake_cache["CMAKE_CACHE_MAJOR_VERSION:INTERNAL"]
        minor = cmake_cache["CMAKE_CACHE_MINOR_VERSION:INTERNAL"]
        patch = cmake_cache["CMAKE_CACHE_PATCH_VERSION:INTERNAL"]
        return f"{major}.{minor}.{patch}"

    return "Failed to find version in cache file"


def cmake_binary_path_from(cmake_cache: dict) -> str:
    """With this we should be able to guarantee future cmake
    calls by us if needed for whatever reason"""
    path_key = "CMAKE_COMMAND:INTERNAL"
    if path_key in cmake_cache:
        return cmake_cache[path_key]
    return "Failed to find cmake binary path in cache file"


def local_file_format(filename):
    return filename[6:-4]


def create_file(f):
    with open(f, "a") as file:
        file.truncate(0)


def ifndef_define(var, val, buf, comment=None):
    buf += f"#ifndef {var} \n"
    if comment:
        buf += f"// {comment} \n"
    buf += f"#define {var} {val}\n"
    buf += f"#endif /* {var} */ \n\n"
    return buf


def config_mpi_version(txtfile, h_file):
    with open(txtfile) as f:
        version = f.readline()
    buf = "\n"
    if any([s in version for s in ["Open MPI", "OpenMPI"]]):
        buf = ifndef_define(
            "OMPI_SKIP_MPICXX", 1, buf, "avoids default including mpicxx"
        )
    if any([s in version for s in ["MPICH"]]):
        buf = ifndef_define(
            "MPICH_SKIP_MPICXX", 1, buf, "avoids default including mpicxx"
        )
    with open(h_file, "w") as f:
        f.write(buf)


def gen_system_file():
    out_file = GENERATED_CONFIGS["system"]
    cmake_cache: dict = parse_cmake_cache_file()

    settings = SystemSettings(
        cmake_binary=cmake_binary_path_from(cmake_cache),
        cmake_version=cmake_version_from(cmake_cache),
        hdf5_version=file_string_or("PHARE_HDF5_version.txt", "HDF5 version failed"),
        mpi_version=file_string_or(
            "PHARE_MPI_Get_library_version.txt", "MPI version failed"
        ),
        uname=exec("uname -a"),
    )
    with open(ROOT_DIR/".phare_config.json", "w") as f:
        json.dump(dataclasses.asdict(settings), f)

    with open(out_file, "w") as f:
        f.write(
            system_cpp_.format(
                settings.uname,
                settings.hdf5_version,
                settings.mpi_version,
                settings.cmake_version,
                settings.cmake_binary,
            )
        )


def config_mpi():
    h_file = GENERATED_CONFIGS["mpi"]
    create_file(h_file)
    operators = dict(MPI_Get_library_version=config_mpi_version)
    for f in FILES:
        local_name = local_file_format(f.name)
        if local_name in operators:
            operators[local_name](f, h_file)


def main():
    config_mpi()
    gen_system_file()


if __name__ == "__main__":
    main()
    print("PHARE configurator generating configs:")
    for k, v in GENERATED_CONFIGS.items():
        print(f"  {v.relative_to(ROOT_DIR)}")

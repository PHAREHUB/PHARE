import os
import json
import subprocess
import dataclasses
from pathlib import Path

FILE_DIR = Path(__file__).resolve().parent
BUILD_DIR = FILE_DIR / "build"
ROOT_DIR = FILE_DIR.parent.parent
DOT_PHARE_DIR = ROOT_DIR / ".phare"
FILES = list(BUILD_DIR.glob("PHARE_*"))
DEF_DIR = ROOT_DIR / "src" / "core" / "def"
GENERATED_CONFIGS = dict(
    mpi=DEF_DIR / "_gen_mpi.hpp",
    system=DEF_DIR / "_gen_sys.hpp",
)


@dataclasses.dataclass
class SystemSettings:
    cmake_binary: str
    cmake_version: str
    cxx_compiler: str
    cxx_compiler_version: str
    hdf5_version: str
    hdf5_is_parallel: str
    mpi_version: str
    python_binary: str
    python_version: str
    uname: str


SYSTEM_CPP_ = """
#include <string_view>
#include <unordered_map>

namespace PHARE {{
struct SystemConfig {{

    constexpr static std::string_view CMAKE_BINARY = R"({})";
    constexpr static std::string_view CMAKE_VERSION = R"({})";
    constexpr static std::string_view CXX_COMPILER = R"({})";
    constexpr static std::string_view CXX_COMPILER_VERSION = R"({})";
    constexpr static std::string_view HDF5_VERSION = R"({})";
    constexpr static std::string_view HDF5_IS_PARALLEL = R"({})";
    constexpr static std::string_view _MPI_VERSION_ = R"({})";
    constexpr static std::string_view PYTHON_BINARY = R"({})";
    constexpr static std::string_view PYTHON_VERSION = R"({})";
    constexpr static std::string_view UNAME = R"({})";

}};

std::unordered_map<std::string, std::string> build_config(){{
  return {{
      {{"CMAKE_BINARY", std::string{{SystemConfig::CMAKE_BINARY}}}},
      {{"CMAKE_VERSION", std::string{{SystemConfig::CMAKE_VERSION}}}},
      {{"CXX_COMPILER", std::string{{SystemConfig::CXX_COMPILER}}}},
      {{"CXX_COMPILER_VERSION", std::string{{SystemConfig::CXX_COMPILER_VERSION}}}},
      {{"HDF5_VERSION", std::string{{SystemConfig::HDF5_VERSION}}}},
      {{"HDF5_IS_PARALLEL", std::string{{SystemConfig::HDF5_IS_PARALLEL}}}},
      {{"MPI_VERSION", std::string{{SystemConfig::_MPI_VERSION_}}}},
      {{"PYTHON_BINARY", std::string{{SystemConfig::PYTHON_BINARY}}}},
      {{"PYTHON_VERSION", std::string{{SystemConfig::PYTHON_VERSION}}}},
      {{"UNAME", std::string{{SystemConfig::UNAME}}}}
  }};
}}

}}

"""


def subprocess_run(cmd, on_error="error message"):
    try:
        proc = subprocess.run(cmd.split(" "), check=True, capture_output=True)
        return (proc.stdout + proc.stderr).decode().strip()
    except Exception:
        return on_error


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


def deduce_mpi_type_from(version):
    mpi_type = "unknown"
    if any([s in version for s in ["Open MPI", "OpenMPI"]]):
        return "OMPI"
    if any([s in version for s in ["MPICH"]]):
        return "MPICH"
    return "unknown"


def config_mpi_version(txtfile, h_file):
    with open(txtfile) as f:
        version = f.readline()
    buf = "\n"
    mpi_type = deduce_mpi_type_from(version)
    if mpi_type == "OMPI":
        buf = ifndef_define(
            "OMPI_SKIP_MPICXX", 1, buf, "avoids default including mpicxx"
        )
    elif mpi_type == "MPICH":
        buf = ifndef_define(
            "MPICH_SKIP_MPICXX", 1, buf, "avoids default including mpicxx"
        )
    with open(h_file, "w") as f:
        f.write(buf)
    return version


def try_compiler_dash_v(compiler):
    return subprocess_run(f"{compiler} -v", "compiler does not support '-v' argument")


def get_python_version(py3_binary):
    return subprocess_run(
        f"{py3_binary} -V", "python3 interpreter does not support '-V' argument"
    )


def gen_system_file():
    out_file = GENERATED_CONFIGS["system"]
    cmake_cache: dict = parse_cmake_cache_file()

    settings = SystemSettings(
        cmake_binary=os.environ["CMAKE_COMMAND"],
        cmake_version=cmake_version_from(cmake_cache),
        cxx_compiler=os.environ["CMAKE_CXX_COMPILER"],
        cxx_compiler_version=try_compiler_dash_v(os.environ["CMAKE_CXX_COMPILER"]),
        hdf5_version=file_string_or("PHARE_HDF5_version.txt", "HDF5 version failed"),
        hdf5_is_parallel=file_string_or(
            "PHARE_HDF5_is_parallel.txt", "HDF5 is parallel check failed"
        ),
        mpi_version=file_string_or(
            "PHARE_MPI_Get_library_version.txt", "MPI version failed"
        ),
        python_binary=os.environ["PYTHON_EXECUTABLE"],
        python_version=get_python_version(os.environ["PYTHON_EXECUTABLE"]),
        uname=subprocess_run("uname -a"),
    )
    with open(DOT_PHARE_DIR / "build_config.json", "w") as f:
        json.dump(dataclasses.asdict(settings), f)

    with open(out_file, "w") as f:
        f.write(
            SYSTEM_CPP_.format(
                settings.cmake_binary,
                settings.cmake_version,
                settings.cxx_compiler,
                settings.cxx_compiler_version,
                settings.hdf5_version,
                settings.hdf5_is_parallel,
                settings.mpi_version,
                settings.python_binary,
                settings.python_version,
                settings.uname,
            )
        )


def config_mpi():
    h_file = GENERATED_CONFIGS["mpi"]
    create_file(h_file)
    operators = dict(MPI_Get_library_version=config_mpi_version)
    results = {}
    for f in FILES:
        local_name = local_file_format(f.name)
        if local_name in operators:
            results[local_name] = operators[local_name](f, h_file)
    return results


def write_local_cmake_file(mpi_results):
    with open("local.cmake", "w") as file:
        file.write(
            """cmake_minimum_required (VERSION 3.20.1)
project(configured_phare)
message("")
message("!!PHARE CONFIGURATED!!")
message("")

"""
        )

        mpi_type = deduce_mpi_type_from(mpi_results["MPI_Get_library_version"])
        if mpi_type == "OMPI":
            # work around for https://github.com/open-mpi/ompi/issues/10761#issuecomment-1236909802
            file.write(
                """set (PHARE_MPIRUN_POSTFIX "${PHARE_MPIRUN_POSTFIX} -q --bind-to none")
                """
            )


def main():
    DOT_PHARE_DIR.mkdir(exist_ok=True, parents=True)
    mpi_results = config_mpi()
    gen_system_file()
    write_local_cmake_file(mpi_results)


if __name__ == "__main__":
    main()
    print("PHARE configurator generating configs:")
    for k, v in GENERATED_CONFIGS.items():
        print(f"  {v.relative_to(ROOT_DIR)}")

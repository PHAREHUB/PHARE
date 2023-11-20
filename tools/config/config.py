import os
from pathlib import Path

FILE_DIR = Path(__file__).resolve().parent
ROOT_DIR = FILE_DIR.parent.parent
FILES = [f for f in Path(FILE_DIR).glob("build/PHARE_*")]

GENERATED_CONFIGS = dict(mpi=ROOT_DIR / "src" / "core" / "def" / "_gen_mpi.hpp")


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


if __name__ == "__main__":
    main()
    print("PHARE configurator generating configs:")
    for k, v in GENERATED_CONFIGS.items():
        print(f"  {v.relative_to(ROOT_DIR)}")

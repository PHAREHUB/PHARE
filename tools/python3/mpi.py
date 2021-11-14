from tools.python3 import run


def mpirun(cmd, procs, capture_output=False, **kwargs):
    return run(f"mpirun -n {procs} {cmd}", capture_output=capture_output, **kwargs)

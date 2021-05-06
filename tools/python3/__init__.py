def decode_bytes(input):
    return input.decode("ascii", errors="ignore")


def run(cmd, shell=True, capture_output=True, check=False, print_cmd=True, **kwargs):
    """
    https://docs.python.org/3/library/subprocess.html
    """
    import subprocess

    if print_cmd:
        print(f"running: {cmd}")
    try:
        return subprocess.run(cmd, shell=shell, capture_output=capture_output, check=check, **kwargs)
    except subprocess.CalledProcessError as e: # only triggers on failure if check=True
        print(f"run failed with error: {e}\n\t{e.stdout}\n\t{e.stderr} ")
        raise RuntimeError(decode_bytes(e.stderr))

def run_mp(cmds):
    """
    spawns len(cmds) threads running all commands at once and waiting for results
    https://docs.python.org/3/library/concurrent.futures.html
    """
    import concurrent.futures

    with concurrent.futures.ThreadPoolExecutor(max_workers=len(cmds)) as executor:
        jobs = [executor.submit(run, cmd) for cmd in cmds]
        results = []
        for future in concurrent.futures.as_completed(jobs):
            try:
                results += [future.result()]
            except Exception as exc:
                print("run_mp generated an exception: %s" % exc)
        return results


def binary_exists_on_path(bin):
    """
    https://linux.die.net/man/1/which
    """
    return run(f"which {bin}").returncode == 0


def scan_dir(path, files_only=False, dirs_only=False, drop=[]):
    import os

    assert os.path.exists(path)
    checks = [
        lambda entry: not files_only or (files_only and entry.is_file()),
        lambda entry: not dirs_only or (dirs_only and entry.is_dir()),
        lambda entry: entry.name not in drop,
    ]
    return [
        entry.name
        for entry in os.scandir(path)
        if all([check(entry) for check in checks])
    ]

import contextlib
@contextlib.contextmanager
def pushd(new_cwd):
    import os
    cwd = os.getcwd()
    os.chdir(new_cwd)
    try:
        yield
    finally:
        os.chdir(cwd)

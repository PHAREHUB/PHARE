def decode_bytes(input):
    return input.decode("ascii", errors="ignore")

class RunTimer:
    def __init__(self, cmd, shell=True, capture_output=True, check=False, print_cmd=True, **kwargs):
        import time
        import subprocess
        self.cmd = cmd
        start = time.time()
        self.run = subprocess.run(self.cmd, shell=shell, capture_output=capture_output, check=check, **kwargs)
        self.t = time.time() - start
        self.stdout = self.run.stdout
        self.stderr = self.run.stderr



def run(cmd, shell=True, capture_output=True, check=False, print_cmd=True, **kwargs):
    """
    https://docs.python.org/3/library/subprocess.html
    """
    import subprocess

    if print_cmd:
        print(f"running: {cmd}")
    try:
        return RunTimer(cmd, shell=shell, capture_output=capture_output, check=check, **kwargs)
    except subprocess.CalledProcessError as e: # only triggers on failure if check=True
        what = f"run failed with error: {e}"
        print(what)
        if capture_output:
            raise RuntimeError(decode_bytes(e.stderr))
        raise RuntimeError(what)

def run_mp(cmds, N_CORES=None, **kwargs):
    """
    spawns N_CORES threads (default=len(cmds)) running commands and waiting for results
    https://docs.python.org/3/library/concurrent.futures.html
    """
    import concurrent.futures

    if N_CORES is None:
        N_CORES = len(cmds)

    with concurrent.futures.ThreadPoolExecutor(max_workers=N_CORES) as executor:
        jobs = [executor.submit(run, cmd, **kwargs) for cmd in cmds]
        results = []
        for future in concurrent.futures.as_completed(jobs):
            try:
                proc = future.result()
                results += [proc]
                if future.exception() is not None:
                    raise future.exception()
                print(proc.cmd, f"finished in {proc.t:.2f} seconds")
            except Exception as exc:
                if kwargs.get("check", False):
                    executor.shutdown(wait=False, cancel_futures=True)
                    raise exc
                else:
                    print(f"run_mp generated an exception: {exc}")
        return results


def find_on_path(file):
    import os
    for dir in os.environ["PATH"].split(os.pathsep):
        full = os.path.join(dir, file)
        if os.path.exists(full):
            return full
    return ""


def binary_exists_on_path(bin):
    """
    https://linux.die.net/man/1/which
    """
    return len(find_on_path(bin))


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


from tools.python3 import decode_bytes, run
import subprocess


def current_branch():
    return decode_bytes(run("git branch --show-current").stdout)


def hashes(N=1):
    return decode_bytes(run(f"git log -{N} --pretty=format:%h").stdout).splitlines()


def log(N=10, use_short=False):
    """ enjoy """
    short = '--pretty=format:"%h%x09%an%x09%ad%x09%s"' if use_short else ""
    return decode_bytes(run(f"git log -{N} {short}").stdout)


def branch_exists(branch):
    try:
        run(f"git show-branch {branch}", check=True)
    except subprocess.CalledProcessError as e:
        return False  # exit failure means branch does not exist
    return True


def checkout(branch, create=False, recreate=False):
    if recreate:
        delete_branch(branch)
        create = True

    if create and not branch_exists(branch):
        run(f"git checkout -b {branch}")
    else:
        run(f"git checkout {branch}")

import fcntl
import os
import sys
import time


def acquire(slots, lock_dir):
    os.makedirs(lock_dir, exist_ok=True)
    fds = [
        os.open(os.path.join(lock_dir, f"slot.{i}"), os.O_RDWR | os.O_CREAT, 0o644)
        for i in range(slots)
    ]
    while True:
        for fd in fds:
            try:
                fcntl.flock(fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
            except OSError:
                continue
            os.set_inheritable(fd, True)
            for other in fds:
                if other != fd:
                    os.close(other)
            return fd
        time.sleep(0.05)


def main(argv):
    if len(argv) < 4:
        print("usage: compile_slot.py <slots> <lock_dir> <command> [args...]", file=sys.stderr)
        return 2
    try:
        slots = int(argv[1])
    except ValueError:
        slots = 0
    if slots < 1:
        print("compile_slot: slots must be an integer >= 1", file=sys.stderr)
        return 2
    acquire(slots, argv[2])
    command = argv[3]
    os.execvp(command, [command, *argv[4:]])


if __name__ == "__main__":
    sys.exit(main(sys.argv))

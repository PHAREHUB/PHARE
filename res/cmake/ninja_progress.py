import os
import re
import select
import signal
import subprocess
import sys
import time

STATUS = re.compile(rb"\[(\d+)/(\d+)\]")
PIP = re.compile(r"pip(\d+(\.\d+)?)?")
VERBOSE = re.compile(r"--verbose|-[A-Za-z]*v[A-Za-z]*")
OFF = ("", "0", "false", "no", "off")
LABEL = "Compiling PHARE "
BAR = 30
REDRAW = 0.1
REFRESH = 1.0


def parse_status(line):
    match = STATUS.match(line)
    if match is None or int(match.group(2)) == 0:
        return None
    return int(match.group(1)), int(match.group(2))


def format_duration(seconds):
    hours, rest = divmod(int(seconds), 3600)
    minutes, seconds = divmod(rest, 60)
    if hours:
        return f"{hours}h{minutes:02d}m"
    if minutes:
        return f"{minutes}m{seconds:02d}s"
    return f"{seconds}s"


def render(done, total, elapsed, width, color=False):
    text = f" {done}/{total} compiled · {format_duration(elapsed)}"
    bar = min(BAR, width - len(LABEL) - len(text))
    if bar < 5:
        return (LABEL + text.lstrip())[: max(width, 0)]
    filled = bar * done // total
    label = f"\x1b[1;36m{LABEL}\x1b[0m" if color else LABEL
    return label + "█" * filled + "░" * (bar - filled) + text


def frontend_verbosity(argv, environ):
    names = [os.path.basename(arg) for arg in argv[:2]]
    if not names:
        return None
    if names[0] == "uv":
        return any(VERBOSE.fullmatch(arg) for arg in argv[1:])
    is_pip = PIP.fullmatch(names[0]) or (
        names[0].startswith("python")
        and (argv[1:3] == ["-m", "pip"] or (len(names) > 1 and PIP.fullmatch(names[1])))
    )
    if not is_pip:
        return None
    if environ.get("PIP_VERBOSE", "").lower() not in OFF:
        return True
    return any(VERBOSE.fullmatch(arg) for arg in argv[1:])


def frontend_is_verbose(chain, environ):
    for argv in chain:
        verbose = frontend_verbosity(argv, environ)
        if verbose is not None:
            return verbose
    return False


def _proc_parent(pid, proc_root):
    try:
        with open(f"{proc_root}/{pid}/cmdline", "rb") as f:
            argv = [
                arg.decode(errors="replace") for arg in f.read().split(b"\0") if arg
            ]
        with open(f"{proc_root}/{pid}/stat") as f:
            ppid = int(f.read().rsplit(")", 1)[1].split()[1])
    except (OSError, IndexError, ValueError):
        return None
    return argv, ppid


def _ps_parent(pid):
    try:
        out = subprocess.run(
            ["ps", "-ww", "-o", "ppid=", "-o", "args=", "-p", str(pid)],
            capture_output=True,
            text=True,
            check=False,
        ).stdout.split()
    except OSError:
        return None
    if not out or not out[0].isdigit():
        return None
    return out[1:], int(out[0])


def process_chain(pid, proc_root="/proc"):
    use_proc = os.path.isdir(proc_root)
    for _ in range(64):
        if pid <= 1:
            return
        found = _proc_parent(pid, proc_root) if use_proc else _ps_parent(pid)
        if found is None:
            return
        argv, pid = found
        yield argv


def targets_build_dir(build_dir, args, cwd):
    directory = cwd
    for i, arg in enumerate(args):
        if arg.startswith("-t"):
            return False
        if arg == "-C" and i + 1 < len(args):
            directory = os.path.join(directory, args[i + 1])
        elif arg.startswith("-C"):
            directory = os.path.join(directory, arg[2:])
    return os.path.realpath(directory) == os.path.realpath(build_dir)


def open_terminal(environ):
    if environ.get("TERM") == "dumb":
        return None
    try:
        fd = os.open("/dev/tty", os.O_WRONLY | os.O_NOCTTY)
    except OSError:
        return None
    if not os.isatty(fd):
        os.close(fd)
        return None
    return fd


class TerminalProgress:
    def __init__(self, fd, color):
        self.fd = fd
        self.color = color
        self.started = False
        self.last = 0.0
        self.shown = None

    def _write(self, data):
        try:
            os.write(self.fd, data)
        except OSError:
            pass

    def _width(self):
        try:
            return os.get_terminal_size(self.fd).columns or 80
        except OSError:
            return 80

    def draw(self, done, total, elapsed, now):
        final = done == total and self.shown != (done, total)
        if not final and now - self.last < REDRAW:
            return
        self.last = now
        self.shown = (done, total)
        percent = 100 * done // total
        line = render(done, total, elapsed, self._width() - 1, self.color)
        head = b"" if self.started else b"\x1b[22;0t\x1b[?25l"
        self.started = True
        self._write(
            head
            + b"\x1bD\x1bM\x1b7\x1bD\r"
            + line.encode()
            + b"\x1b[J\x1b8"
            + f"\x1b]0;PHARE build {percent}%\x07\x1b]9;4;1;{percent}\x07".encode()
        )

    def close(self):
        if self.started:
            self._write(b"\x1b7\x1bD\r\x1b[J\x1b8\x1b]9;4;0;0\x07\x1b[23;0t\x1b[?25h")
        os.close(self.fd)


def relay(command, progress):
    proc = subprocess.Popen(command, stdout=subprocess.PIPE)
    signal.signal(signal.SIGINT, signal.SIG_IGN)
    source = proc.stdout.fileno()
    out = sys.stdout.buffer
    start = time.monotonic()
    pending = b""
    status = None
    try:
        while True:
            ready, _, _ = select.select([source], [], [], REFRESH)
            if ready:
                data = os.read(source, 65536)
                if not data:
                    break
                out.write(data)
                out.flush()
                *lines, pending = (pending + data).split(b"\n")
                for line in lines:
                    status = parse_status(line) or status
            if status is not None:
                now = time.monotonic()
                progress.draw(*status, now - start, now)
        code = proc.wait()
    finally:
        progress.close()
    return 128 - code if code < 0 else code


def main(argv):
    build_dir, ninja, args = argv[1], argv[2], argv[3:]
    environ = os.environ
    fd = None
    if (
        environ.get("PHARE_NO_PROGRESS", "").lower() in OFF
        and not sys.stdout.isatty()
        and targets_build_dir(build_dir, args, os.getcwd())
        and not frontend_is_verbose(process_chain(os.getppid()), environ)
    ):
        fd = open_terminal(environ)
    if fd is None:
        os.execv(ninja, [ninja, *args])
    return relay([ninja, *args], TerminalProgress(fd, not environ.get("NO_COLOR")))


if __name__ == "__main__":
    sys.exit(main(sys.argv))

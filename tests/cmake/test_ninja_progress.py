import fcntl
import importlib.util
import os
import pty
import re
import shutil
import signal
import struct
import subprocess
import sys
import tempfile
import termios
import threading
import time
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "res/cmake/ninja_progress.py"
MODULE = ROOT / "res/cmake/ninja_progress.cmake"
CMAKE = os.environ.get("CMAKE_COMMAND") or "cmake"
TIMEOUT = 60
CLEAR = b"\x1b7\x1bD\r\x1b[J\x1b8"
CLOSE = CLEAR + b"\x1b]9;4;0;0\x07\x1b[23;0t\x1b[?25h"
FAKE_UV = (
    "import fcntl, signal, subprocess, sys, termios\n"
    "fcntl.ioctl(int(sys.argv[1]), termios.TIOCSCTTY, 0)\n"
    "signal.signal(signal.SIGINT, lambda *_: None)\n"
    "sys.exit(subprocess.call(sys.argv[sys.argv.index('--') + 1 :]))\n"
)


def load_script():
    spec = importlib.util.spec_from_file_location("ninja_progress", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def fake_ninja(directory, body):
    path = Path(directory) / "ninja"
    path.write_text(f"#!{sys.executable}\nimport sys, time\n{body}\n")
    path.chmod(0o755)
    return path


def emit(*lines, pause=0.3):
    return "\n".join(
        f"sys.stdout.write({line!r}); sys.stdout.flush(); time.sleep({pause})"
        for line in lines
    )


def clean_env(**extra):
    env = {
        k: v
        for k, v in os.environ.items()
        if k not in ("PHARE_NO_PROGRESS", "NO_COLOR", "PIP_VERBOSE")
    }
    env["TERM"] = "xterm-256color"
    env.update(extra)
    return env


def run_on_terminal(
    build_dir,
    ninja,
    ninja_args=(),
    cwd=None,
    env=None,
    uv_flags=(),
    interrupt_on=None,
    output_on_terminal=False,
):
    master, slave = pty.openpty()
    fcntl.ioctl(slave, termios.TIOCSWINSZ, struct.pack("HHHH", 24, 100, 0, 0))
    proc = subprocess.Popen(
        ["uv", "-c", FAKE_UV, str(slave), "sync", *uv_flags, "--", sys.executable]
        + [str(SCRIPT), str(build_dir), str(ninja), *ninja_args],
        executable=sys.executable,
        cwd=cwd or build_dir,
        env=env or clean_env(),
        stdin=subprocess.DEVNULL,
        stdout=slave if output_on_terminal else subprocess.PIPE,
        stderr=subprocess.PIPE,
        start_new_session=True,
        pass_fds=(slave,),
    )
    os.close(slave)
    screen = bytearray()

    def read_terminal():
        while True:
            try:
                data = os.read(master, 4096)
            except OSError:
                return
            if not data:
                return
            screen.extend(data)

    reader = threading.Thread(target=read_terminal, daemon=True)
    reader.start()
    try:
        if interrupt_on is not None:
            deadline = time.monotonic() + TIMEOUT
            while interrupt_on not in screen:
                if time.monotonic() > deadline:
                    raise AssertionError(
                        f"never drew {interrupt_on!r}: {bytes(screen)!r}"
                    )
                time.sleep(0.02)
            os.killpg(proc.pid, signal.SIGINT)
        out, _ = proc.communicate(timeout=TIMEOUT)
        reader.join(TIMEOUT)
    finally:
        if proc.poll() is None:
            os.killpg(proc.pid, signal.SIGKILL)
            proc.communicate()
        os.close(master)
    return bytes(screen), out, proc.returncode


class FormatTest(unittest.TestCase):
    def setUp(self):
        self.script = load_script()

    def test_parses_ninja_status(self):
        parse = self.script.parse_status
        self.assertEqual(parse(b"[12/345] Building CXX object a.o\n"), (12, 345))
        self.assertIsNone(parse(b"FAILED: a.o\n"))
        self.assertIsNone(parse(b"  [1/2] indented\n"))
        self.assertIsNone(parse(b"[0/0] nothing\n"))

    def test_formats_durations(self):
        fmt = self.script.format_duration
        self.assertEqual(fmt(12.9), "12s")
        self.assertEqual(fmt(372), "6m12s")
        self.assertEqual(fmt(3725), "1h02m")

    def test_wide_terminal_shows_bar(self):
        line = self.script.render(412, 1730, 372, 99)
        self.assertEqual(
            line, "Compiling PHARE " + "█" * 7 + "░" * 23 + " 412/1730 compiled · 6m12s"
        )

    def test_no_remaining_time_estimate(self):
        for done in (1, 87, 412, 1729, 1730):
            self.assertTrue(
                self.script.render(done, 1730, 30, 99).endswith("compiled · 30s")
            )

    def test_narrow_terminal_drops_bar(self):
        line = self.script.render(412, 1730, 372, 45)
        self.assertEqual(line, "Compiling PHARE 412/1730 compiled · 6m12s")

    def test_never_wider_than_terminal(self):
        for width in range(120):
            self.assertLessEqual(len(self.script.render(412, 1730, 372, width)), width)

    def test_color_only_wraps_label(self):
        plain = self.script.render(412, 1730, 372, 99)
        colored = self.script.render(412, 1730, 372, 99, color=True)
        self.assertTrue(colored.startswith("\x1b[1;36mCompiling PHARE \x1b[0m"))
        self.assertEqual(re.sub(r"\x1b\[[0-9;]*m", "", colored), plain)


class FrontendTest(unittest.TestCase):
    def setUp(self):
        self.script = load_script()

    def verbosity(self, *argv, env=None):
        return self.script.frontend_verbosity(list(argv), env or {})

    def test_uv(self):
        self.assertTrue(self.verbosity("/usr/bin/uv", "sync", "-v"))
        self.assertTrue(self.verbosity("uv", "pip", "install", "-vv", "."))
        self.assertTrue(self.verbosity("uv", "sync", "--verbose"))
        self.assertFalse(self.verbosity("uv", "sync", "--reinstall-package", "pyphare"))

    def test_pip(self):
        self.assertTrue(
            self.verbosity("/v/bin/python3", "/v/bin/pip", "install", "-v", ".")
        )
        self.assertTrue(
            self.verbosity("python3", "-m", "pip", "install", "--verbose", ".")
        )
        self.assertTrue(self.verbosity("pip3.12", "install", "-Uv", "."))
        self.assertFalse(self.verbosity("pip", "install", "-V"))
        self.assertFalse(self.verbosity("pip", "install", "."))

    def test_pip_verbose_environment(self):
        self.assertTrue(self.verbosity("pip", "install", ".", env={"PIP_VERBOSE": "1"}))
        self.assertFalse(
            self.verbosity("pip", "install", ".", env={"PIP_VERBOSE": "0"})
        )
        self.assertFalse(self.verbosity("uv", "sync", env={"PIP_VERBOSE": "1"}))

    def test_other_processes_are_not_frontends(self):
        self.assertIsNone(self.verbosity("cmake", "--build", "build", "-v"))
        self.assertIsNone(
            self.verbosity("/usr/bin/python3", "_in_process.py", "build_wheel")
        )
        self.assertIsNone(self.verbosity("ninja", "-v"))
        self.assertIsNone(self.verbosity())

    def test_closest_frontend_decides(self):
        chain = [
            ["cmake", "--build", "b"],
            ["pip", "install", "."],
            ["uv", "run", "-v"],
        ]
        self.assertFalse(self.script.frontend_is_verbose(chain, {}))
        self.assertTrue(self.script.frontend_is_verbose(chain[2:], {}))
        self.assertFalse(self.script.frontend_is_verbose([["bash"]], {}))


class ProcessChainTest(unittest.TestCase):
    def setUp(self):
        self.script = load_script()

    def write_process(self, root, pid, ppid, argv, comm="x"):
        directory = root / str(pid)
        directory.mkdir()
        (directory / "cmdline").write_bytes(
            b"\0".join(a.encode() for a in argv) + b"\0"
        )
        (directory / "stat").write_text(f"{pid} ({comm}) S {ppid} {pid} {pid} 0 -1\n")

    def test_walks_parents_until_init(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            self.write_process(root, 300, 200, ["make", "x"], comm="a) S 7 (b")
            self.write_process(root, 200, 1, ["uv", "sync", "-v"])
            self.write_process(root, 1, 0, ["init"])
            chain = list(self.script.process_chain(300, str(root)))
        self.assertEqual(chain, [["make", "x"], ["uv", "sync", "-v"]])

    def test_vanished_process_ends_chain(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            self.write_process(root, 300, 200, ["make"])
            self.assertEqual(
                list(self.script.process_chain(300, str(root))), [["make"]]
            )

    @unittest.skipUnless(
        shutil.which("ps") and os.path.isdir("/proc"), "needs ps and /proc"
    )
    def test_ps_fallback_matches_proc(self):
        with tempfile.TemporaryDirectory() as tmp:
            missing = str(Path(tmp) / "missing")
            from_ps = list(self.script.process_chain(os.getpid(), missing))
        from_proc = list(self.script.process_chain(os.getpid()))
        self.assertEqual(from_ps[0], from_proc[0])
        self.assertEqual(len(from_ps), len(from_proc))


class BuildDirectoryTest(unittest.TestCase):
    def setUp(self):
        self.script = load_script()
        self.tmp = tempfile.TemporaryDirectory()
        self.build = Path(self.tmp.name) / "build"
        (self.build / "CMakeFiles/TryCompile-x").mkdir(parents=True)

    def tearDown(self):
        self.tmp.cleanup()

    def targets(self, args, cwd=None):
        return self.script.targets_build_dir(
            str(self.build), list(args), str(cwd or self.build)
        )

    def test_main_build(self):
        self.assertTrue(self.targets([]))
        self.assertTrue(self.targets(["-j", "8", "all"]))
        self.assertTrue(self.targets(["-C", "."]))
        self.assertTrue(self.targets(["-C", str(self.build)], cwd=self.tmp.name))
        self.assertTrue(self.targets([f"-C{self.build}"], cwd=self.tmp.name))

    def test_other_directories(self):
        self.assertFalse(self.targets([], cwd=self.build / "CMakeFiles/TryCompile-x"))
        self.assertFalse(
            self.targets(["-C", str(self.build / "CMakeFiles/TryCompile-x")])
        )

    def test_tools(self):
        self.assertFalse(self.targets(["-t", "recompact"]))
        self.assertFalse(self.targets(["-trestat"]))


class RelayTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.build = Path(self.tmp.name) / "build"
        self.build.mkdir()

    def tearDown(self):
        self.tmp.cleanup()

    def ninja(self, body):
        return fake_ninja(self.tmp.name, body)

    def test_without_terminal_output_and_exit_code_pass_through(self):
        ninja = self.ninja(
            emit("[1/2] a\n", "noise\n", "[2/2] b\n", "tail") + "\nsys.exit(3)"
        )
        proc = subprocess.run(
            [sys.executable, str(SCRIPT), str(self.build), str(ninja)],
            cwd=self.build,
            env=clean_env(),
            capture_output=True,
            start_new_session=True,
            timeout=TIMEOUT,
            check=False,
        )
        self.assertEqual(proc.stdout, b"[1/2] a\nnoise\n[2/2] b\ntail")
        self.assertEqual(proc.returncode, 3)

    def test_draws_on_terminal(self):
        ninja = self.ninja(
            emit("[1/2] a\n", "noise\n", "[2/2] b\n", "tail") + "\nsys.exit(3)"
        )
        screen, out, code = run_on_terminal(self.build, ninja)
        self.assertEqual(out, b"[1/2] a\nnoise\n[2/2] b\ntail")
        self.assertEqual(code, 3)
        self.assertTrue(screen.startswith(b"\x1b[22;0t\x1b[?25l"), screen)
        self.assertIn(
            b"\x1bD\x1bM\x1b7\x1bD\r\x1b[1;36mCompiling PHARE \x1b[0m", screen
        )
        self.assertIn("1/2 compiled · 0s".encode(), screen)
        self.assertIn("2/2 compiled · 0s".encode(), screen)
        self.assertIn(b"\x1b]0;PHARE build 50%\x07\x1b]9;4;1;50\x07", screen)
        self.assertIn(b"\x1b]0;PHARE build 100%\x07\x1b]9;4;1;100\x07", screen)
        self.assertTrue(screen.endswith(CLOSE), screen)

    def test_last_step_is_drawn_even_right_after_another(self):
        ninja = self.ninja(emit("[1/2] a\n", pause=0.02) + "\n" + emit("[2/2] b\n"))
        screen, _, code = run_on_terminal(self.build, ninja)
        self.assertEqual(code, 0)
        self.assertIn(b"2/2 compiled", screen)

    def test_elapsed_time_advances_while_waiting(self):
        ninja = self.ninja(emit("[1/2] a\n", pause=2.6) + "\n" + emit("[2/2] b\n"))
        screen, _, code = run_on_terminal(self.build, ninja)
        self.assertEqual(code, 0)
        self.assertIn("1/2 compiled · 2s".encode(), screen)

    def test_no_color(self):
        ninja = self.ninja(emit("[1/2] a\n"))
        screen, _, _ = run_on_terminal(self.build, ninja, env=clean_env(NO_COLOR="1"))
        self.assertIn(b"\rCompiling PHARE ", screen)
        self.assertNotIn(b"\x1b[1;36m", screen)

    def test_silent_until_first_status(self):
        ninja = self.ninja(emit("ninja: no work to do.\n"))
        screen, out, code = run_on_terminal(self.build, ninja)
        self.assertEqual((screen, out, code), (b"", b"ninja: no work to do.\n", 0))

    def assertSilent(self, **kwargs):
        ninja = self.ninja(emit("[1/2] a\n"))
        screen, out, code = run_on_terminal(self.build, ninja, **kwargs)
        self.assertEqual((screen, out, code), (b"", b"[1/2] a\n", 0))

    def test_silent_when_frontend_is_verbose(self):
        self.assertSilent(uv_flags=("-v",))

    def test_silent_when_disabled(self):
        self.assertSilent(env=clean_env(PHARE_NO_PROGRESS="1"))

    def test_silent_on_dumb_terminal(self):
        self.assertSilent(env=clean_env(TERM="dumb"))

    def test_silent_for_tools(self):
        self.assertSilent(ninja_args=("-t", "recompact"))

    def test_silent_for_other_build_directories(self):
        other = Path(self.tmp.name) / "CMakeScratch"
        other.mkdir()
        self.assertSilent(cwd=other)

    def test_ninja_keeps_terminal_when_its_output_goes_there(self):
        ninja = self.ninja(
            "print(sys.stdout.isatty(), flush=True)\n" + emit("[1/2] a\n")
        )
        screen, _, code = run_on_terminal(self.build, ninja, output_on_terminal=True)
        self.assertEqual((screen, code), (b"True\r\n[1/2] a\r\n", 0))

    def test_interrupt_clears_the_line(self):
        ninja = self.ninja(emit("[1/3] a\n", pause=TIMEOUT))
        started = time.monotonic()
        screen, _, code = run_on_terminal(
            self.build, ninja, interrupt_on=b"1/3 compiled"
        )
        self.assertLess(time.monotonic() - started, TIMEOUT / 2)
        self.assertEqual(code, 130)
        self.assertTrue(screen.endswith(CLOSE), screen)


@unittest.skipUnless(shutil.which("ninja"), "needs ninja")
class LauncherTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.source = Path(self.tmp.name) / "source"
        self.build = Path(self.tmp.name) / "build"
        self.source.mkdir()
        self.ninja = shutil.which("ninja")
        (self.source / "CMakeLists.txt").write_text(
            "cmake_minimum_required(VERSION 3.20)\n"
            "project(probe NONE)\n"
            f'include("{MODULE}")\n'
            "add_custom_target(hello ALL COMMAND ${CMAKE_COMMAND} -E echo hello-from-build)\n"
        )

    def tearDown(self):
        self.tmp.cleanup()

    def cmake(self, *args):
        proc = subprocess.run(
            [CMAKE, *args],
            capture_output=True,
            text=True,
            timeout=TIMEOUT,
            env=clean_env(),
            start_new_session=True,
            check=False,
        )
        self.assertEqual(proc.returncode, 0, proc.stdout + proc.stderr)
        return proc.stdout

    def configure(self, *defines):
        self.cmake(
            "-S",
            str(self.source),
            "-B",
            str(self.build),
            "-G",
            "Ninja",
            f"-DPHARE_PROJECT_DIR={ROOT}",
            f"-DPython_EXECUTABLE={sys.executable}",
            *defines,
        )

    def cached(self, name):
        text = (self.build / "CMakeCache.txt").read_text()
        match = re.search(rf"^{name}:[A-Z]+=(.*)$", text, re.MULTILINE)
        return match and match.group(1)

    def test_plain_cmake_keeps_ninja(self):
        self.configure(f"-DCMAKE_MAKE_PROGRAM={self.ninja}")
        self.assertEqual(self.cached("CMAKE_MAKE_PROGRAM"), self.ninja)
        self.assertFalse((self.build / "phare_ninja").exists())

    def test_pip_build_goes_through_launcher(self):
        self.configure(f"-DCMAKE_MAKE_PROGRAM={self.ninja}", "-DSKBUILD=2")
        launcher = self.build / "phare_ninja"
        self.assertEqual(self.cached("CMAKE_MAKE_PROGRAM"), str(launcher))
        self.assertEqual(self.cached("PHARE_NINJA"), self.ninja)
        self.assertTrue(os.access(launcher, os.X_OK))
        self.assertIn(f'"{SCRIPT}" "{self.build}" "{self.ninja}"', launcher.read_text())
        self.assertIn("hello-from-build", self.cmake("--build", str(self.build)))

    def test_reconfigure_keeps_real_ninja(self):
        self.configure(f"-DCMAKE_MAKE_PROGRAM={self.ninja}", "-DSKBUILD=2")
        self.configure("-DSKBUILD=2")
        self.assertEqual(self.cached("PHARE_NINJA"), self.ninja)
        self.assertIn("hello-from-build", self.cmake("--build", str(self.build)))


if __name__ == "__main__":
    unittest.main()

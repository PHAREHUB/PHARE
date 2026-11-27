"""include-what-you-use for -DwithIWYU=ON, set as CMAKE_CXX_INCLUDE_WHAT_YOU_USE (res/cmake/iwyu.cmake)

CMake runs this per translation unit as: iwyu.py <options> -- <iwyu args> <compile args>
iwyu's report for the translation unit is kept in <build>/iwyu/, then all kept reports are applied
by fix_include to copies of the reported files, and <build>/iwyu.diff is the diff of those copies,
a patch for `git apply`, empty if there is nothing to change. Translation units outside --only
directories (dependencies, subprojects) are skipped, as are all but the first compile command of
a source in --compile_commands: the others are the same source with other options (e.g. the
cpp_simulator.cpp permutations) and would mostly repeat its report.

Most of src/ is headers no translation unit is named after, iwyu only reports on a translation
unit and its own header. So when the translation unit is --headers_with, every header in the
--only directories is also checked on its own (-x c++-header), --jobs at a time. Headers in the
--only directory of --headers_with take its flags, it has every include path. Other headers (tests,
benches) take the flags, from --compile_commands, of a translation unit that includes them,
directly or through other headers: their targets add include paths (gtest, benchmark, pybind11,
the test's own directories). A header that doesn't compile on its own, or that no compiled
translation unit includes (e.g. tests of a disabled directory), is listed in <build>/iwyu.errors.

Test targets get src/ through -isystem, iwyu would then spell PHARE's headers with <>: -isystem
paths in the --only directories are passed to iwyu as -I.

fix_include doesn't order what it adds as the project does. --clang_format formats the #include
runs of each changed copy with the project's .clang-format (only those lines), so the diff is
empty exactly when the files both include what they use and have their includes formatted.
"""

import argparse
import concurrent.futures
import difflib
import fcntl
import json
import os
import re
import shlex
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path


def main():
    split = sys.argv.index("--")
    args = _parse_args(sys.argv[1:split])
    rest = sys.argv[split + 1 :]

    project = Path(args.project).resolve()
    source = Path(rest[rest.index("-c") + 1]).resolve()
    output = rest[rest.index("-o") + 1]
    if not _is_checked(source, project, args.only):
        return 0
    if not _is_first(source, output, args.compile_commands):
        return 0
    rest = _project_paths_not_system(rest, project, args)

    reports = Path(args.build) / "iwyu"
    reports.mkdir(parents=True, exist_ok=True)

    # the object file is unique per target, a source can be compiled for several
    report = reports / (output.replace("/", "__") + ".out")
    # non-zero exit when it has suggestions
    proc = subprocess.run([args.iwyu, *rest], check=False, capture_output=True, text=True)
    _write(report, proc.stdout + proc.stderr)

    if args.headers_with and source == (project / args.headers_with).resolve():
        headers(args, project, reports, rest)

    # parallel builds: each translation unit wrote its report before taking the lock,
    # so the last one to take it sees every report
    with open(reports / ".lock", "w") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        _write(Path(args.build) / "iwyu.diff", diff(args, project, reports))
    return 0  # never fail the compile, the diff is the result


def headers(args, project, reports, rest):
    """every header in the --only directories through iwyu on its own, with the flags of
    --headers_with or of the nearest translation unit (module docstring)"""
    iwyu_args, own = _split_iwyu_args(rest)
    own_root = _only_root(project / args.headers_with, project, args.only)
    by_source = _flags_by_source(project, args)
    includers = _includers(project, args.only)
    files = sorted(
        f for o in args.only for f in (project / o).rglob("*") if f.suffix in (".h", ".hpp")
    )

    def check(header):
        rel = header.relative_to(project)
        if header.is_relative_to(own_root):
            cwd, flags = None, own
        else:
            found = _includer_flags(header, by_source, includers)
            if not found:
                return f"{rel} (no compiled translation unit includes it)"
            cwd, flags = found
        proc = subprocess.run(  # non-zero exit when it has suggestions
            [args.iwyu, *iwyu_args, *flags, "-x", "c++-header", str(header)],
            check=False,
            capture_output=True,
            text=True,
            cwd=cwd,
        )
        out = proc.stdout + proc.stderr
        _write(reports / ("header__" + str(rel).replace("/", "__") + ".out"), out)
        return str(rel) if "error:" in out else None

    with concurrent.futures.ThreadPoolExecutor(args.jobs) as pool:
        errors = [e for e in pool.map(check, files) if e]
    _write(Path(args.build) / "iwyu.errors", "".join(e + "\n" for e in errors))


def diff(args, project, reports):
    """all reports applied to copies of the reported files, as a patch relative to the project"""
    text = "".join(f.read_text() for f in sorted(reports.glob("*.out")))
    reported = {
        Path(m.group(1)).resolve()
        for m in re.finditer(r"^(.+) should (?:add|remove) these lines:$", text, re.M)
    }
    reported = sorted(f for f in reported if f.is_file() and _is_checked(f, project, args.only))
    if not reported:
        return ""

    with tempfile.TemporaryDirectory(dir=args.build) as tmp:
        for f in reported:
            copy = Path(tmp) / f.relative_to(project)
            copy.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(f, copy)

        # reports name files by absolute path, point them at the copies
        text = text.replace(str(project) + os.sep, tmp + os.sep)
        only_re = "^" + re.escape(tmp) + "/(" + "|".join(map(re.escape, args.only)) + ")/"
        subprocess.run(
            [args.fix_include, *args.fix_args.split(), "--only_re", only_re],
            check=False,  # non-zero exit when there are edits
            input=text,
            capture_output=True,
            text=True,
        )

        patch = []
        for f in reported:
            rel = f.relative_to(project)
            before, after = f.read_text(), (Path(tmp) / rel).read_text()
            if after == before:
                continue
            if args.clang_format:
                after = format_includes(args.clang_format, project, Path(tmp) / rel)
            patch += difflib.unified_diff(
                before.splitlines(keepends=True),
                after.splitlines(keepends=True),
                f"a/{rel}",
                f"b/{rel}",
            )
    return "".join(patch)


def format_includes(clang_format, project, path):
    """clang-format only the #include runs of the file (blank lines inside a run included, and
    any leading blank lines fix_include leaves when it replaces the first include)"""
    lines = path.read_text().splitlines()
    ranges, i = [], 0
    while i < len(lines):
        if not _include_kind(lines[i]):
            i += 1
            continue
        start = end = i
        while i < len(lines) and (_include_kind(lines[i]) or not lines[i].strip()):
            if _include_kind(lines[i]):
                end = i
            i += 1
        ranges.append((start + 1, end + 1))
    if ranges and all(not line.strip() for line in lines[: ranges[0][0] - 1]):
        ranges[0] = (1, ranges[0][1])
    if ranges:
        style = f"--style=file:{project / '.clang-format'}"
        lines_args = [f"--lines={a}:{b}" for a, b in ranges]
        subprocess.run([clang_format, style, *lines_args, "-i", str(path)], check=True)
    return path.read_text()


def _include_kind(line):
    """'<' or '"' for an #include line, else None"""
    m = re.match(r"\s*#\s*include\s*([<\"])", line)
    return m.group(1) if m else None


def _split_iwyu_args(rest):
    """the iwyu options of the command (-Xiwyu <arg> pairs and the driver mode CMake may add),
    and the compile args without output (-o) and source (-c)"""
    iwyu_args, compile_args, i = [], [], 0
    while i < len(rest):
        if rest[i] == "-Xiwyu":
            iwyu_args += rest[i : i + 2]
            i += 2
            continue
        (iwyu_args if rest[i].startswith("--driver-mode=") else compile_args).append(rest[i])
        i += 1
    return iwyu_args, _without_output_and_source(compile_args)


def _without_output_and_source(compile_args):
    compile_args = list(compile_args)
    for opt in ("-o", "-c"):
        if opt in compile_args:
            i = compile_args.index(opt)
            del compile_args[i : i + 2]
    return compile_args


def _project_paths_not_system(compile_args, project, args):
    """-isystem as -I for the --only directories, in the project and the build: iwyu spells what
    it finds through a system path with <>, PHARE's own headers are spelled with """""
    own = [d / o for d in (project, Path(args.build).resolve()) for o in args.only]
    out, i = [], 0
    while i < len(compile_args):
        arg = compile_args[i]
        if arg == "-isystem" and i + 1 < len(compile_args):
            path, step = compile_args[i + 1], 2
        elif arg.startswith("-isystem") and len(arg) > len("-isystem"):
            path, step = arg[len("-isystem") :], 1
        else:
            out.append(arg)
            i += 1
            continue
        if any(Path(path).resolve().is_relative_to(d) for d in own):
            out.append("-I" + path)
        else:
            out += compile_args[i : i + step]
        i += step
    return out


def _flags_by_source(project, args):
    """{source: (working directory, compile args)} of the first compile command of each checked
    translation unit in --compile_commands"""
    try:
        entries = json.loads(Path(args.compile_commands).read_text())
    except (OSError, ValueError):
        return {}
    by_source = {}
    for entry in entries:
        source = Path(entry["directory"], entry["file"]).resolve()
        if source in by_source or not _is_checked(source, project, args.only):
            continue
        command = entry.get("arguments") or shlex.split(entry["command"])
        compile_args = _without_output_and_source(command[1:])  # without the compiler
        by_source[source] = (
            entry["directory"],
            _project_paths_not_system(compile_args, project, args),
        )
    return by_source


def _includers(project, only):
    """{file: files that #include it}, for the files in the --only directories, the include
    resolved next to the including file or from the project root"""
    includers = {}
    for o in only:
        for f in (project / o).rglob("*"):
            if f.suffix not in (".h", ".hpp", ".ipp", ".cpp") or not f.is_file():
                continue
            for name in re.findall(r'^\s*#\s*include\s*[<"]([^>"]+)[>"]', f.read_text(), re.M):
                for candidate in (f.parent / name, project / name):
                    if candidate.is_file():
                        includers.setdefault(candidate.resolve(), set()).add(f.resolve())
                        break
    return includers


def _includer_flags(header, by_source, includers):
    """the flags of the translation unit nearest the header in the include graph that includes
    it, directly or through other headers, preferring one in the header's directory: its target
    has the include paths the header needs (gtest, benchmark, pybind11, test directories)"""
    seen, level = {header}, [header]
    while level:
        found = sorted(f for f in level if f in by_source)
        if found:
            same_dir = [f for f in found if f.parent == header.parent]
            return by_source[(same_dir or found)[0]]
        level = sorted({i for f in level for i in includers.get(f, ()) if i not in seen})
        seen.update(level)
    return None


def _only_root(path, project, only):
    return next((project / o for o in only if path.is_relative_to(project / o)), None)


def _is_first(source, output, compile_commands):
    """whether this is the first compile command of the source, any if there is no database"""
    try:
        entries = json.loads(Path(compile_commands).read_text())
    except (OSError, ValueError):
        return True
    for entry in entries:
        if Path(entry["directory"], entry["file"]).resolve() == source:
            return entry.get("output", output) == output
    return True


def _is_checked(path, project, only):
    return any(path.is_relative_to(project / o) for o in only)


def _parse_args(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument("--iwyu", required=True)
    parser.add_argument("--fix_include", required=True)
    parser.add_argument("--fix_args", default="")
    parser.add_argument("--clang_format", default="")
    parser.add_argument("--project", required=True)
    parser.add_argument("--build", required=True)
    parser.add_argument("--compile_commands", default="")
    parser.add_argument("--headers_with", default="")
    parser.add_argument("--jobs", type=int, default=4)
    parser.add_argument("--only", action="append", required=True)
    return parser.parse_args(argv)


def _write(path, text):
    tmp = path.with_name(path.name + f".{os.getpid()}.tmp")
    tmp.write_text(text)
    os.replace(tmp, path)


if __name__ == "__main__":
    sys.exit(main())

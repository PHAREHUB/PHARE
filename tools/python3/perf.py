
import os
from tools.python3 import run, run_mp


def version():
    # validated on perf version: 5.19
    proc = run("perf -v", shell=True, capture_output=True)
    if " " not in proc or "." not in proc:
        raise ValueError("Unparsable result from 'perf -v'")
    return [int(digit) for digit in proc.split(" ").split(".")]


def check(force_kernel_space = False):
    """ perf can require some system config / read the error if thrown """
    kernel_space_opt = "a" if force_kernel_space else ""
    run(f"perf stat -{kernel_space_opt}d sleep 1", shell=True, capture_output=True, check=True)
    record("ls", [], "/tmp/perf_record_check.dat")


def parse_key(key, force_kernel_space):
    user_space_postfix = ":u"
    if key.endswith(user_space_postfix):
        if force_kernel_space:
            raise RuntimeError(f"Userspace event found {key}")
        return key[:-len(user_space_postfix)]
    return key


def parse_stat_csv(file, force_kernel_space = False):
    import csv

    comments_lines = 2  # validate
    row_val_idx, row_id_idx = 0, 2
    with open(file, newline="") as csvfile:
        [next(csvfile) for i in range(comments_lines)]  # skip headers
        return {
            parse_key(row[row_id_idx], force_kernel_space): row[row_val_idx]
            for row in csv.reader(csvfile, delimiter=",")
        }


# https://perf.wiki.kernel.org/index.php/Tutorial
# http://www.brendangregg.com/perf.html
# or run "perf list"
def events_str(events):
    if len(events) == 0:
        return ""
    return f"-e {events if isinstance(events, str) else ','.join(events)}"


def out_str(output_file):
    return "" if output_file is None else f"-o {os.path.relpath(output_file)}"


def stat_cmd(exe, events, output_file, options=""):
    return f"perf stat -x , {options} {out_str(output_file)} {events_str(events)} {exe}"


def stat(exe, events, output_file=None):
    return run(stat_cmd(exe, events, output_file), check=True)


def record_cmd(exe, events, output_file, options=""):
    return f"perf record {options} {out_str(output_file)} {events_str(events)} {exe}"


def record(exe, events, output_file=None):
    return run(record_cmd(exe, events, output_file), check=True)


def stat_mp(exe, events, output_files):
    return run_mp([stat_cmd(exe, events, out) for out in output_files])

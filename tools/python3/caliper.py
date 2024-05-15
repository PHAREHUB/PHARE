#
#
#

from pathlib import Path
from tools.python3 import extend_env_path, run

SCRIPT_DIR = Path(__file__).resolve().parent

# /home/deegan/git/phare/master/build/subprojects/caliper/src/tools/cali-query
PROJECT_ROOT_DIR = SCRIPT_DIR.parent.parent
DEFAULT_CALI_QUERY_PATH = (
    f"{PROJECT_ROOT_DIR}/build/subprojects/caliper/src/tools/cali-query"
)


# https://hatchet.readthedocs.io/en/latest/basic_tutorial.html#introduction
def hatchet_on_caliper(cali_file="phare_outputs/recorder.0.cali"):
    import hatchet as ht
    import subprocess

    cali_file = "region_profile.cali"

    # Setup desired cali query.
    # cali_query = "cali-query"
    # grouping_attribute = "function"
    # default_metric = "sum(sum#time.duration),inclusive_sum(sum#time.duration)"
    # query = "select function,%s group by %s format json-split" % (
    #     default_metric,
    #     grouping_attribute,
    # )

    # # Use ``cali-query`` here to produce the json-split stream.
    # with extend_env_path(DEFAULT_CALI_QUERY_PATH):
    #     cali_json = subprocess.Popen(
    #         [cali_query, "-q", query, cali_file], stdout=subprocess.PIPE
    #     )

    # # Use hatchet's ``from_caliper`` API with the resulting json-split.
    # # The result is stored into Hatchet's GraphFrame.
    # gf = ht.GraphFrame.from_caliper(cali_json.stdout)

    # # Printout the DataFrame component of the GraphFrame.
    # print(gf.dataframe)

    # # Printout the graph component of the GraphFrame.
    # # Use "time (inc)" as the metric column to be displayed
    # print(gf.tree(metric_column="time (inc)"))

    # # Setup desired cali query.
    # grouping_attribute = "function"
    # default_metric = "sum(sum#time.duration),inclusive_sum(sum#time.duration)"
    # query = "select function,%s group by %s format json-split" % (
    #     default_metric,
    #     grouping_attribute,
    # )

    # # Use hatchet's ``from_caliper`` API with the path to the cali file and the
    # # query. This API will internally run ``cali-query`` on this file to
    # # produce a json-split stream. The result is stored into Hatchet's
    # # GraphFrame.
    # with extend_env_path(DEFAULT_CALI_QUERY_PATH):
    #     gf = ht.GraphFrame.from_caliper(cali_file, query)

    # # Printout the DataFrame component of the GraphFrame.
    # print(gf.dataframe)

    # # Printout the graph component of the GraphFrame.
    # # Use "time (inc)" as the metric column to be displayed
    # print(gf.tree(metric_column="time (inc)"))

    grouping_attribute = "function"
    default_metric = "sum(sum#time.duration),inclusive_sum(sum#time.duration)"
    query = "select function,%s group by %s format json-split" % (
        default_metric,
        grouping_attribute,
    )

    query = "SELECT function,time.duration ORDER BY time.duration FORMAT json-split"
    with extend_env_path(DEFAULT_CALI_QUERY_PATH):
        gf = ht.GraphFrame.from_caliper(cali_file, query)
    print(gf.dataframe)


def hatchet_region_profile_inclusive(time_per_fn_per_rank, ranks):
    fn_ts = {r: {} for r in range(ranks)}

    for rank in range(ranks):
        for k, v in time_per_fn_per_rank[rank].items():
            # print(k, v)
            if "/" not in k:  # is root scope
                fn_ts[rank][k] = v
            else:
                bits = k.split("/")
                fn_ts[rank][bits[0]] += v

    return fn_ts


def hatchet_region_profile(json_file="trace.json", ranks=2, inclusive=False):
    import json

    time_per_rank = {r: 0 for r in range(ranks)}
    time_per_fn_per_rank = {r: {} for r in range(ranks)}
    time_per_fn_per_ts_per_rank = {r: [] for r in range(ranks)}

    ts_idx = 0
    post_init = False
    with open(json_file, "r") as f:
        for d in json.load(f):
            if "event.end#region" not in d:
                continue
            rank = d["mpi.rank"]
            time = d["time"]
            time_per_rank[rank] += time
            if "path" not in d:
                continue
            fn = d["path"]

            if post_init and fn == "Simulator::advance":
                print("post_init")
                ts_idx += 1
            elif not post_init and fn == "Simulator::advance":
                post_init = True

            if fn not in time_per_fn_per_rank[rank]:
                time_per_fn_per_rank[rank][fn] = 0
            time_per_fn_per_rank[rank][fn] += time

            if post_init:
                if ts_idx == len(time_per_fn_per_ts_per_rank[rank]):
                    time_per_fn_per_ts_per_rank[rank].append({})
                time_per_fn_per_ts_per_rank[rank][ts_idx][fn] = time

    if inclusive:
        time_per_fn_per_rank = hatchet_region_profile_inclusive(
            time_per_fn_per_rank, ranks
        )

    print(len(time_per_fn_per_ts_per_rank[0]))
    # for i, d in enumerate(time_per_fn_per_ts_per_rank[0]):
    #     for k, v in d.items():
    #         print(i, k, v)

    print("time_per_rank", time_per_rank)


if __name__ == "__main__":
    hatchet_on_caliper()

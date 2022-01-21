
import rapidjson

""" # sample Run data
    {
      "name": "SolveFixture<2, 1, 1>/_2_1_1_push/20/400/50",
      "family_index": 0,
      "per_family_instance_index": 791,
      "run_name": "SolveFixture<2, 1, 1>/_2_1_1_push/20/400/50",
      "run_type": "iteration",
      "repetitions": 1,
      "repetition_index": 0,
      "threads": 1, # ignored
      "iterations": 2,
      "real_time": 481763362.8845215,
      "cpu_time": 439014839.00003767,
      "time_unit": "ns"
    }
"""


def process_key(json, keys, i, manips, x):
    if i not in manips:
        json[keys[i]] = x
    else:
        v = manips[i]
        if v is not None:
            json[keys[i]] = manips[i](x)



def ETL(json, keys, key_manips = {}, split_adder = lambda x: {}):
    json2 = []
    for run in json["benchmarks"]:
        json2 += [run]

        bits = run["name"].split("/")
        len_bits = len(bits)
        name_split_parts = len_bits - len(keys)

        for ki, i in enumerate(range(name_split_parts, len_bits)):
            process_key(json2[-1], keys, ki, key_manips, bits[i])

        for k, v in split_adder(bits).items():
            json2[-1][k] = v

    return json2

class split_adder:
    def __init__(self, before_keys, val_keys):
        self.before_keys = before_keys
        self.val_keys = val_keys
        self.vals = {}

    def __call__(self, bits):
        for i, k in enumerate(bits[1][1:].split("_")):
            self.vals[self.val_keys[i]] = k
        return self.vals

care_keys = ["ig0", "ig1", "n_parts","cells", "threads"]
del_keys = ["ig0","ig1","family_index","per_family_instance_index","run_name","run_type"]
adder = split_adder(care_keys, ["dim", "interp", "version"])

key_manips= {
  0 : lambda x : None,
  1 : lambda x : None,
  2 : lambda x : int(x),
  3 : lambda x : int(x),
  4 : lambda x : int(x.split(":")[1])}

def main():
    import pandas as pd

    with open("gbench.json") as file:
        json = ETL(rapidjson.loads(file.read()), care_keys, key_manips, adder)
        for j in json:
          for k in del_keys: del j[k]
        with open("gbench2.json", 'w') as out:
            rapidjson.dump(json, out)

    with open("gbench2.json") as file:
        df = pd.DataFrame(rapidjson.loads(file.read()))

        import seaborn as sns
        sns.set(rc={'figure.figsize':(11, 4)})

        df["real_time"].plot(linewidth=.5)

        from pandasgui import show
        show(df)

        import matplotlib.pyplot as plt
        plt.show()

if __name__ == "__main__":
    main()



def plot_scatter(runs):
    # https://stackoverflow.com/a/44985520
    import numpy as np
    import matplotlib.pyplot as plt

    nresults = runs.n
    labels = list(range(nresults))

    str_len  = len(str(nresults))
    y_labels = []
    run_list = []

    for thread_i, ppc_list in runs.data.items():
        for ppc_i, cells_list,  in ppc_list.items():
            for cells_i, run_t,  in cells_list.items():
                y_labels += [run_t[1]["name"]]
                run_list += [run_t[1]]

    axis_legs = [f"{str(i).zfill(str_len)}_{k}" for i, k in enumerate(y_labels)]

    width = 0.4
    fig, ax = plt.subplots(figsize=(35, 15))
    for i, k in enumerate(run_list):
        data = np.array(k["real_time"])
        x = i + (width - width / 2.0)
        ax.scatter(x, data, s=25)
        mean = data.mean()
        ax.plot([i - width / 2.0, i + width / 2.0], [mean, mean], color="k")

    ax.set_xticks(labels)
    ax.set_xticklabels(labels)
    lgd = ax.legend(
        [axis_leg for axis_leg in axis_legs], loc="center left", bbox_to_anchor=(1, 0.5)
    )
    filename = f"gbench_scatter.png"
    out = f"{filename}"
    fig.savefig(out, bbox_extra_artists=(lgd,), bbox_inches="tight")

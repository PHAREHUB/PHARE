
# https://hatchet.readthedocs.io/en/latest/basic_tutorial.html#introduction
def hatchet_on_caliper(filepath):
    import hatchet as ht
    grouping_attribute = "function"
    default_metric = "sum(sum#time.duration),inclusive_sum(sum#time.duration)"
    query = "select function,%s group by %s format json-split" % (
        default_metric,
        grouping_attribute,
    )

    query = "SELECT function,time.duration ORDER BY time.duration FORMAT json-split"
    gf = ht.GraphFrame.from_caliper(filepath, query)
    print(gf.dataframe)

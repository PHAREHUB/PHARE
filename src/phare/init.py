#! /usr/bin/env python3


def get_user_inputs(jobname):
    from phare import populateDict
    import importlib
    print(jobname)
    jobmodule = importlib.import_module(jobname)
    populateDict() # in __init__.py

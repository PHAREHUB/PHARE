def get_user_inputs(jobname):
    import importlib
    from . import populateDict
    import pyphare.pharein as _init_

    _init_.PHARE_EXE = True
    print(jobname)
    jobmodule = importlib.import_module(jobname)  # lgtm [py/unused-local-variable]
    if jobmodule is None:
        raise RuntimeError("failed to import job")
    populateDict()

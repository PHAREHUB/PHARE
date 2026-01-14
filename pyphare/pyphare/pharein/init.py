def get_user_inputs(jobname):
    import sys
    import importlib

    import pyphare.pharein as _init_

    from . import populateDict

    try:
        _init_.PHARE_EXE = True
        jobmodule = importlib.import_module(jobname)  # lgtm [py/unused-local-variable]
        if jobmodule is None:
            raise RuntimeError("failed to import job")
        populateDict()

    except Exception as e:
        import traceback

        print(f"Exception caught in pharein/init::get_user_inputs: \n{e}")
        print(traceback.format_exc())
        sys.exit(1)
    except ...:
        print(f"UNKNOWN Exception caught in pharein/init::get_user_inputs")
        sys.exit(1)

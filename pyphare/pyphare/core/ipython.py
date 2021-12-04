



import weakref
h5_weak_refs = {}

def clean_dead_h5_weak_refs():
    for h5_path, h5_weak_ref in h5_weak_refs.copy().items():
        if h5_weak_ref["ref"]() is None:
            del h5_weak_refs[h5_path]

# https://stackoverflow.com/a/40222538/795574
try:
    def exit_register(fun, *args, **kwargs):
        def callback():
            fun()
            ip.events.unregister('post_execute', callback)
        ip.events.register('post_execute', callback)

    ip = get_ipython() # should fail if not running under ipython

    def add_h5_weak_ref(h5Filepath, h5File):
        import pathlib, datetime, time
        clean_dead_h5_weak_refs()

        h5_weak_refs[h5Filepath] = {
          "ref": weakref.ref(h5File),
          "access": int(time.time()),
          "path": pathlib.Path(h5Filepath)
        }

except NameError:

    def add_h5_weak_ref(h5Filepath, h5File):
        pass # empty function for normal (not jupyter) case

    # forwarding import to keep functionality even if ipython is not in use
    #  it's possible this isn't even needed, and we just want to execute the callback for only if
    #    we're running under ipython, in that case, "exit_register" -> "ipython_exit_register" and
    #     replace the next line with def ipython_exit_register(fun, *args, **kwargs): pass
    from atexit import register as exit_register # lgtm [py/unused-import]




@exit_register
def close_h5_files_refs():
    import pathlib, datetime, time
    for h5_weak_ref_dict in list(h5_weak_refs.values()):
        h5_ref = h5_weak_ref_dict["ref"]()
        if h5_ref is not None:
            mtime = int(h5_weak_ref_dict["path"].stat().st_mtime)
            if mtime > h5_weak_ref_dict["access"] and h5_ref.__bool__():
                h5_ref.close()

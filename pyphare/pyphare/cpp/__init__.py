import importlib
import json
import os

_libs = {}


def _simulator_id_parts(sim):
    """The build permutation fields, in the order the res/sim files list them."""
    parts = [str(sim.ndim)]

    if sim.interp_order:
        parts += [str(sim.interp_order), str(sim.refined_particle_nbr)]
        parts += [sim.model_layout, sim.allocator]
    elif sim.mhd_timestepper:  # MHD only
        parts += [sim.model_layout, sim.allocator]

    if sim.mhd_timestepper:
        hall_active = "true" if sim.hall else "false"
        parts += [
            sim.reconstruction,
            sim.limiter,
            sim.riemann,
            hall_active,
        ]

    return parts


def simulator_id(sim):
    return "_".join(_simulator_id_parts(sim))


def permutation_line(sim):
    """The res/sim permutation file line that builds this simulation's module.
    Not derivable from simulator_id by splitting: field values contain '_'
    (SSPRK4_5)."""
    return ",".join(_simulator_id_parts(sim))


def built_simulator_ids():
    """Permutations present in the loaded build, or None if pybindlibs itself is
    not importable (no build on PYTHONPATH at all)."""
    import pkgutil

    try:
        import pybindlibs
    except ModuleNotFoundError:
        return None

    return sorted(
        name[len("cpp_") :]
        for _, name, _ in pkgutil.iter_modules(pybindlibs.__path__)
        if name.startswith("cpp_") and name != "cpp_etc"
    )


def _no_such_module_message(sim, mod_str):
    built = built_simulator_ids()
    if built is None:
        listing = "    pybindlibs not importable -- no build on PYTHONPATH"
    elif built:
        listing = "\n".join(f"    {name}" for name in built)
    else:
        listing = "    (none)"
    return (
        f"No module named '{mod_str}'.\n"
        f"This run needs the build permutation\n"
        f"    {permutation_line(sim)}\n"
        f"Add it to a res/sim permutation file (res/sim/all.txt) and rebuild.\n"
        f"Permutations currently built:\n{listing}"
    )


def cpp_lib(sim):
    global _libs

    mod_id = f"cpp_{simulator_id(sim)}"
    mod_str = f"pybindlibs.{mod_id}"
    if mod_str not in _libs:
        try:
            _libs[mod_str] = importlib.import_module(mod_str)
        except ModuleNotFoundError as e:
            if not mpi_initialized() or mpi_size() == 1:
                from .cmake import build

                build(mod_id)
                return cpp_lib(sim)
            raise ModuleNotFoundError(_no_such_module_message(sim, mod_str)) from e
    return _libs[mod_str]


def cpp_etc_lib():
    lib = importlib.import_module("pybindlibs.cpp_etc")
    if "cpp_etc" not in _libs:
        _libs["cpp_etc"] = lib
        _check_model_layouts(lib)
    return lib


def _check_model_layouts(lib):
    """supported_model_layouts() is a python copy, it must match the C++ list"""
    cpp_layouts = [str(l).split(".")[-1] for l in lib.supported_layouts()]
    assert supported_model_layouts() == cpp_layouts, (
        f"layout mismatch: {supported_model_layouts()} != {cpp_layouts}"
    )


def build_config():
    return cpp_etc_lib().phare_build_config()


def build_config_as_json():
    return json.dumps(build_config())


def splitter_type(sim):
    return cpp_lib(sim).Splitter


def split_pyarrays_fn(sim):
    return cpp_lib(sim).split_pyarray_particles


def mpi_rank():
    return cpp_etc_lib().mpi_rank()


def mpi_size():
    return cpp_etc_lib().mpi_size()


def mpi_barrier():
    return cpp_etc_lib().mpi_barrier()


def mpi_initialized():
    return cpp_etc_lib().mpi_initialized()


def mpi_hostname():
    return cpp_etc_lib().mpi_hostname()


def mpi_is_root_rank_for_host():
    """collective on first call, every rank must make the first call"""
    return cpp_etc_lib().mpi_is_root_rank_for_host()


def is_rank0():
    """works before MPI is initialized, or without the C++ module, falling back to rank 0"""
    try:
        if mpi_initialized():
            return mpi_rank() == 0
    except ImportError:
        # missing module or mpi not initialized
        ...
    # Open MPI, PMIx, MPICH/Intel MPI (hydra), slurm srun
    envs = ["OMPI_COMM_WORLD_RANK", "PMIX_RANK", "PMI_RANK", "SLURM_PROCID"]
    for env in envs:
        if env in os.environ:
            return int(os.environ[env]) == 0
    return True  # FALL BACK ALWAYS RANK 0


def print_rank0(*args, **kwargs):
    if is_rank0():
        print(*args, **kwargs)


def layout_mode_enum():
    return cpp_etc_lib().LayoutMode()


def supported_model_layouts():
    """python copy of the C++ list, usable without a build (e.g. at test import)
    checked against C++ on first load of cpp_etc"""
    from pyphare.pharein.simulation import valid_model_layouts

    return list(valid_model_layouts)


def supported_mhd_model_layouts():
    """model layouts of MHD only simulations, untiled (AoS) or tiled (AoSTS)"""
    from pyphare.pharein.simulation import valid_mhd_model_layouts

    return list(valid_mhd_model_layouts)


def simulation_to_simopts(sim):
    lib = cpp_etc_lib()
    opts = lib.SimOpts()
    opts.dimension = sim.ndim
    opts.interp_order = sim.interp_order
    opts.nbRefinedPart = sim.refined_particle_nbr
    opts.layout_mode = getattr(lib.LayoutMode, sim.model_layout)
    opts.alloc_mode = getattr(lib.AllocatorMode, sim.allocator)

    if getattr(sim, "mhd_timestepper", None):
        opts.time_integrator_type = getattr(lib.TimeIntegratorType, sim.mhd_timestepper)
        opts.reconstruction_type = getattr(lib.ReconstructionType, sim.reconstruction)
        opts.slope_limiter_type = getattr(lib.SlopeLimiterType, sim.limiter)
        opts.riemann_solver_type = getattr(lib.RiemannSolverType, sim.riemann)
        opts.Hall = sim.hall
        opts.Resistivity = sim.res
        opts.HyperResistivity = sim.hyper_res

    return opts

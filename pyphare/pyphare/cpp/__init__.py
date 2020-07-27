

def splitter_type(dim, interp, n_particles):
    from pybindlibs import cpp
    return getattr(cpp,
        "Splitter_" + str(dim) + "_" + str(interp)+ "_" + str(n_particles),
    )

def create_splitter(dim, interp, n_particles):
    return splitter_type(dim, interp, n_particles)()


def split_pyarrays_fn(dim, interp, n_particles):
    from pybindlibs import cpp
    return getattr(cpp,
        "split_pyarray_particles_" + str(dim) + "_" + str(interp)+ "_" + str(n_particles),
    )

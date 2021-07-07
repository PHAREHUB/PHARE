

def timestamps_with_step(sim, dump_step):
    import numpy as np
    nbr_dump_step = sim.final_time / dump_step
    return dump_step * np.arange(nbr_dump_step)

def all_timestamps(sim):
    import numpy as np
    nbr_dump_step = int(sim.final_time / sim.time_step) + 1
    return sim.time_step * np.arange(nbr_dump_step)


def dump_all_diags(pops=[], flush_every=100):
    import pyphare.pharein as ph, numpy as np

    sim = ph.global_vars.sim

    timestamps = all_timestamps(sim)

    for quantity in ["density", "bulkVelocity"]:
        ph.FluidDiagnostics(
            quantity=quantity,
            write_timestamps=timestamps,
            compute_timestamps=timestamps,
            flush_every=flush_every,
        )

    for pop in pops:
        for quantity in ["density", "flux"]:
          ph.FluidDiagnostics(
              quantity=quantity,
              write_timestamps=timestamps,
              compute_timestamps=timestamps,
              flush_every=flush_every,
              population_name=pop
          )

        for quantity in ['domain', 'levelGhost', 'patchGhost']:
            ph.ParticleDiagnostics(
                quantity=quantity,
                compute_timestamps=timestamps,
                write_timestamps=timestamps,
                flush_every=flush_every,
                population_name=pop
            )

    for quantity in ["E", "B"]:
        ph.ElectromagDiagnostics(
            quantity=quantity,
            write_timestamps=timestamps,
            compute_timestamps=timestamps,
            flush_every=flush_every,
        )



def dump_all_diags(pops=[]):
    import phare.pharein as ph, numpy as np

    sim = ph.globals.sim

    timestamps = np.arange(0, sim.time_step+3)*sim.time_step


    for type in ["density", "bulkVelocity"]:
        ph.FluidDiagnostics(
            diag_type=type,
            write_timestamps=timestamps,
            compute_timestamps=timestamps,
        )

    for pop in pops:
        for type in ["density", "flux"]:
          ph.FluidDiagnostics(
              diag_type=type,
              write_timestamps=timestamps,
              compute_timestamps=timestamps,
              population_name=pop
          )

        for type in ['domain', 'levelGhost', 'patchGhost']:
            ph.ParticleDiagnostics(
                diag_type=type,
                compute_timestamps=timestamps,
                write_timestamps=timestamps,
                population_name=pop
            )

    for type in ["E", "B"]:
        ph.ElectromagDiagnostics(
            diag_type=type,
            write_timestamps=timestamps,
            compute_timestamps=timestamps,
        )

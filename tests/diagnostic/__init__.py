

def dump_all_diags(pops=[]):
    import phare.pharein as ph, numpy as np

    sim = ph.globals.sim

    timestamps = np.arange(0, sim.time_step+3)*sim.time_step

    for type in ["density", "bulkVelocity"]:
        ph.FluidDiagnostics(
            diag_type=type,
            write_timestamps=timestamps,
            compute_timestamps=timestamps,
            start_iteration=0,
            last_iteration=990,
        )

    for pop in pops:
        for type in ["density", "flux"]:
          ph.FluidDiagnostics(
              diag_type=type,
              write_timestamps=timestamps,
              compute_timestamps=timestamps,
              start_iteration=0,
              last_iteration=990,
              population_name=pop
          )

        for type in ['domain', 'levelGhost', 'patchGhost']:
            ph.ParticleDiagnostics(
                diag_type=type,
                compute_timestamps=timestamps,
                write_timestamps=timestamps,
                start_iteration=0,
                last_iteration=90,
                population_name=pop
            )

    for type in ["E", "B"]:
        ph.ElectromagDiagnostics(
            diag_type=type,
            write_timestamps=timestamps,
            compute_timestamps=timestamps,
            start_teration=0,
            last_iteration=990
        )

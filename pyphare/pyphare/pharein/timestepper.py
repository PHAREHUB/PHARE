"""
Time-step resolution and validation for pharein.Simulation.

'time_step' (the public Simulation constructor option) may be:
  - a scalar (or absent)                 -> ConstantTimeStepper
  - {'mode': 'constant', 'value': <dt>}  -> ConstantTimeStepper ('value' optional)
  - {'mode': 'adaptive', 'cfl_wave': <w>, 'cfl_diffusive': <d>}
                                         -> AdaptiveTimeStepper
                                            ('cfl_diffusive' defaults to 'cfl_wave')

resolve_time_stepper(**kwargs) is the single entry point used by Simulation's checker();
it returns the TimeStepper instance stored as Simulation.time_stepper.
"""

from abc import ABC
from dataclasses import dataclass

from ..core import phare_utilities


def _restart_start_time(restart_options):
    return (restart_options or {}).get("restart_time", 0)


@dataclass
class TimeStepper(ABC):
    """
    Base class holding the timing parameters common to every time-stepping mode.
    """

    mode = None

    start_time: float
    final_time: float

    def within_simulation_duration(self, time_period):
        raise NotImplementedError

    def populate_dict(self, dp):
        """Mirror the public `time_step` dict shape (mode + per-mode params) on the C++ side.

        `dp` is a DictPopulator (see pharein.initialize): an object exposing
        add_string/add_double/add_int/add_enum_int/... - passed in rather than imported, to
        avoid a circular import between this module and pharein.initialize.
        """
        dp.add_enum_int("simulation/time_step/mode", "TimeStepType", self.mode)
        dp.add_double("simulation/final_time", self.final_time)


@dataclass
class ConstantTimeStepper(TimeStepper):
    mode = "constant"
    time_step: float
    time_step_nbr: int

    def within_simulation_duration(self, time_period):
        return time_period[0] >= 0 and time_period[1] < self.time_step_nbr

    def populate_dict(self, dp):
        super().populate_dict(dp)
        dp.add_double("simulation/time_step/value", self.time_step)
        dp.add_int("simulation/time_step_nbr", self.time_step_nbr)


@dataclass
class AdaptiveTimeStepper(TimeStepper):
    mode = "adaptive"

    cfl_wave: float = None
    cfl_diffusive: float = None

    def within_simulation_duration(self, time_period):
        raise NotImplementedError(
            "adaptive dt has no fixed time_step_nbr - duration is checked at runtime by the C++ side"
        )

    def populate_dict(self, dp):
        super().populate_dict(dp)
        # dt is computed each step on the C++ side; bound the run by final_time
        dp.add_double("simulation/time_step/cfl_wave", self.cfl_wave)
        dp.add_double("simulation/time_step/cfl_diffusive", self.cfl_diffusive)


# ------------------------------------------------------------------------------


def _resolve_constant_time_step(start_time, time_step, time_step_nbr, final_time):
    final_and_dt = (
        final_time is not None and time_step is not None and time_step_nbr is None
    )
    nsteps_and_dt = (
        time_step_nbr is not None and time_step is not None and final_time is None
    )
    final_and_nsteps = (
        final_time is not None and time_step_nbr is not None and time_step is None
    )

    if sum([final_and_dt, final_and_nsteps, nsteps_and_dt]) != 1:
        raise ValueError(
            "Error: Specify either 'final_time' and 'time_step' or 'time_step_nbr' and 'time_step'"
            + " or 'final_time' and 'time_step_nbr'"
        )

    if final_time is None:
        final_time = start_time + time_step * time_step_nbr

    total_time = final_time - start_time
    if total_time < 0:
        raise RuntimeError("Simulation time cannot be negative - review inputs")

    if phare_utilities.fp_equal(total_time, 0):
        return ConstantTimeStepper(start_time, final_time, 1.0, 0)

    if final_and_dt:
        # keep dt bit-exact, round to the nearest whole step count, and derive final_time
        # from it (identical to C++ finalTime_ = start + dt*nbr). Flooring instead would
        # leave the user's raw final_time slightly above nbr*dt, so a dump grid built at
        # multiples of dt up to final_time overshoots it. See master check_time().
        time_step_nbr = int(round(total_time / time_step))
        final_time = start_time + time_step_nbr * time_step
    elif final_and_nsteps:
        time_step = total_time / time_step_nbr
    # else nsteps_and_dt: time_step and time_step_nbr are already both given

    return ConstantTimeStepper(start_time, final_time, time_step, time_step_nbr)


def _resolve_dict_time_step(ts, *, start_time, final_time, time_step_nbr):
    valid_modes = ("constant", "adaptive")
    mode = ts.get("mode")
    if mode not in valid_modes:
        raise ValueError(
            f"Error: time_step dict requires 'mode' in {valid_modes}, got {mode!r}"
        )

    def _check_keys(allowed):
        extra = set(ts) - allowed
        if extra:
            raise ValueError(
                f"Error: invalid time_step keys for mode '{mode}': {sorted(extra)}, "
                f"allowed {sorted(allowed)}"
            )

    if mode == "constant":
        _check_keys({"mode", "value"})
        return _resolve_constant_time_step(
            start_time=start_time,
            time_step=ts.get("value"),
            time_step_nbr=time_step_nbr,
            final_time=final_time,
        )

    # adaptive: dt is recomputed each step from a CFL constraint. 'time_step_nbr' is not
    # allowed (imposing a step count with a variable dt is out of scope for now).
    _check_keys({"mode", "cfl_wave", "cfl_diffusive"})
    if "cfl_wave" not in ts:
        raise ValueError("Error: adaptive time_step requires 'cfl_wave'")
    if time_step_nbr is not None:
        raise ValueError(
            "Error: adaptive time_step is incompatible with a constant 'time_step' / 'time_step_nbr'"
        )
    if final_time is None:
        raise ValueError("Error: adaptive time_step requires 'final_time'")

    cfl_wave = ts["cfl_wave"]
    cfl_diffusive = ts.get("cfl_diffusive", cfl_wave)
    if cfl_wave <= 0:
        raise ValueError("Error: adaptive time_step 'cfl_wave' must be > 0")
    if cfl_diffusive <= 0:
        raise ValueError("Error: adaptive time_step 'cfl_diffusive' must be > 0")
    if final_time - start_time < 0:
        raise RuntimeError("Simulation time cannot be negative - review inputs")

    return AdaptiveTimeStepper(start_time, final_time, cfl_wave, cfl_diffusive)


def resolve_time_stepper(**kwargs):
    """
    Resolve the public 'time_step' / 'time_step_nbr' / 'final_time' Simulation options
    (plus 'restart_options' for the start time) into a validated TimeStepper.

    'time_step' may also directly be an already-resolved TimeStepper instance, in which case
    it is returned as-is (used by callers that build one themselves, e.g. tools/bench).
    """
    ts = kwargs.get("time_step")
    if isinstance(ts, TimeStepper):
        return ts

    start_time = _restart_start_time(kwargs.get("restart_options"))
    final_time = kwargs.get("final_time")
    time_step_nbr = kwargs.get("time_step_nbr")

    if isinstance(ts, dict):
        return _resolve_dict_time_step(
            ts, start_time=start_time, final_time=final_time, time_step_nbr=time_step_nbr
        )

    return _resolve_constant_time_step(
        start_time=start_time,
        time_step=ts,
        time_step_nbr=time_step_nbr,
        final_time=final_time,
    )

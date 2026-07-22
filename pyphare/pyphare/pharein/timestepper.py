"""
Time-step resolution and validation for pharein.Simulation.

'time_step' (the public Simulation constructor option) may be:
  - a scalar (or absent)                              -> ConstantTimeStepper
  - {'mode': 'constant', 'value': <dt>}                -> ConstantTimeStepper ('value' optional)
  - {'mode': 'adaptive', 'cfl': <c>, 'fourier': <f>}    -> AdaptiveTimeStepper ('fourier' defaults to 'cfl')

resolve_time_stepper(**kwargs) is the single entry point used by Simulation's checker()
pipeline; it returns the TimeStepper instance stored as Simulation.time_stepper.
"""

from abc import ABC
from dataclasses import dataclass
from typing import Optional

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
    time_step: Optional[float]
    time_step_nbr: Optional[int]

    # hard coded in C++ MultiPhysicsIntegrator::getMaxFinerLevelDt
    nSubcycles = 4

    def __post_init__(self):
        self.level_time_steps = None
        self.level_step_nbr = None

    def resolve_levels(self, max_nbr_levels):
        """Populate per-level dt / coarse-step-count arrays used for AMR subcycling."""
        raise NotImplementedError

    def within_simulation_duration(self, time_period):
        raise NotImplementedError

    def populate_dict(self, dp):
        """Mirror the public `time_step` dict shape (mode + per-mode params) on the C++ side.

        `dp` is a dict_populator() (see pharein.initialize.general): an object exposing
        add_string/add_double/add_int/... - passed in rather than imported, to avoid a
        circular import between this module and pharein.initialize.general.
        """
        dp.add_string("simulation/time_step/mode", self.mode)


@dataclass
class ConstantTimeStepper(TimeStepper):
    mode = "constant"

    def resolve_levels(self, max_nbr_levels):
        step_diff = 1 / self.nSubcycles
        self.level_time_steps = [
            self.time_step * (step_diff**ilvl) for ilvl in range(max_nbr_levels)
        ]
        self.level_step_nbr = [
            self.nSubcycles**ilvl * self.time_step_nbr for ilvl in range(max_nbr_levels)
        ]

    def within_simulation_duration(self, time_period):
        return time_period[0] >= 0 and time_period[1] < self.time_step_nbr

    def populate_dict(self, dp):
        super().populate_dict(dp)
        dp.add_double("simulation/time_step/value", self.time_step)
        dp.add_int("simulation/time_step_nbr", self.time_step_nbr)


@dataclass
class AdaptiveTimeStepper(TimeStepper):
    mode = "adaptive"

    cfl: float = None
    fourier: float = None

    def resolve_levels(self, max_nbr_levels):
        pass  # dt (and per-level step counts) are unknown ahead of the run

    def within_simulation_duration(self, time_period):
        raise NotImplementedError(
            "adaptive dt has no fixed time_step_nbr - duration is checked at runtime by the C++ side"
        )

    def populate_dict(self, dp):
        super().populate_dict(dp)
        # dt is computed each step on the C++ side; bound the run by final_time
        dp.add_double("simulation/time_step/cfl", self.cfl)
        dp.add_double("simulation/time_step/fourier", self.fourier)
        dp.add_double("simulation/final_time", self.final_time)


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
        return ConstantTimeStepper(start_time, final_time, 0, 0)

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
    _check_keys({"mode", "cfl", "fourier"})
    if "cfl" not in ts:
        raise ValueError("Error: adaptive time_step requires 'cfl'")
    if time_step_nbr is not None:
        raise ValueError(
            "Error: adaptive time_step is incompatible with a constant 'time_step' / 'time_step_nbr'"
        )
    if final_time is None:
        raise ValueError("Error: adaptive time_step requires 'final_time'")

    cfl = ts["cfl"]
    fourier = ts.get("fourier", cfl)
    if cfl <= 0:
        raise ValueError("Error: adaptive time_step 'cfl' must be > 0")
    if fourier <= 0:
        raise ValueError("Error: adaptive time_step 'fourier' must be > 0")
    if final_time - start_time < 0:
        raise RuntimeError("Simulation time cannot be negative - review inputs")

    return AdaptiveTimeStepper(start_time, final_time, None, None, cfl, fourier)


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

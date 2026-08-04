#
#

from dataclasses import dataclass, field


@dataclass
class LoadBalancer:
    # Levels other than L0 are never explicitly rebalanced by this class - they get
    # rebalanced as a side effect of regridding, whenever that naturally happens.
    # Only L0 is rebalanced "on demand" by the logic these options configure.

    # which way load is assessed - applies to every level, not just L0
    mode: str = field(default_factory=lambda: "nppc")

    # acceptable imbalance essentially - used both for L0's explicit rebalance
    # decision and as the threshold applied to every other level when it is
    # regridded. currently a single value shared by all levels; there is no
    # per-level override
    tol: float = field(default_factory=lambda: 0.05)

    # L0 only: whether to force a rebalance check of L0 at the very first step
    on_init: bool = field(default_factory=lambda: True)

    # L0 only: everything below this line governs *when* L0's explicit rebalance
    # check runs - none of it affects other levels

    # if auto, other values are not used if active
    auto: bool = field(default_factory=lambda: False)
    next_rebalance_backoff_multiplier: int = field(default_factory=lambda: 2)
    next_rebalance: int = field(default_factory=lambda: 200)
    max_next_rebalance: int = field(default_factory=lambda: 1000)

    # if !auto these values are used if active
    every: int = field(default_factory=lambda: None)

    # internal, allows not registering object for default init
    _register: bool = field(default_factory=lambda: True)

    def __post_init__(self):
        if self.auto and self.every:
            raise RuntimeError("LoadBalancer cannot work with both 'every' and 'auto'")

        if self.every is None:
            self.auto = True
            self.every = 0  # python3 -> c++ doesn't understand 'None'

        allowed_modes = [
            "nppc",  # count particles per rank
            "homogeneous",  # count cells per rank
        ]

        if self.mode not in allowed_modes:
            raise RuntimeError(f"LoadBalancer mode '{self.mode}' is not valid")

        if self._register:
            from . import global_vars as gv

            if not gv.sim:
                raise RuntimeError(
                    "LoadBalancer cannot be registered as no simulation exists"
                )
            if gv.sim.load_balancer:
                raise RuntimeError("LoadBalancer is already registered to simulation")
            gv.sim.load_balancer = self

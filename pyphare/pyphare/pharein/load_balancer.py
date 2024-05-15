#
#

from dataclasses import dataclass, field
from . import global_vars as gv


@dataclass
class LoadBalancer:
    # whether or not load balancing is performed
    active: bool = field(default_factory=lambda: False)

    # which way load is assessed
    mode: str = field(default_factory=lambda: "nppc")

    # acceptable imbalance essentially
    tol: float = field(default_factory=lambda: 0.05)

    # whether to rebalance/check imbalance on init
    on_init: bool = field(default_factory=lambda: True)

    # if auto, other values are not used if active
    auto: bool = field(default_factory=lambda: True)
    next_rebalance_backoff_multiplier: int = field(default_factory=lambda: 2)
    next_rebalance: int = field(default_factory=lambda: 200)
    max_next_rebalance: int = field(default_factory=lambda: 1000)

    # if !auto these values are used if active
    every: int = field(default_factory=lambda: 1)

    # internal, allows not registering object for default init
    _register: bool = field(default_factory=lambda: True)

    def __post_init__(self):
        allowed_modes = [
            "nppc",  # count particles per rank
            "homogeneous",  # count cells per rank
        ]

        if self.mode not in allowed_modes:
            raise RuntimeError(f"LoadBalancer mode '{self.mode}' is not valid")

        if self._register:
            if not gv.sim:
                raise RuntimeError(
                    f"LoadBalancer cannot be registered as no simulation exists"
                )
            if gv.sim.load_balancer:
                raise RuntimeError(f"LoadBalancer is already registered to simulation")
            gv.sim.load_balancer = self

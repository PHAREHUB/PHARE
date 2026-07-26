#
#

from dataclasses import dataclass, field
from . import global_vars as gv


@dataclass
class LoadBalancer:
    """
    LoadBalancer controls how simulation work (particles and cells) is
    redistributed across MPI ranks as the simulation runs, to keep ranks
    evenly loaded on an evolving, adaptively-refined domain. Optional - if
    none is created, load balancing runs with the defaults below.

    Can be created any time after :class:`~pyphare.pharein.Simulation`, at
    most once per simulation.

    ``mode`` shapes how load is measured wherever patches get redistributed
    across ranks, at any level. ``tol`` does too, but is also reused as the
    threshold for a separate, proactive rebalancing of the coarsest level
    (L0) between regrids. ``auto``, ``every`` and ``on_init`` control only
    that L0 proactive rebalancing - finer levels are redistributed across
    ranks only as part of their own AMR regrid, on their own cadence, never
    on this one.

    **Usage example:**

    .. code-block:: python

        from pyphare.pharein import LoadBalancer

        LoadBalancer(mode="nppc", tol=0.05, every=200)

    **Parameters**:

        * **active** (``bool``), default=True, whether the proactive L0 rebalancing described below is performed at all.
        * **mode** (``str``), {"nppc" (default), "homogeneous"}, how load is assessed per rank: "nppc" counts particles per rank, "homogeneous" counts cells per rank. Applies at every level, whenever that level's patches are redistributed across ranks.
        * **tol** (``float``), default=0.05, acceptable relative imbalance between ranks. Used both hierarchy-wide (only used with `mode="nppc"`) and as the threshold that triggers a proactive L0 rebalance between regrids.
        * **on_init** (``bool``), default=True, whether L0 is checked/rebalanced right after the simulation is initialized.
        * **auto** (``bool``), default=True unless `every` is given, whether L0's proactive rebalancing cadence is chosen automatically (starting at `next_rebalance` steps, backing off by `next_rebalance_backoff_multiplier` up to `max_next_rebalance`) rather than at a fixed cadence. Mutually exclusive with `every`.
        * **every** (``int``), if set, proactively rebalance L0 every `every` steps instead of automatically. Mutually exclusive with `auto`.
    """

    # whether or not load balancing is performed
    active: bool = field(default_factory=lambda: True)

    # which way load is assessed
    mode: str = field(default_factory=lambda: "nppc")

    # acceptable imbalance essentially
    tol: float = field(default_factory=lambda: 0.05)

    # whether to rebalance/check imbalance on init
    on_init: bool = field(default_factory=lambda: True)

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
            if not gv.sim:
                raise RuntimeError(
                    "LoadBalancer cannot be registered as no simulation exists"
                )
            if gv.sim.load_balancer:
                raise RuntimeError("LoadBalancer is already registered to simulation")
            gv.sim.load_balancer = self

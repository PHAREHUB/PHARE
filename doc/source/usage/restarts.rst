
========
Restarts
========

``restart_options``, passed to :class:`~pyphare.pharein.Simulation` (see
:doc:`simulation`), controls two things: writing checkpoint files during a
run, and resuming a later run from one of them.

Writing checkpoints
---------------------

.. code-block:: python

    import numpy as np

    Simulation(
        # ...
        restart_options={
            "dir": "restarts",
            "mode": "conserve",
            "timestamps": np.array([1.0, 2.0, 3.0]),
        },
    )

- ``mode="conserve"`` (default) keeps any existing checkpoint files;
  ``mode="overwrite"`` replaces them.
- ``timestamps`` are simulation times, and must line up with the
  simulation's ``time_step``.
- ``elapsed_timestamps`` writes a checkpoint every time that much wall-clock
  time has passed instead (in seconds, or as ``datetime.timedelta``) - handy
  on a cluster with a wall-time limit, since it doesn't require knowing in
  advance which simulation time will be reached.
- ``keep_last=N`` automatically deletes older checkpoint directories,
  keeping only the ``N`` most recent - useful to bound disk usage on a long
  run that checkpoints often.

.. tip::

   On supercomputers, jobs are typically killed after a fixed wall time
   (e.g. 24h), but you can rarely predict which simulation time the run
   will have reached by then. Prefer ``elapsed_timestamps`` over
   ``timestamps`` in that case, e.g. at 10h and 23h, so a checkpoint is
   guaranteed to be dumped with enough time left before the job is
   killed - rather than picking simulation timestamps that risk falling
   right after the job dies.

Resuming from a checkpoint
-----------------------------

.. code-block:: python

    Simulation(
        # ...
        restart_options={
            "dir": "restarts",
            "mode": "conserve",
            "restart_time": "auto",  # or an explicit time, e.g. 2.0
        },
    )

``restart_time="auto"`` picks the latest available checkpoint in ``dir``.
All other ``Simulation`` parameters (domain, resolution, ``max_nbr_levels``,
...) must match the run that produced the checkpoint.

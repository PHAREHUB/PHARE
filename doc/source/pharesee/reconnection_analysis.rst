
Reconnection analysis
=====================

Specialized tools for analyzing 2D magnetic reconnection simulations.
These methods are only applicable to 2D runs that include magnetic field
diagnostics (``EM_B``).


Computing the vector potential
-------------------------------

The z-component of the magnetic vector potential :math:`A_z` satisfies
:math:`\mathbf{B} = \nabla \times \mathbf{A}`, so in 2D:

.. math::

   B_x = \frac{\partial A_z}{\partial y}, \qquad
   B_y = -\frac{\partial A_z}{\partial x}

:math:`A_z` is computed by numerically integrating the magnetic field on the
finest available AMR level.

.. code-block:: python

    Az, (xn, yn) = run.GetMagneticFlux(time, interp="nearest")

**Signature**::

    Run.GetMagneticFlux(time, interp="nearest", xn=None, yn=None, Xn=None, Yn=None)

**Returns** ``(Az, (xn, yn))``:

- ``Az`` — :class:`numpy.ndarray` of shape ``(len(xn), len(yn))`` containing
  :math:`A_z` values on the finest-level uniform grid.
- ``xn`` — 1-D :class:`numpy.ndarray` of x-coordinates.
- ``yn`` — 1-D :class:`numpy.ndarray` of y-coordinates.

The optional ``xn``, ``yn``, ``Xn``, ``Yn`` parameters allow reusing
pre-computed coordinate grids across multiple time steps, which avoids
rebuilding the meshgrid on every call (see :ref:`reconnection-worked-example`
and :meth:`~pyphare.pharesee.run.Run.GetReconnectionRate` for how this is
exploited internally).


Finding X-points
-----------------

An X-point is a saddle point of :math:`A_z`.
``FindPrimaryXPoint`` locates the X-point closest to the domain center by:

1. Computing the magnitude of :math:`\nabla A_z`.
2. Retaining candidate cells whose gradient magnitude is below the 5th
   percentile (near-zero gradient).
3. Among those candidates, selecting the cell with the most negative
   determinant of the Hessian of :math:`A_z` (strongest saddle character).

.. code-block:: python

    x_xpoint, y_xpoint, idx = run.FindPrimaryXPoint(Az, xn, yn)

**Signature**::

    Run.FindPrimaryXPoint(Az, xn, yn)

**Parameters**:

- ``Az`` — 2-D :class:`numpy.ndarray`, the vector potential array returned by
  :meth:`~pyphare.pharesee.run.Run.GetMagneticFlux`.
- ``xn``, ``yn`` — 1-D coordinate arrays for the two axes.

**Returns** a 3-tuple ``(x_xpoint, y_xpoint, idx)``:

- ``x_xpoint`` — float, x-coordinate of the primary X-point.
- ``y_xpoint`` — float, y-coordinate of the primary X-point.
- ``idx`` — 2-tuple of ints, array index ``(ix, iy)`` locating the X-point in
  ``Az``.


Reconnection rate
------------------

The reconnection rate is the time derivative of the magnetic flux threading
the X-point, :math:`\mathrm{d}A_z/\mathrm{d}t\big|_{\text{X-point}}`.

.. code-block:: python

    times = [0.0, 10.0, 20.0, 30.0]
    times_centered, rates, flux_at_xpoint, xpoint_trajectory = run.GetReconnectionRate(times)

**Signature**::

    Run.GetReconnectionRate(times, interp="nearest")

**Parameters**:

- ``times`` — list or array of simulation times at which to evaluate the flux.
  Must contain at least two entries.
- ``interp`` — interpolation method forwarded to
  :meth:`~pyphare.pharesee.run.Run.GetMagneticFlux` (default ``"nearest"``).

**Returns** a 4-tuple ``(times_centered, rates, flux_at_xpoint, xpoint_trajectory)``:

- ``times_centered`` — :class:`numpy.ndarray` of length ``len(times) - 1``,
  mid-points of consecutive time intervals.
- ``rates`` — :class:`numpy.ndarray` of length ``len(times) - 1``,
  finite-difference reconnection rate :math:`\Delta A_z / \Delta t` at each
  interval.
- ``flux_at_xpoint`` — :class:`numpy.ndarray` of length ``len(times)``,
  value of :math:`A_z` at the primary X-point at each input time.
- ``xpoint_trajectory`` — :class:`numpy.ndarray` of shape
  ``(len(times), 2)``, X-point position ``[x, y]`` at each input time.

The coordinate grids are built once from the finest level at ``times[0]`` and
reused for all subsequent times, so the method is efficient for long time
series.


.. _reconnection-worked-example:

Worked example
--------------

The following script performs a complete reconnection analysis on a 2D Harris
current sheet run.  It assumes the simulation has been configured and executed
with magnetic field diagnostics enabled (see :doc:`../tutorials/harris_2d`).

.. code-block:: python

    import numpy as np
    import matplotlib.pyplot as plt

    from pyphare.pharesee.run import Run

    # --- 1. Load data ----------------------------------------------------------
    run = Run("/data/scratch/phare_jobs/harris2d")

    # --- 2. Compute Az at several times ----------------------------------------
    snap_times = [0.0, 5.0, 10.0, 15.0, 20.0]

    fig_contour, axes = plt.subplots(1, len(snap_times), figsize=(16, 4), sharey=True)

    for ax, t in zip(axes, snap_times):
        Az, (xn, yn) = run.GetMagneticFlux(t, interp="nearest")

        # --- 3. Find X-point ---------------------------------------------------
        x_xp, y_xp, idx = run.FindPrimaryXPoint(Az, xn, yn)

        # --- 4. Plot Az contours with X-point marked ---------------------------
        Xn, Yn = np.meshgrid(xn, yn, indexing="ij")
        ax.contour(Xn, Yn, Az, levels=30, colors="steelblue", linewidths=0.8)
        ax.plot(x_xp, y_xp, "r+", markersize=10, markeredgewidth=2, label="X-point")
        ax.set_title(f"t = {t:.1f}")
        ax.set_xlabel("x")

    axes[0].set_ylabel("y")
    axes[-1].legend()
    fig_contour.suptitle(r"$A_z$ contours — Harris sheet")
    fig_contour.tight_layout()
    fig_contour.savefig("Az_contours.png", dpi=150)

    # --- 5. Compute and plot reconnection rate vs time -------------------------
    rate_times = np.linspace(0.0, 20.0, 41)   # dense time series
    times_c, rates, flux, traj = run.GetReconnectionRate(rate_times.tolist())

    fig_rate, (ax_flux, ax_rate) = plt.subplots(2, 1, figsize=(8, 6), sharex=True)

    ax_flux.plot(rate_times, flux, "o-", markersize=3)
    ax_flux.set_ylabel(r"$A_z$ at X-point")

    ax_rate.plot(times_c, rates, "s-", markersize=3, color="tomato")
    ax_rate.axhline(0, color="k", linewidth=0.5)
    ax_rate.set_ylabel(r"Reconnection rate $\mathrm{d}A_z/\mathrm{d}t$")
    ax_rate.set_xlabel("time")

    fig_rate.suptitle("Reconnection rate — Harris sheet")
    fig_rate.tight_layout()
    fig_rate.savefig("reconnection_rate.png", dpi=150)
    plt.show()

The script produces two figures:

- **Az_contours.png** — :math:`A_z` field lines at each snapshot with the
  primary X-point marked in red.
- **reconnection_rate.png** — time evolution of :math:`A_z` at the X-point
  (top) and the reconnection rate (bottom).

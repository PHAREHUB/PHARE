=======================
External magnetic field
=======================

**Goal:** give a run an external (background) magnetic field B0, a field you
prescribe and that PHARE does not evolve, and write it to disk.

.. warning::

   No solver uses B0 yet: declaring it does not change the plasma or the
   magnetic field of the run. PHARE computes it, keeps it up to date and
   writes it.


Minimal example
---------------

The block goes right after ``Simulation`` and before the model:

.. code-block:: python

    import pyphare.pharein as ph

    sim = ph.Simulation(
        # ... a 2D grid and time, as in the first simulation ...
    )

    ph.DipoleExternalField(position=(10.0, 10.0), moment=(0.0, 1.0), radius=1.0)

    # ... ph.MaxwellianFluidModel(...) and ph.ElectronModel(...), or ph.MHDModel(...) ...

    ph.ElectromagDiagnostics(quantity="B0", write_timestamps=[0.0])

The run writes ``EM_B0.h5``, read back with ``run.GetB0(0.0)`` (see
:doc:`/data_analysis/reading_data`).


What happens
------------

Both the :term:`hybrid model` and :term:`MHD` accept the block. Without it
there is no B0: it is neither allocated nor written, and asking for its
diagnostics is an error.

B0 is never given directly: PHARE takes a vector potential :math:`\mathbf{A}_0`
and computes :math:`\mathbf{B}_0 = \nabla\times\mathbf{A}_0` with the curl of
the grid, on every patch, ghost nodes included. The divergence of B0 on the
grid is then zero to round-off, whatever the potential.

B0 is computed

* at the start of the run, on every :term:`level`;
* on a level each time it is created or regridded: it is computed again from
  the potential there, not interpolated from the coarser level;
* after each step of a level, at the new time, when the potential depends on
  time;
* at a restart, from the block of the script that restarts: the block can be
  changed, added or removed between a run and its restart.


Two ways to give the field
--------------------------

**A dipole.** ``DipoleExternalField`` is a static magnetic dipole, in 2D or 3D:

.. code-block:: python

    ph.DipoleExternalField(position=(10.0, 10.0), moment=(0.0, 1.0), radius=1.0)

Inside ``radius`` around ``position`` the field is uniform, which removes the
singularity at the center; ``radius=0`` is the point dipole.

**Your own potential.** ``UserDefinedExternalField`` takes the components of
:math:`\mathbf{A}_0` as functions of the coordinates. A component you do not
give is zero:

.. code-block:: python

    def az(x, y):
        return 0.3 * y + 0.2 * x      # B0 = (0.3, -0.2, 0)

    ph.UserDefinedExternalField(a0z=az)

For a field that changes in time, the functions take the time as a last
argument, and the time derivative of the potential is given too:

.. code-block:: python

    def az(x, y, t):
        return x * y * (1.0 + t)

    def daz_dt(x, y, t):
        return x * y

    ph.UserDefinedExternalField(a0z=az, da0z_dt=daz_dt)

The functions receive NumPy arrays (all the points of a patch at once) and
return an array of the same shape.

B0 must carry no current: the solvers take the current from B1 only and drop
the terms of :math:`\nabla\times\mathbf{B}_0`. ``UserDefinedExternalField``
checks it when it is declared: on 10 points per direction across the domain,
and on 10 instants from the start to the end of the run for a time-dependent
potential, it computes :math:`\nabla\times\mathbf{B}_0` by finite differences
of the potential and requires it below
:math:`10^{-6}\,\max|\mathbf{A}_0|/L^2`, with :math:`L` the smallest extent of
the domain. In 1D this means a uniform B0, that is a linear potential.

To remove the singularity of a dipole inside the domain, make the field uniform
inside a radius, as ``DipoleExternalField`` does: :math:`\mathbf{m}\times\mathbf{r}/\max(r, R)^3`
in 3D (:math:`\max(r, R)^2` in 2D). A smoothed potential such as
:math:`\mathbf{m}\times\mathbf{r}/(r^2 + R^2)^{3/2}` carries current everywhere
and is rejected. The uniform core still carries a current on its surface
:math:`r = R`: a sample point within a few finite-difference steps of it
(about :math:`0.03\,L`) is rejected too.


Options
-------

``DipoleExternalField`` (all required, by keyword):

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Option
     - Meaning
   * - ``position``
     - center of the dipole, one component per dimension
   * - ``moment``
     - dipole moment, one component per dimension, not zero
   * - ``radius``
     - finite and >= 0; the field is uniform inside it, 0 is the point dipole

``UserDefinedExternalField`` (by keyword):

.. list-table::
   :header-rows: 1
   :widths: 25 55 20

   * - Option
     - Meaning
     - Default
   * - ``a0x``, ``a0y``, ``a0z``
     - components of the potential, functions of the coordinates
       (``a0z(x, y)`` in 2D), or of the coordinates and the time
       (``a0z(x, y, t)``). At least one is required. In 1D ``a0x`` does not
       contribute to the field and must stay ``None``
     - ``None`` (zero)
   * - ``da0x_dt``, ``da0y_dt``, ``da0z_dt``
     - time derivative of the potential, functions of the coordinates and the
       time. At least one is required when the potential depends on time; not
       accepted otherwise
     - ``None`` (zero)


Common mistakes
---------------

* ``RuntimeError: An external field must be declared before the model`` — move
  the block between ``Simulation`` and ``MaxwellianFluidModel`` /
  ``MHDModel``.
* ``RuntimeError: An external field is already declared`` — a run has one
  external field block; add the fields in one potential.
* ``RuntimeError: Error: 'B0' diagnostics require an external field`` —
  ``ElectromagDiagnostics(quantity="B0")`` without an external field block.
* ``ValueError: Error: a dipole external field makes no sense in 1D`` — use
  ``UserDefinedExternalField`` in 1D.
* ``ValueError: Error: a time-dependent potential requires one of da0x_dt,
  da0y_dt, da0z_dt`` — the potential takes the time, give its time derivative.
* ``ValueError: Error: da0x_dt, da0y_dt, da0z_dt given but the potential does
  not depend on time`` — add the time as last argument of the potential, or
  drop the derivative.
* ``ValueError: Error: the external field is not curl-free at t=...`` — the
  curl of the potential's B0 carries current at that time and point; fix the
  potential.
* ``ValueError: Error: the external field is not finite at t=..., at (...)`` or
  ``the derivatives of the external field are not finite`` — the potential is
  singular at or near that point; remove the singularity (see above).


Reference
---------

.. autoclass:: pyphare.pharein.DipoleExternalField

.. autoclass:: pyphare.pharein.UserDefinedExternalField

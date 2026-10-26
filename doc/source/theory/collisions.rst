
Ion-ion collisions
==================

By default PHARE solves the Vlasov equation for each ion population. Coulomb collisions
between ions can be added: the Landau-Fokker-Planck collision operator is then sampled by
binary Monte-Carlo collisions between the macro-ions of each cell, as in Takizuka & Abe (1977)
and Nanbu (1997). Electron-ion collisions are not modelled; the fluid electrons only feel a
resistivity, see the Ohm's law terms.

.. code-block:: python

    sim = ph.Simulation(
        ...,
        collisions={"nu0": ph.collision_nu0(n0_cm3=5.0, B0_nT=5.0, coulomb_log=20.0)},
    )


Normalization
-------------

In PHARE units (time in :math:`\Omega_p^{-1}`, velocities in :math:`V_A`, densities in
:math:`n_0`, masses in :math:`m_p`, charges in :math:`e`), every ion-ion collision rate is
proportional to the single dimensionless number

.. math::

   \nu_0 = \ln\Lambda\, \frac{e^4 n_0}{4\pi\varepsilon_0^2 m_p^2 V_A^3 \Omega_p}
         = \ln\Lambda\, \frac{(c/V_A)^4}{4\pi n_0 d_i^3}

which is the ``nu0`` key. ``collision_nu0`` computes it from :math:`n_0`, :math:`B_0` and the
Coulomb logarithm, taken constant. In the solar wind at 1 AU (:math:`n_0 = 5` cm\ :sup:`-3`,
:math:`B_0 = 5` nT), :math:`\nu_0/\ln\Lambda = 2.15\times 10^{-8}`: collisions are negligible
there, and ``nu0`` is usually set larger to model a collisional regime.


Algorithm
---------

Collisions act on the domain particles of each cell, between the particle push and the moment
deposit of the corrector step, so they change velocities only, every ``every`` time steps of
each AMR level, with that level's time step. Coarse cells covered by a finer level are collided
too: their particles feed the finer level's border and new refined regions.

1. Populations with the same mass and charge form one species group.
2. Particles are paired within each group (consecutive particles of a random order, an odd group
   starting with a three-particle triangle at half strength) and between each pair of groups
   (every particle of the larger group meets one of the smaller, used in turn).
3. A pair of relative speed :math:`u` and reduced mass :math:`m_{ab}` is deflected according to
   the scattering parameter

   .. math::

      s = \nu_0 \frac{(Z_a Z_b)^2\, n_\mathrm{eff}\, \Delta t}{m_{ab}^2\, u^3}

   with Nanbu's cumulative law (default, ``"nanbu"``, using the fit of Pérez et al. 2012),
   Takizuka & Abe's gaussian (``"takizuka_abe"``) or Bobylev's fixed angle (``"bobylev"``).
4. Macro-particles of unequal weights (from refinement, or populations of different densities):
   :math:`n_\mathrm{eff}` uses the larger weight of the pair and each partner is deflected with
   the probability :math:`w_\mathrm{other}/w_\mathrm{max}` (Higginson et al. 2020, in the form
   of Angus et al. 2025). This conserves momentum and energy on average only;
   ``exact_conservation`` (default on) then restores them in each cell (Angus et al. 2025).

The random draws are functions of the seed, the level, the time step, the cell and the
particles, so they do not depend on the domain decomposition, and a restart replays them.
Runs on different numbers of MPI ranks still differ at round-off in the moments, which changes
the pairing: they are different statistical realizations of the same plasma.


Accuracy and cost
-----------------

The relaxation rates converge to those of the Landau operator as :math:`\Delta t \to 0`. On the
isotropization of a bi-Maxwellian, the rate error is roughly :math:`0.23\sqrt{s_{th}}` for
Nanbu's law, with :math:`s_{th}` the value of :math:`s` at the thermal relative speed; about
1.5 times that for Takizuka & Abe, half for Bobylev, and twice that with weights spread over two
decades. A warning is printed at initialization above :math:`s_{th} = 0.02` (about 3%).
``every`` multiplies :math:`s`, so it divides the cost by ``every`` and multiplies the error by
about :math:`\sqrt{\mathtt{every}}`. During the run, a second warning is printed once if more
than 1% of the pairs reach :math:`s > 1`, which for Maxwellians happens around
:math:`s_{th} \approx 0.04`, or fail to restore their energy.

Collisions add noise to the fields: with :math:`\nu_0 = 0.1`, the energy of the field
fluctuations of a Maxwellian plasma grows by 5 to 25%, and the electric field near a coarse-fine
boundary by 13 to 20%, scaling like the particle noise. No heating from placing the collisions
after the push was detected; the expected effect (Vay et al. 2025, for explicit PIC) is about
:math:`10^{-4}` of the energy over the runs tried, below their noise.

A collision step costs about 60 ns per particle (90 ns with two species), whatever the
dimension: about 2.5 particle pushes in 1D, less than one in 2D and 3D. ``every`` larger than 1
divides that cost.


Limitations
-----------

The Coulomb logarithm is a constant. Collisions require a constant time step, and act on hybrid
levels only, not on MHD levels.


References
----------

- Takizuka T., Abe H. (1977), J. Comput. Phys. 25, 205, doi:10.1016/0021-9991(77)90099-7
- Nanbu K. (1997), Phys. Rev. E 55, 4642, doi:10.1103/PhysRevE.55.4642
- Bobylev A. V., Nanbu K. (2000), Phys. Rev. E 61, 4576, doi:10.1103/PhysRevE.61.4576
- Pérez F. et al. (2012), Phys. Plasmas 19, 083104, doi:10.1063/1.4742167
- Higginson D. P., Holod I., Link A. (2020), J. Comput. Phys. 413, 109450,
  doi:10.1016/j.jcp.2020.109450
- Angus J. R. et al. (2025), J. Comput. Phys. 531, 113927, doi:10.1016/j.jcp.2025.113927
- Vay J.-L. et al. (2025), Phys. Rev. E 111, 025306, doi:10.1103/PhysRevE.111.025306
- Wang C. et al. (2008), J. Comput. Phys. 227, 4308, doi:10.1016/j.jcp.2007.12.027

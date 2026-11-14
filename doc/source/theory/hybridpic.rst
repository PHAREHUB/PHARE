
====================
The Hybrid PIC model
====================

The Hybrid formalism consists in modeling the plasma as a
combination of constituents treated with different physical models.
This usually means that ions and electrons are treated differently.
A rather complete description of the different ways a code can be "hybrid"
is given in The Hybrid Multiscale Simulation Technology by A.S. Lipatov.
In astrophysical and space applications, the main application domains of
PHARE, "hybrid" usually means that ions are considered at the kinetic level
while electrons are considered as a fluid. This is the case for PHARE and
this is what we mean by "hybrid" on this page.


The ion equations
-----------------

The hybrid model consists in evolving in space :math:`\mathbf{r}` and time
:math:`t`, the velocity distribution function :math:`f_p` of each ion
populations p under the influence of the electric :math:`\mathbf{E}`
and magnetic field :math:`\mathbf{B}`.
This is done by solving the Vlasov equation for all ion populations when
collisions are negligible.


.. math::
   \frac{\partial f_p}{\partial t} + \mathbf{v}\cdot \frac{\partial f_p}{\partial \mathbf{r}} + \frac{q_p}{m_p}(\mathbf{\mathbf{E} + \mathbf{v}\times\mathbf{B}})\cdot \frac{\partial f_p}{\partial \mathbf{v}} = 0
   :label: eqvlasov


Having the new distribution everywhere at :math:`t+\Delta t`, it is easy
to calculate the ion moments as the sum of the moments
of all populations. Namely, for the ion density :math:`n_i` and bulk
velocity :math:`\mathbf{u_i}`


.. math::
   :nowrap:

   \begin{align}
   n_i(\mathbf{r},t) &= \sum_p \int f_p(\mathbf{r}, \mathbf{v}, t) d\mathbf{v} \\
   \mathbf{u}_i(\mathbf{r},t) &= \frac{1}{n_i}\sum_p \int \mathbf{v} f_p(\mathbf{r}, \mathbf{v}, t) d\mathbf{v}
   \end{align}


The electron momentum equation
------------------------------
What about the electrons? Remember? They are assumed to behave as a fluid.
This is wrong of course in collisionless systems since nothing makes
the density, velocity etc. of the electrons to depend on purely local physics
as collisions would in a "real" fluid. But that's an approximation the hybrid
formalism makes to simplify the physics (and make simulation lighter) compared
to the fully kinetic system and that is already a much more realistic way of
modeling the plasma and say, single fluid magnetohydrodynamics. Now there are subtleties.
The electron momentum equation is:


.. math::
     m_en_e \frac{d\mathbf{u_e}}{dt} = -\nabla\cdot \mathbf{P_e} - e n_e(\mathbf{E} + \mathbf{u_e}\times\mathbf{B})



Electromagnetic field equations
-------------------------------
"Treating electrons as a fluid", you probably think we solve that equation,
in contrast to the Vlasov equation we used for ions. Well not really...
But let's say we did. Now we would have to wonder where the magnetic field
and electric field would come from. For the magnetic field, the answer is easy.
We just use the Maxwell-Faraday equation:


.. math::
   \frac{\partial \mathbf{B}}{\partial t} = -\nabla\times\mathbf{E}



What about the electric field now? There is all the trick of Hybrid codes.
We actually do not solve the electron momentum equation directly to get the
new electron fluid momentum. Instead we make assumptions on the electron
fluid, and use that momentum equation to calculate the electric field !
Thus, the momentum equation is re-written:


.. math::
   \mathbf{E} = -\mathbf{u_e}\times\mathbf{B} - \frac{1}{en_e}\nabla\cdot \mathbf{P_e}  +\frac{m_e}{e}\frac{d\mathbf{u_e}}{dt}



Quasineutrality
---------------

At this point, the equation for the electric field still has unknowns.
The most obvious perhaps is :math:`n_e` the electron particle density.
This is where the hybrid formalism makes the assumption that at the scale
we solve the equations, the plasma is quasineutral (:term:`quasineutrality`): the electron charge density
balances the total ion charge density,

.. math::
   e\,n_e = \sum_p q_p n_p \equiv e\,n

so we have only one variable :math:`n`, the plasma density, and it is known from the ions.
Quasineutrality also gives the electron bulk velocity. The current density is the sum of
the ion and electron currents, :math:`\mathbf{j} = \sum_p q_p n_p \mathbf{u}_p - e n \mathbf{u_e}`, so

.. math::
   \mathbf{u_e} = \frac{1}{en}\left(\sum_p q_p n_p \mathbf{u}_p - \mathbf{j}\right)

where :math:`q_p n_p \mathbf{u}_p` is the charge flux of population :math:`p`. This is what
PHARE computes. With a single population of protons (:math:`q_p = e`, :math:`n_p = n`) it
reduces to the familiar

.. math::
   \mathbf{u_e} = \mathbf{u_i} - \frac{\mathbf{j}}{en}

With several populations the ion bulk velocity :math:`\mathbf{u_i}` is weighted by mass,
not by charge, so this short form must not be used.

The total current density is obtained from the static Maxwell-Ampere equation,
meaning we neglect the displacement current:


.. math::
   \mu_0 \mathbf{j} = \nabla\times\mathbf{B}


The electric field is now equal to

.. math::

   \mathbf{E} = -\mathbf{u_e}\times\mathbf{B} - \frac{1}{en}\nabla\cdot \mathbf{P_e}  +\frac{m_e}{e}\frac{d\mathbf{u_e}}{dt}




Massless electrons
------------------

The next assumption usually made in Hybrid codes, that is also made in PHARE,
is that the spatial and time scales at which we are interested in are much larger
and longer that the scales at which the electron bulk inertia matters. The electrons
being so light compare to even protons, that it is mostly ok to neglect the last term of the equation above,
which now reads:



.. math::

    \mathbf{E} = -\mathbf{u_e}\times\mathbf{B} - \frac{1}{en}\nabla\cdot \mathbf{P_e}



Electron closure
----------------

Since we do not have an electron distribution function in hand, the pressure
is not known a priori. Hybrid codes thus have to come with a so-called closure
equation which role is to give us the pressure
everywhere at time t, based on some assumption on the system. Usually,
unless in very specific conditions, there is no rigorous way of getting such equation
and most hybrid code assume a closure that is "reasonable" and above all
"simple" to calculate.

Perhaps the simplest and most used electron closure is the isothermal one.
This simply say that the electron pressure :math:`P_e` is given by the product of
the density by some scalar constant that we call "the electron temperature".


.. math::

    P_e= nT_e


Dissipative terms
-----------------

Using above equations to calculate the electric field would result in
current sheets to collapse at grid scale in the  absence of an intrinsic
dissipation scale in the system. Two ways are typically employed in
Hybrid codes to include such a dissipation. Joule resistivity well known to
be used already in MHD codes. It is a simple term :math:`\eta \mathbf{j}` to
add on the right hand side of the electric field equation. This term adds
diffusion of magnetic flux. However there is no scale at which this terms
dominate over the electron ideal term :math:`-\mathbf{u_e}\times\mathbf{B}`,
unless :math:`\eta` is so large that ion scale structures are diffused away too.

Another term that can be employed is the so-called hyper-resistivity
(sometimes called hyper-viscosity) that takes the form :math:`-\nu\nabla^2\mathbf{j}`
In contrast to classical resistivity, this terms (due to the second order
derivative) comes with an intrinsic scale at which it is dominant over electron
convection term and efficiently adds sub-ion scale dissipation.

PHARE include these two terms and the electric field is obtained via :


.. math::

    \mathbf{E} = -\mathbf{u_e}\times\mathbf{B} - \frac{1}{en}\nabla P_e +\eta\mathbf{j} - \nu\nabla^2\mathbf{j}

By default :math:`\nu` is a constant. With ``hyper_mode="spatial"`` it is multiplied by a
factor that grows with :math:`|\mathbf{B}|/n` and follows the square of the mesh size on
refined levels (on level :math:`\ell` it never drops below :math:`\nu/4^{\ell}`), so that the dissipation stays at the
grid scale everywhere (see :doc:`/user_guide/electrons`).


Numerical scheme
----------------

PHARE solves the ion equations with the Particle-In-Cell method (:doc:`pic`) and the
field equations on a staggered grid (:doc:`/numerics/grid`). One time step closes the loop
between them:

.. plot::
   :include-source: False
   :caption: The loop of one step, in 1D. 1 · each particle adds to the moments (density,
      flux) of its nearest nodes. 2 · Ohm's and Faraday's laws give E and B on the grid (here
      the transverse components: Ey, Ez on nodes, By, Bz at cell centres). 3 · E and B, read
      at each particle, move it.

   from phare_draw import *

   T = ["Deposit → moments", "Ohm, Faraday → fields", "Push"]
   fig, ax = strip(3, titles=T, phone=True, aspect=False, height=1.3)
   parts = [1.3, 2.55, 3.35]
   for k, a in enumerate(ax):
       a.plot([0, 4], [0, 0], color=C["light"], lw=1)
       for x in range(5):
           mark(a, "primal", x, 0, color=C["ion"] if k == 0 else C["E"], size=8)
       if k == 1:
           for x in range(4):
               mark(a, "dual", x + 0.5, 0, color=C["B"], size=7)
           text(a, 0, 0.42, "Ey, Ez", color=C["E"], bold=True, ha="left"); text(a, 0.5, -0.42, "By, Bz", color=C["B"], bold=True, ha="left")
       if k != 1:
           particles(a, [(p, 0.35) for p in parts], size=40)
       if k == 0:
           for p in parts:
               for x in (int(p), int(p) + 1):
                   arrow(a, (p, 0.3), (x, 0.08), color=C["grey"], lw=1, ms=8)
       if k == 2:
           for p, d in zip(parts, (0.35, -0.25, 0.3)):
               arrow(a, (p, 0.55), (p + d, 0.55), color=C["ion"], lw=1.8, ms=10)
       a.set_xlim(-0.3, 4.3); a.set_ylim(-0.6, 0.75)
   finish(fig, ax)

The order of these operations within a step (a :term:`predictor-predictor-corrector <PPC>` scheme) is on
:doc:`/numerics/time_integration`; the particle push and the deposit are on
:doc:`/numerics/particles`.





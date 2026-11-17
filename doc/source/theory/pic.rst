
==============================
The Particle-In-Cell formalism
==============================


There are two ways to solve the Vlasov equation for ion populations. It can be
solved calculating eulerian derivatives, i.e. discretizing velocity and spatial
dimensions and solving the equation at those fixed locations. A hybrid code that solves its
ions this way is called a "Vlasov hybrid code" (PHARE is not one). This approach is
generally complex and requires lots of computational
resources. The other way consists in adopting a Lagrangian viewpoint. That is,
cutting the initial distribution function in N weighted bins and follow the dynamics
of those chunks in phase space. The little pieces of distributions are called
"macro-particles". Decomposing the distribution function of the population into the
contribution of macro-particles in the following way is the base of the so-called
"Particle-in-Cell" (PIC) method.


.. math::
   f_p(\mathbf{r}, \mathbf{v}, t) = \sum_m^{N_p} w_m S(\mathbf{r} - \mathbf{r}_m(t))\delta(\mathbf{v}-\mathbf{v}_m(t))



where :math:`\mathbf{r}_m` and :math:`\mathbf{v}_m` are the position and velocity
of the :math:`m_{th}` macro-particle. :math:`w_m` represents the weight of that
macro-particle, i.e. how much it counts in the evaluation of :math:`f_p`.
:math:`\delta` is the Dirac function, which says that a macro-particle represent
a specific velocity in the distribution. In contrast, the function :math:`S` is
a finite support function representing the "shape" of the macro-particle in the
spatial dimension. This function tells us how far from the macro-particle a local
distribution sees its influence. In PHARE we use b-spline functions to model :math:`S`.
PHARE uses b-splines of the first, second and third order.
The higher the order  the further a macro-particle influences the distribution,
but the longer it takes to compute it.

The shape is chosen with ``interp_order`` (1, 2 or 3) in the ``Simulation`` block. It is
also the weight with which a particle contributes to the grid nodes around it
(:term:`deposition`) and with which it reads the fields there.

.. plot::
   :include-source: False
   :caption: The particle shapes :math:`S(x)` for ``interp_order`` 1, 2 and 3 (B-splines
      of degree 1, 2, 3), for a particle at :math:`x = 0`. They cover 2, 3 and 4 cells.
      Each has unit area, so a particle always deposits its full weight.

   import numpy as np
   import matplotlib.pyplot as plt

   def bspline(x, order):
       a = np.abs(x)
       if order == 1:
           return np.clip(1 - a, 0, None)
       if order == 2:
           return np.where(a < 0.5, 0.75 - a**2,
                           np.where(a < 1.5, 0.5 * (1.5 - a)**2, 0.0))
       return np.where(a < 1, 2 / 3 - a**2 + a**3 / 2,
                       np.where(a < 2, (2 - a)**3 / 6, 0.0))

   x = np.linspace(-2.5, 2.5, 1001)
   fig, ax = plt.subplots(figsize=(5, 3.4))
   for order, c in zip((1, 2, 3), ("#2a7ab9", "#e07a1f", "#3a9d5d")):
       ax.plot(x, bspline(x, order), color=c, lw=2.5, label=f"order {order}")
   ax.set_xlabel("distance to the particle (cells)", fontsize=14)
   ax.set_xticks(range(-2, 3))
   ax.set_yticks([])
   ax.tick_params(labelsize=14)
   ax.legend(fontsize=13, frameon=False)
   for side in ("top", "right", "left"):
       ax.spines[side].set_visible(False)
   fig.tight_layout()



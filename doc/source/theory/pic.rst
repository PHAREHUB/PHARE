
==============================
The Particle-In-Cell formalism
==============================


There are two ways to solve the Vlasov equation for ion populations. It can be
solved calculating eulerian derivatives, i.e. discretizing velocity and spatial
dimensions and solving the equation at those fixed locations. This is called a
"Vlasov Hybrid code". It is generally complex and require lots of computational
resources. The other way consists in adopting a Lagrangian viewpoint. That is,
cutting the initial distribution function in N weighted bins and follow the dynamics
of those chunks in phase space. The little pieces of distributions are called
"macro-particles". Decomposing the distribution function of the population into the
contribution of macro-particles in the following way is the base of the so-called
"Particle-in-Cell" (PIC) method.


.. math::
   \begin{equation}
   f_p(\mathbf{r}, \mathbf{v}, t) = \sum_m^{N_p} w_m S(\mathbf{r} - \mathbf{r}_m(t))\delta(\mathbf{v}-\mathbf{v}_m(t))
   \end{equation}



where :math:`\mathbf{r}_m` and :math:`\mathbf{v}_m` are the position and velocity
of the :math:`m_{th}` macro-particle. :math:`w_m` represents the weight of that
macro-particle, i.e. how much it counts in the evaluation of :math:`f_p`.
:math:`\delta` is the Dirac function, which says that a macro-particle represent
a specific velocity in the distribution. In contrast, the function :math:`S` is
a finite support function representing the "shape" of the macro-particle in the
spatial dimension. This function tells us how far from the macro-particle a local
distribution sees its influence. In PHARE we use b-splinefunctions to model :math:`S`.
PHARE uses b-splines of the first, second and third order.
The higher the order  the further a macro-particle influences the distribution,
but the longer it takes to compute it.



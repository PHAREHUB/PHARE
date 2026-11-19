========================================================
Hybrid time integration: predictor, predictor, corrector
========================================================

**Goal:** know how the hybrid solver advances the fields and the particles over one time
step (:term:`PPC` for short). How to choose the time step is on :doc:`time_integration`.

The difficulty in a hybrid code is that the electric field depends on the ion moments at
the *end* of the step, which are only known after the particles have moved, and the
particles move in that same field. PHARE solves this with two predictions and one
correction. Starting from :math:`\mathbf{B}^n`, :math:`\mathbf{E}^n` and the ion moments
at time :math:`n`:

.. plot::
   :include-source: False
   :caption: Red: fields; green: particles and ion moments. Each stage starts again from Bⁿ and
      advances with a mid-step (averaged) field; only the second push is kept.

   from phare_draw import *

   T = ["Predictor 1, then a trial push",
        "Predictor 2, then the real push",
        "Corrector"]
   rows = [("B*, E*  (Ohm with nⁿ)", "ñⁿ⁺¹: predicted, particles discarded", "Faraday with Eⁿ"),
           ("B**, E**  (Ohm with ñⁿ⁺¹)", "nⁿ⁺¹: particles kept", "Faraday with average 1"),
           ("Bⁿ⁺¹, Eⁿ⁺¹  (Ohm with nⁿ⁺¹)", None, "Faraday with average 2")]
   fig, ax = strip(3, titles=T, aspect=False, height=1.6)
   for k, (a, (flab, push, far)) in enumerate(zip(ax, rows)):
       for x, lab in [(0, "tⁿ"), (5, "tⁿ+½"), (10, "tⁿ⁺¹")]:
           a.plot([x, x], [-0.95, -0.8], color=C["grey"], lw=1)
           text(a, x, -1.2, lab, color=C["grey"])
       a.plot([0, 10], [-0.875, -0.875], color=C["light"], lw=1)
       mark(a, "primal", 0, 1, color=C["grey"]); text(a, 0, 0.65, "Bⁿ", color=C["grey"])
       arrow(a, (0.3, 1), (9.6, 1), color=C["B"], lw=2.2)
       mark(a, "primal", 10, 1, color=C["B"]); text(a, 10, 1.45, flab, color=C["B"], bold=True, ha="right")
       text(a, 0.3, 1.25, far, color=C["B"], ha="left", va="bottom", fs=10)
       if push:
           mark(a, "dual", 5, 1, color=C["E"]); text(a, 5, 0.6, f"average {k + 1} of tⁿ and tⁿ⁺¹ fields", color=C["E"])
           mark(a, "primal", 0, 0, color=C["ion"])
           arrow(a, (0.3, 0), (9.6, 0), color=C["ion"], lw=2.2)
           mark(a, "primal", 10, 0, color=C["ion"])
           text(a, 10, -0.4, push, color=C["ion"], ha="right")
           text(a, 0.3, -0.4, "push with the average", color=C["ion"], ha="left")
       a.set_xlim(-0.6, 10.6); a.set_ylim(-1.45, 1.8)
   finish(fig, ax)

In equations, with :math:`\mathcal{O}(\mathbf{B}, \text{moments})` the electric field from
Ohm's law (:doc:`/theory/hybridpic`):

.. math::

   \begin{aligned}
   \mathbf{B}^{*} &= \mathbf{B}^n - \Delta t\,\nabla\times\mathbf{E}^n, &
   \mathbf{E}^{*} &= \mathcal{O}(\mathbf{B}^{*}, n^n, \mathbf{F}^n) \\
   \text{push with } &\tfrac12(\mathbf{E}^n + \mathbf{E}^{*}),\ \tfrac12(\mathbf{B}^n + \mathbf{B}^{*})
   & &\rightarrow\ \tilde n^{n+1}, \tilde{\mathbf{F}}^{n+1} \\
   \mathbf{B}^{**} &= \mathbf{B}^n - \tfrac{\Delta t}{2}\nabla\times(\mathbf{E}^n + \mathbf{E}^{*}), &
   \mathbf{E}^{**} &= \mathcal{O}(\mathbf{B}^{**}, \tilde n^{n+1}, \tilde{\mathbf{F}}^{n+1}) \\
   \text{push with } &\tfrac12(\mathbf{E}^n + \mathbf{E}^{**}),\ \tfrac12(\mathbf{B}^n + \mathbf{B}^{**})
   & &\rightarrow\ n^{n+1}, \mathbf{F}^{n+1} \\
   \mathbf{B}^{n+1} &= \mathbf{B}^n - \tfrac{\Delta t}{2}\nabla\times(\mathbf{E}^n + \mathbf{E}^{**}), &
   \mathbf{E}^{n+1} &= \mathcal{O}(\mathbf{B}^{n+1}, n^{n+1}, \mathbf{F}^{n+1})
   \end{aligned}

Here :math:`n` is the ion charge density and :math:`\mathbf{F}` the ion charge flux; Ohm's law
uses the version of it deposited on the positions of E, so that it conserves energy
(:doc:`hybrid_ohm`). The
first push only serves to predict the moments; its particles are thrown away. The second
push is the real one.

.. seealso:: How this works inside: :doc:`/dev/ppc_order`.

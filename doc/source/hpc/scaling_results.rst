===============
Scaling results
===============

**Goal:** show how PHARE's run time changes with the number of MPI ranks,
for reference setups you can reproduce.

.. admonition:: In preparation
   :class: in-preparation

   This page will give, for each reference setup: the input script, the
   machine and compiler, and the measured results below.

Strong scaling
--------------

Same problem, more ranks: how much faster does it run?

.. admonition:: In preparation
   :class: in-preparation

   Time per step versus number of ranks, with the ideal line.

Weak scaling
------------

Problem size grows with the number of ranks: does the time per step stay
constant?

.. admonition:: In preparation
   :class: in-preparation

   Time per step versus number of ranks at fixed work per rank.

Where the time goes
-------------------

.. admonition:: In preparation
   :class: in-preparation

   Breakdown of a step (particle push, field solve, ghost exchanges,
   regrid, diagnostics) from the scope timers of
   :doc:`measuring_performance`.

How to reproduce
----------------

.. admonition:: In preparation
   :class: in-preparation

   Scripts, job files and the commands to plot the results.

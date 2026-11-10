# Adaptive Mesh Refinement

## Patch based approach

## Recursive time integration

## Regridding

Regridding is the process of changing the geometry of an existing mesh level in the AMR hierarchy,
by changing the number and geometry of some patches of that level. This typically occurs when the
needs for refinement evolve with the solution over time.

Concretely, regridding means removing an entire level from the AMR hierarchy and replacing it by a
new one. The so-called "new" level is initialized by copying values from the old level wherever they
overlap, and by refining values from coarser levels where they do not.

Technically, regridding involves different methods depending on whether particle data or field data
is refined. The specific operation performed also depends on the type of field being refined: for
instance, the magnetic field refinement must preserve the divergence-free character of the field.

## Field refinement

## Particle refinement

## Field coarsening

Coarsening, also sometimes called "restriction", is the process of projecting the solution existing
on a given level onto the region it overlaps on the next coarser level. Coarsening is used so that
the solution on the coarse level "feels" the fine solution in the overlapped region, which is
assumed to be of better quality.

## Fields at level boundaries

## Particle at level boundaries

(tagging-for-refinement)=

## Tagging for refinement

When refined regions are not given by the user, they are found from the solution itself. At
initialization and whenever the hierarchy is regridded, a criterion is evaluated in every cell of a
level. The cells where it is met are tagged, then grown by a few cells and
{ref}`clustered <clustering>` into the patches of the next finer level. This section defines the
criteria, whose syntax is described in {doc}`/usage/amr_tagging`.

### Common principles

A criterion computes an indicator $I(\vb{x})$ in each cell from the values of one quantity $q$, and
the cell is tagged if $I(\vb{x}) > \varepsilon$ for at least one of the tagged quantities,
$\varepsilon_q$ being the threshold given for quantity $q$.

Below, $\vb{u}(\vb{x})$ denotes the vector of the cell-centred components of the quantity in the
cell centered on $\vb{x}$ (a single value for a scalar quantity) and $\left\|\cdot\right\|$ its
Euclidean norm, the absolute value for a scalar quantity. $u(\vb{x})$ denotes the cell-centred value
of one component, and $\vu{e}_p$ the unit offset along Cartesian direction $p$.

An indicator reads the values up to a few cells away from the tagged cell, its stencil reach. When
this reach exceeds the number of ghost cells of the fields, the corresponding outermost layers of
cells of each patch are not tagged. Debug builds print a warning at startup when this happens.

### Löhner's criterion

Löhner's error estimator[^lohner1987] measures the change of slope from one cell to the next,
normalized by the slopes themselves. A noise filter, weighted by $\epsilon_r$, ignores the
variations that are small compared with the local magnitude of the quantity, and $\epsilon_a$ is an
absolute floor, in code units. $\epsilon_r$ and $\epsilon_a$ are the parameters `reltol` (default
0.01) and `abstol` (default $10^{-8}$). PHARE provides the two multidimensional forms of the
literature, which both lie between 0 and $1/(1 + \epsilon_r)$.

#### Per direction

`lohner-per-direction` is the form of PLUTO[^mignone2012]. Along each direction $d$,

$$
\begin{aligned}
n_d &= \left\|\vb{u}(\vb{x} + \vu{e}_d) - 2\vb{u}(\vb{x}) + \vb{u}(\vb{x} - \vu{e}_d)\right\|, \\
D_d &= \left\|\vb{u}(\vb{x} + \vu{e}_d) - \vb{u}(\vb{x})\right\| + \left\|\vb{u}(\vb{x}) - \vb{u}(\vb{x} - \vu{e}_d)\right\| \\
    &\quad + \epsilon_r \left(\left\|\vb{u}(\vb{x} + \vu{e}_d)\right\| + 2\left\|\vb{u}(\vb{x})\right\|
      + \left\|\vb{u}(\vb{x} - \vu{e}_d)\right\|\right) + \epsilon_a
\end{aligned}
$$

and the directions are combined in a single ratio,

$$
I = \sqrt{\frac{\sum_d n_d^2}{\sum_d D_d^2}}
$$

The stencil reach is one cell.

#### With cross derivatives

`lohner-cross-derivatives` is the form of FLASH[^flash4ug] [^fryxell2000], which includes all the
cross derivatives and reads points two cells apart. With the slope over two cells
$\Delta_p \vb{u}(\vb{y}) = \vb{u}(\vb{y} + \vu{e}_p) - \vb{u}(\vb{y} - \vu{e}_p)$, every ordered
pair of directions $(p, q)$, $p = q$ included, contributes

$$
\begin{aligned}
n_{pq} &= \left\|\Delta_p \vb{u}(\vb{x} + \vu{e}_q) - \Delta_p \vb{u}(\vb{x} - \vu{e}_q)\right\|, \\
D_{pq} &= \left\|\Delta_p \vb{u}(\vb{x} + \vu{e}_q)\right\| + \left\|\Delta_p \vb{u}(\vb{x} - \vu{e}_q)\right\|
          + \epsilon_r S_{pq} + \epsilon_a
\end{aligned}
$$

where $S_{pq}$ is the sum of $\left\|\vb{u}\right\|$ over the four points
$\vb{x} \pm \vu{e}_q \pm \vu{e}_p$. The pairs are weighted by the mesh sizes $\Delta x_p$ of the
patch,

$$
I = \sqrt{\frac{\sum_{p,q} \left(w_{pq}\, n_{pq}\right)^2}{\sum_{p,q} \left(w_{pq}\, D_{pq}\right)^2}},
\qquad w_{pq} = \frac{1}{\Delta x_p \Delta x_q}
$$

The stencil reach is two cells.

### Default criterion

In 2D and 3D, along each direction $d$, and considering a scalar quantity $u$.

$$
I_d = \frac{\left|u(\vb{x} + 2\vu{e}_d) - u(\vb{x})\right|}{1 + \left|u(\vb{x} + \vu{e}_d) - u(\vb{x})\right|}
$$

and the largest value over the directions is kept, with a stencil reach of two cells. Vector
components are treated as independant scalar quantities in 2D and 3D.

In 1D, with $\bar u(\vb{x})$ the average of the five values $u(\vb{x} + k\vu{e}_x)$,
$k = -2, \dots, 2$,

$$
I = \frac{\left|\bar u(\vb{x} + \vu{e}_x) - \bar u(\vb{x})\right|}{1 + \left|\bar u(\vb{x})\right|}
  = \frac{\left|u(\vb{x} + 3\vu{e}_x) - u(\vb{x} - 2\vu{e}_x)\right|}{5\,\left(1 + \left|\bar u(\vb{x})\right|\right)}
$$

with a stencil reach of three cells. For a vector quantity in 1D, the indicators $I^c$ of its
components are combined into $\sqrt{\sum_c (I^c)^2}$, which is compared with $\varepsilon$.

(clustering)=

## Clustering

(tiling)=

### Tiling

(berger-rigoutsos)=

### Berger-Rigoutsos

[^lohner1987]: R. Löhner, *An adaptive finite element scheme for transient problems in CFD*,
    Computer Methods in Applied Mechanics and Engineering 61 (1987),
    <https://doi.org/10.1016/0045-7825(87)90098-3>

[^mignone2012]: A. Mignone et al., *The PLUTO code for adaptive mesh computations in astrophysical
    fluid dynamics*, The Astrophysical Journal Supplement Series 198 (2012), 7, eq.
    50-51, <https://arxiv.org/abs/1110.0740>

[^flash4ug]: FLASH4 user guide, §8.6.3.1 Refinement Criteria, eq. 8.3-8.4,
    <https://flash.rochester.edu/site/flashcode/user_support/flash4_ug_4p62/node59.html>

[^fryxell2000]: B. Fryxell et al., *FLASH: an adaptive mesh hydrodynamics code for modeling
    astrophysical thermonuclear flashes*, The Astrophysical Journal Supplement Series
    131 (2000), 273-334, <https://iopscience.iop.org/article/10.1086/317361/fulltext/>

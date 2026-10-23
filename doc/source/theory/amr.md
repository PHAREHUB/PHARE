# Adaptive Mesh Refinement

## Patch based approach

## Recursive time integration

## Regridding

Regridding is the process of changing the geometry of an existing mesh level in the AMR hierarchy,
by changing the number and geometry of some patches of that level. This typically occurs when the
needs for refinement evolve with the solution over time.

Concretely, regridding means removing an entire level from the AMR hierarchy and replacing it by a
new one. The so-called "new" level is initialized by copying values from the old level wherever
they overlap, and by refining values from coarser levels where they do not.

Technically, regridding involves different methods depending on whether particle data or field
data is refined. The specific operation performed also depends on the type of field being refined:
for instance, the magnetic field refinement must preserve the divergence-free character of the
field.

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

A criterion computes an indicator $I$ in each cell from the values of one quantity, and the cell is
tagged if $I > \tau$ for at least one of the tagged quantities, $\tau$ being the threshold given for
that quantity. The indicator of a vector quantity is computed on each component separately, and the
cell is tagged if any of them exceeds $\tau$. The default criterion in 1D is the exception: it
combines the indicators of the components into one value before comparing it with $\tau$.

The criteria read cell-centred values. A quantity that is not cell-centred, such as a magnetic field
component stored on cell faces, is averaged over the two neighbouring values in each direction
where it is not cell-centred.

In 2D and 3D, the indicator is computed along each direction separately, from the values of the
cells aligned with the tagged cell, and the largest of these values is kept. Below, $a_i$ denotes
the cell-centred value of one component in cell $i$ along the direction considered.

An indicator reads the values up to a few cells away from the tagged cell, its stencil reach. When
this reach exceeds the number of ghost cells of the fields, the corresponding outermost layers of
cells of each patch are not tagged. Debug builds print a warning at startup when this happens.

### Löhner's criterion

This criterion is the error estimator of Löhner[^lohner1987], in the form used by the FLASH
code[^fryxell2000]. It is computed along each direction separately, without the cross derivatives
that the multidimensional form of FLASH combines with the other directions into a single ratio:

$$
I_i = \frac{|a_{i+1} - 2a_i + a_{i-1}|}
           {|a_{i+1} - a_i| + |a_i - a_{i-1}|
            + \epsilon_r \left(|a_{i+1}| + 2|a_i| + |a_{i-1}|\right) + \epsilon_a}
$$

where $\epsilon_r$ and $\epsilon_a$ are the parameters `reltol` (default 0.02) and `abstol`
(default $10^{-8}$). The numerator is the change of slope from one cell to the next, normalized by
the slopes themselves. The term in $\epsilon_r$ is a noise filter, which ignores the variations
that are small compared with the local value of the quantity, and $\epsilon_a$ is an absolute
floor, in code units. The stencil reach is one cell.

### Default criterion

In 2D and 3D, along each direction,

$$
I_i = \frac{|a_{i+2} - a_i|}{1 + |a_{i+1} - a_i|}
$$

with a stencil reach of two cells.

In 1D, with $\bar a_i$ the average of the five values $a_{i-2}$ to $a_{i+2}$,

$$
I_i = \frac{|\bar a_{i+1} - \bar a_i|}{1 + |\bar a_i|}
    = \frac{|a_{i+3} - a_{i-2}|}{5\,(1 + |\bar a_i|)}
$$

with a stencil reach of three cells. For a vector quantity, the indicators $I^c$ of its components
are combined into $\sqrt{\sum_c (I^c)^2}$, which is compared with $\tau$.

[^lohner1987]: R. Löhner, *An adaptive finite element scheme for transient problems in CFD*,
    Computer Methods in Applied Mechanics and Engineering 61 (1987),
    <https://doi.org/10.1016/0045-7825(87)90098-3>

[^fryxell2000]: B. Fryxell et al., *FLASH: an adaptive mesh hydrodynamics code for modeling
    astrophysical thermonuclear flashes*, The Astrophysical Journal Supplement Series 131 (2000),
    273-334

(clustering)=

## Clustering

(tiling)=

### Tiling

(berger-rigoutsos)=

### Berger-Rigoutsos

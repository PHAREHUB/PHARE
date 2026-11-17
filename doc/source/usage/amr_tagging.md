# AMR tagging

With `refinement="tagging"`, refined regions are not given by the user but found from the solution
itself. At initialization and whenever the hierarchy is regridded, a criterion is evaluated in every
cell of a level. The cells where it is met are tagged, grown by `tag_buffer` cells, and clustered
into the patches of the next finer level. Refined regions therefore follow the structures of the
solution as it evolves.

This page gives the syntax. The criteria are defined in {ref}`tagging-for-refinement`.

## Enabling tagging

```python
import pyphare.pharein as ph

ph.Simulation(
    # ...
    refinement="tagging",
    max_nbr_levels=3,
    tag_buffer=1,
    tagging_threshold=0.1,
)
```

- `refinement="tagging"` selects tagging. The default, `"boxes"`, uses the fixed regions given by
  `refinement_boxes` instead.
- `max_nbr_levels` must be given: it defaults to 1, that is no refinement at all.
- `tag_buffer` (default 1) is the number of cells by which each tagged cell is grown in every
  direction before clustering.
- `tagging_threshold` (default 0.1) is the threshold applied to the magnetic field when the
  `tagging` keyword below gives no quantities.

Without the `tagging` keyword, the default criterion is applied to the magnetic field with threshold
`tagging_threshold`. The other keywords shaping the hierarchy (`clustering`, `nesting_buffer`,
`smallest_patch_size`, `largest_patch_size`) are described with
{py:class}`~pyphare.pharein.Simulation` in {doc}`simulation`.

## Choosing the criterion and the quantities

The `tagging` keyword is a dict whose keys are all optional:

```python
ph.Simulation(
    # ...
    refinement="tagging",
    max_nbr_levels=3,
    tagging={
        "method": "lohner-per-direction",
        "quantities": {"B": {"threshold": 0.5, "reltol": 0.02}, "rho": 0.5},
    },
)
```

- `method` is the criterion: `"default"` (the default), `"lohner-per-direction"` (Löhner's estimator
  in the form of PLUTO), `"lohner-cross-derivatives"` (Löhner's estimator in the form of FLASH, with
  cross derivatives) or `"gradient-scale"` (cell size over the gradient length of the quantity).
- `quantities` maps each tagged quantity to its threshold and the parameters of the method for this
  quantity, as a dict `{"threshold": ..., <parameter>: ...}`. A number alone is the threshold, the
  parameters then taking their defaults. A cell is tagged if the criterion is met for at least one
  quantity, each quantity being compared to its own threshold. Without `quantities`, the magnetic
  field is tagged with threshold `tagging_threshold`, except with `"gradient-scale"`, whose required
  `ref_magnitude` needs explicit quantities. A list of `(name, threshold or dict)` pairs is also
  accepted, a name given twice being an error.

The parameters of each method, given per quantity next to its threshold:

| Method              | Parameter       | Default   | Meaning                                                |
| ------------------- | --------------- | --------- | ------------------------------------------------------ |
| `"default"`         |                 |           | no parameter                                           |
| both Löhner methods | `reltol`        | $0.01$    | weight of the noise filter, $\geq 0$                   |
| both Löhner methods | `ref_magnitude` | $10^{-8}$ | reference magnitude added to the denominator, $\geq 0$ |
| `"gradient-scale"`  | `ref_magnitude` | required  | reference magnitude added to the magnitude, $> 0$      |

`ref_magnitude` is in the units of the quantity. An unknown parameter is an error.

Giving `tagging` without `refinement="tagging"` is an error.

### Taggable quantities

| Name                      | Hybrid | MHD | Quantity                              |
| ------------------------- | ------ | --- | ------------------------------------- |
| `B`                       | yes    | yes | magnetic field                        |
| `Bx`, `By`, `Bz`          | yes    | yes | one component of the magnetic field   |
| `rho` or `density`        | no     | yes | mass density                          |
| `rhoV` or `momentum`      | no     | yes | momentum density                      |
| `rhoVx`, `rhoVy`, `rhoVz` | no     | yes | one component of the momentum density |
| `Etot`                    | no     | yes | total energy density                  |

The electric field, the current density, the velocities, the pressures and the hybrid ion densities
are not available for tagging yet: they are rejected with an error.

### Choosing a threshold

The criteria are not on the same scale, so a threshold does not carry over from one method to
another.

- `"lohner-per-direction"` and `"lohner-cross-derivatives"` are dimensionless and lie between 0 and
  $1/(1 + \mathrm{reltol})$. Their thresholds are typically a few tenths.
- `"default"` measures the variations of the quantity in code units, absolutely when they are small
  compared with 1 and relatively when they are large.
- `"gradient-scale"` is the cell size over the gradient length of the quantity, so its threshold is
  the target ratio $\Delta x / L$. A smooth profile of length $L$ gives
  $\Delta x / L \cdot \|\vb{u}\| / (\|\vb{u}\| + u_\mathrm{ref})$, $u_\mathrm{ref}$ being
  `ref_magnitude`, and a jump $\Delta u$ over one cell gives
  $\Delta u / (2(\|\vb{u}\| + u_\mathrm{ref}))$ on both sides of it. Where the quantity is small
  compared with `ref_magnitude`, its variations are measured relative to `ref_magnitude` instead of
  its magnitude, which keeps noise around nulls from being refined.

{ref}`tagging-for-refinement` provides the definition of the criteria.

## Examples

Default criterion on the magnetic field, as when `tagging` is not given:

```python
ph.Simulation(
    # ...
    refinement="tagging",
    max_nbr_levels=3,
    tagging_threshold=0.1,
)
```

Löhner's criterion on the density and on one magnetic field component of an MHD run, each with its
own threshold:

```python
ph.Simulation(
    # ...
    refinement="tagging",
    max_nbr_levels=4,
    tagging={"method": "lohner-per-direction", "quantities": {"rho": 0.4, "Bz": 0.3}},
)
```

Gradient scale on the magnetic field and the density of an MHD run, each with a reference magnitude
in its own units:

```python
ph.Simulation(
    # ...
    refinement="tagging",
    max_nbr_levels=3,
    tagging={
        "method": "gradient-scale",
        "quantities": {
            "B": {"threshold": 0.1, "ref_magnitude": 0.01},
            "rho": {"threshold": 0.2, "ref_magnitude": 0.05},
        },
    },
)
```

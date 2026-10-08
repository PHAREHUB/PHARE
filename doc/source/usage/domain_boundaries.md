# Domain boundaries

By default, every direction of the simulation domain is periodic. A direction becomes non-periodic
as soon as one of its two boundaries is given in the `Simulation` kwarg `domain_boundaries`, which then
sets the behavior of each domain boundary.

```{note}
Non-periodic domain boundaries are currently only supported by the MHD model
(`model_options=["MHDModel"]`).
```

## Specifying boundaries

`domain_boundaries` is a dict whose keys are boundary locations, written `"<direction><side>"` with
direction `x`, `y` or `z` and side `lower` or `upper`:

```text
domain_boundaries={
    "<location>": {"type": "<boundary type>", "<parameter>": <value>, ...},
    ...
}
```

- `type` is required and selects one of the boundary types listed below;
- the other keys are the parameters of that type: all of them are required, and keys that the type
  does not take are rejected.

The following rules are checked when the `Simulation` is created:

- a direction whose two locations are absent from `domain_boundaries` is periodic;
- a direction with one location given is non-periodic, and its other location must be given too;
- locations beyond the simulation dimension (e.g. `zlower` in 2D) are rejected.

The resulting periodicity of each direction is available as the list of booleans
`sim.periodicities`.

## Boundary types

| `type`                     | Use case                  | Parameters                             |
| -------------------------- | ------------------------- | -------------------------------------- |
| `open`                     | outflow                   | none                                   |
| `reflective`               | perfectly conducting wall | none                                   |
| `super-magnetofast-inflow` | super-magnetofast inflow  | `density`, `pressure`, `velocity`, `B` |

Each boundary type is implemented by imposing domain boundary conditions on the conservative
variables, on the magnetic field and, for some types, on the tangential component of the electric
field; see
{doc}`../numerics/domain_boundary_conditions` for how these conditions are discretized.

### Open

Density, momentum and total energy use a zero-gradient (Neumann) condition, and the magnetic field
in the ghost cells is extrapolated by the divergence-free transverse Neumann condition.

### Reflective

An impermeable, perfectly conducting wall:

- density and total energy use a zero-gradient (Neumann) condition;
- momentum is symmetric: zero normal component, zero-gradient tangential components;
- the magnetic field in the ghost cells is extrapolated by the divergence-free transverse Neumann
  condition.

Additionally, the tangential electric field is set to zero on the boundary. Hence the normal
component of $\vb{B}$ verifies

$$
\eval{\pdv{B_n}{t}}_{\Gamma} = - \eval{\left(\curl \vb{E}\right) \cdot \vb{n}}_{\Gamma} = \eval{\left(\pdv{E_{t_1}}{t_2} - \pdv{E_{t_2}}{t_1}\right)}_{\Gamma} = 0
$$

so that $B_n$ is constant in time on the boundary:

$$
\eval{B_n}_{\Gamma}(t) = \eval{B_n}_{\Gamma}(t=0).
$$ (eq:usage_domain_boundaries_reflective_normal_b)

Eq. {eq}`eq:usage_domain_boundaries_reflective_normal_b` is the relation that must be satisfied by
the magnetic field at the surface of a perfect conductor.

### Super-magnetofast inflow

Imposes a uniform inflow state. It is valid when the inflow speed normal to the boundary exceeds the
fast magnetosonic speed, so that all characteristics enter the domain. Parameters:

- `density`, `pressure`: finite positive scalars;
- `velocity`: a 3-vector, or a scalar giving the inward normal speed (its sign is set from the side
  of the boundary);
- `B`: a 3-vector.

```{note}
Only constant values are allowed to be passed as parameters for now.
```

In what follows, $\vb{v}_\text{in}$ and $\vb{B}_\text{in}$ denote the prescribed inflow velocity and
magnetic field.

Density, momentum and total energy are imposed through Dirichlet conditions. The tangential
components of the magnetic field are imposed in the ghost cells by the divergence-free transverse
Dirichlet condition, but the normal component cannot be enforced directly. Instead, the tangential
electric field is set on the boundary as:

$$
\eval{\vb{E}_t}_{\Gamma} = \vb{E}_{\text{in},t}, \qquad \vb{E}_\text{in} \equiv - \vb{v}_\text{in} \cross \vb{B}_\text{in}.
$$ (eq:usage_domain_boundaries_inflow_convective_electric_field)

Hence, according to Faraday's law

$$
\eval{\pdv{B_n}{t}}_{\Gamma} = - \eval{\left(\curl \vb{E}\right) \cdot \vb{n}}_{\Gamma} = \eval{\left( \pdv{E_{t_1}}{t_2} - \pdv{E_{t_2}}{t_1} \right)}_{\Gamma} = 0,
$$

so that $B_n$ remains constant throughout the simulation:

$$
\eval{B_n}_{\Gamma}(t) = \eval{B_n}_{\Gamma}(t=0).
$$ (eq:usage_domain_boundaries_inflow_normal_b)

```{note}
Since $B_n$ keeps its initial value on the boundary
(Eq. {eq}`eq:usage_domain_boundaries_inflow_normal_b`), the initial magnetic field must have the
same normal component as `B` on the whole inflow boundary. This is checked when the `MHDModel` is created, which raises a `ValueError` otherwise. The tangential
components of the initial magnetic field need not match `B`.
```

## Example

A 1D domain with a super-magnetofast inflow at `xlower` and an open outflow at `xupper`:

```python
sim = ph.Simulation(
    ...,
    model_options=["MHDModel"],
    domain_boundaries={
        "xlower": {
            "type": "super-magnetofast-inflow",
            "density": 1.0,
            "pressure": 1.0,
            "velocity": 2.0,
            "B": [0.75, 1.0, 0.0],
        },
        "xupper": {"type": "open"},
    },
)

ph.MHDModel(
    vx=lambda x: 2.0 + 0 * x,
    bx=lambda x: 0.75 + 0 * x,
    by=lambda x: 1.0 + 0 * x,
)
```

Here the fast magnetosonic speed of the inflow state is about 1.7 (with `gamma=5/3`), so the inflow
speed of 2.0 is super-magnetofast. The initial `bx` equals `B[0]` on the inflow boundary, as
required for the normal component of the magnetic field.

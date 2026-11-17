# Domain boundary conditions

## Domain boundary conditions for fields

PHARE imposes domain boundary conditions on fields through the values it sets in the ghost cells
beyond the boundary. Every stencil next to the domain edge reads values beyond it, so the boundary
condition is what decides those {term}`ghost cell` values. How to choose them in an input script is on {doc}`/user_guide/boundary_conditions`;
where each field component sits on the grid (primal or dual) is on {doc}`/numerics/grid`. In what follows:

- $\phi$ is a scalar field;
- $\vb*{\phi}$ is a vector field;
- $\Gamma$ is a portion of the domain boundary;
- $\vb{n}$ is the unit normal vector to $\Gamma$ directed inside the domain;
- $\vb{t}_1$ and $\vb{t}_2$ are unit tangent vectors to $\Gamma$ such that
  $(\vb{n}, \vb{t}_1, \vb{t}_2)$ is a right-handed orthonormal basis;
- $\phi_n = \vb*{\phi} \cdot \vb{n}$ and $\phi_{t_i} = \vb*{\phi} \cdot \vb{t}_i$ are the normal and
  tangential components of $\vb*{\phi}$;
- given a point $\vb{x}$, its mirror image across the boundary $\Gamma$ is denoted $\vb{x}^*$.

When a field is primal in the direction normal to $\Gamma$, some of its nodes lie exactly on
$\Gamma$; such nodes are their own mirror image.

```{eval-rst}
.. plot::
   :include-source: False
   :caption: Mirror images at a lower boundary :math:`\Gamma` (x = 0), with two ghost layers
      (grey). Top: a primal quantity; the node on :math:`\Gamma` is its own mirror, the ghost
      at :math:`-k\,\Delta x` mirrors the node at :math:`+k\,\Delta x`. Bottom: a dual
      quantity; the ghost at :math:`-(k-\frac12)\Delta x` mirrors the value at
      :math:`+(k-\frac12)\Delta x`.

   import matplotlib.pyplot as plt

   fig, axes = plt.subplots(2, 1, figsize=(5, 4), sharex=True)
   for ax, dual, title in zip(axes, (False, True), ("primal", "dual")):
       off = 0.5 if dual else 0.0
       ax.axvspan(-2.4, 0, color="0.88", lw=0)
       ax.axvline(0, color="k", lw=2)
       xs = [k + off for k in range(-2, 3)] if dual else list(range(-2, 3))
       for x in xs:
           ghost = x < 0
           ax.plot(x, 0, "o" if not dual else "s", ms=11,
                   color="0.6" if ghost else "#2a7ab9", zorder=3)
       ghosts = [x for x in xs if x < 0]
       for k, x in enumerate(sorted(ghosts, reverse=True)):
           h = 0.35 + 0.25 * k
           ax.annotate("", xy=(-x, 0.08), xytext=(x, 0.08),
                       arrowprops=dict(arrowstyle="<->", color="#e07a1f",
                                       connectionstyle=f"arc3,rad=-{h}"))
       ax.set_title(title, fontsize=14, loc="left")
       ax.set_ylim(-0.3, 1.15)
       ax.set_yticks([])
       ax.tick_params(labelsize=13)
       for side in ("top", "right", "left"):
           ax.spines[side].set_visible(False)
   axes[0].text(0.05, 0.75, "Γ", fontsize=15)
   axes[1].set_xlabel("x / Δx", fontsize=14)
   axes[1].set_xlim(-2.4, 2.8)
   fig.tight_layout()
```

### Dirichlet boundary condition

The Dirichlet boundary condition directly imposes the value of the field $\phi$ on a boundary:

$$
\eval{\phi}_{\Gamma } = \phi_0
$$ (eq:numerics_bc_dirichlet)

Two discretizations are available. The linear one enforces this condition at second order by
linear extrapolation at the ghost location $\vb{x}_G$, using the physical value at its mirror image
$\vb{x}_G^*$ across the boundary $\Gamma$:

$$
\frac{ \phi(\vb{x}_G) + \phi(\vb{x}_G^*)}{2} \simeq \phi_0   \implies \phi(\vb{x}_G) \simeq 2 \phi_0 - \phi(\vb{x}_G^*)
$$ (eq:numerics_bc_dirichlet_linear_extrapolation)

The constant one sets every ghost value to the imposed value, $\phi(\vb{x}_G) = \phi_0$. It is
first-order accurate on $\Gamma$, but it is exact for a uniform imposed state and does not depend on
the interior solution.

In both cases, a node lying on $\Gamma$ is directly set to $\phi_0$.

### Neumann boundary condition

The Neumann boundary condition imposes the value of the field's normal derivative with respect to
the boundary:

$$
\eval{\grad \phi \cdot \vb{n}}_{\Gamma } = q_0
$$ (eq:numerics_bc_neumann)

A second-order implementation can be obtained by linear extrapolation at the ghost location
$\vb{x}_G$, using the physical value at its mirror image $\vb{x}_G^*$ across the boundary $\Gamma$:

$$
\frac{ \phi(\vb{x}^*_G) - \phi(\vb{x}_G)}{\norm{\vb{x}^*_G - \vb{x}_G}}  \simeq q_0 \implies \phi(\vb{x}_G) \simeq \phi(\vb{x}_G^*) - q_0\norm{\vb{x}^*_G - \vb{x}_G}
$$ (eq:numerics_bc_neumann_linear_extrapolation)

```{note}
Only the zero-gradient case $q_0 = 0$ is currently implemented, for which the ghost value is a copy
of its mirror value: $\phi(\vb{x}_G) = \phi(\vb{x}_G^*)$.
```

### Symmetric boundary condition

For a scalar field, the symmetric boundary condition corresponds to a zero Neumann boundary
condition:

$$
\eval{\grad \phi \cdot \vb{n}}_{\Gamma } = 0.
$$ (eq:numerics_bc_symmetric_scalar)

For a vector field, it enforces a zero Dirichlet condition on the normal component and a zero
Neumann condition on the tangential components:

$$
\left\lbrace
\begin{aligned}
    & \eval{\phi_{n}}_{\Gamma } = 0 \\
    & \eval{\grad \phi_{t_1} \cdot \vb{n}}_{\Gamma } = 0 \\
    & \eval{\grad \phi_{t_2} \cdot \vb{n}}_{\Gamma } = 0 \\
\end{aligned}
\right.
$$ (eq:numerics_bc_symmetric_vector)

The Dirichlet condition on the normal component uses the linear discretization.

### Antisymmetric boundary condition

For a scalar field, the antisymmetric boundary condition corresponds to a zero Dirichlet boundary
condition:

$$
\eval{\phi}_{\Gamma } = 0.
$$ (eq:numerics_bc_antisymmetric_scalar)

For a vector field, it enforces a zero Neumann condition on the normal component and a zero
Dirichlet condition on the tangential components:

$$
\left\lbrace
\begin{aligned}
    & \eval{\grad \phi_{n} \cdot \vb{n}}_{\Gamma } = 0 \\
    & \eval{\phi_{t_1}}_{\Gamma } = 0 \\
    & \eval{\phi_{t_2}}_{\Gamma } = 0 \\
\end{aligned}
\right.
$$ (eq:numerics_bc_antisymmetric_vector)

The Dirichlet conditions on the tangential components use the linear discretization.

### Divergence-free conditions on the magnetic field

For the magnetic field, the boundary condition on the tangential components is complemented by a
condition on the normal component that is compatible with $\div \vb{B} = 0$. On $\Gamma$, the
solenoidal constraint reads

$$
\eval{\pdv{B_n}{n}}_{\Gamma} = - \eval{\nabla_t \cdot \vb{B}_t}_{\Gamma}
= - \eval{\left( \pdv{B_{t_1}}{t_1} + \pdv{B_{t_2}}{t_2} \right)}_{\Gamma},
$$ (eq:numerics_bc_divergence_free_normal_b)

which is a Neumann condition on $B_n$ whose data is set by the tangential components on $\Gamma$. It
is the only Neumann data compatible with a divergence-free field. The value of $B_n$ on $\Gamma$ itself
is not prescribed: it is updated by Faraday's law from the electric field on the boundary.

#### Divergence-free transverse Dirichlet

The tangential components are imposed, and the normal component follows from
Eq. {eq}`eq:numerics_bc_divergence_free_normal_b`:

$$
\left\lbrace
\begin{aligned}
    & \eval{\vb{B}_t}_{\Gamma} = \vb{B}_{0,t} \\
    & \eval{\pdv{B_n}{n}}_{\Gamma} = - \nabla_t \cdot \vb{B}_{0,t}
\end{aligned}
\right.
$$ (eq:numerics_bc_divergence_free_transverse_dirichlet)

The Neumann data of $B_n$ is here entirely determined by the prescribed tangential field.

#### Divergence-free transverse Neumann

The tangential components have a zero normal derivative, and the normal component follows from
Eq. {eq}`eq:numerics_bc_divergence_free_normal_b`:

$$
\left\lbrace
\begin{aligned}
    & \eval{\pdv{\vb{B}_t}{n}}_{\Gamma} = \vb{0} \\
    & \eval{\pdv{B_n}{n}}_{\Gamma} = - \eval{\nabla_t \cdot \vb{B}_t}_{\Gamma}
\end{aligned}
\right.
$$ (eq:numerics_bc_divergence_free_transverse_neumann)

The Neumann data of $B_n$ is here given by the tangential field of the solution on $\Gamma$.

#### Discretization

Both conditions are applied to the magnetic ghost cells at every ghost fill, that is at
initialization, at regrid and at each stage of the time integration. The tangential components are
filled first: by the constant Dirichlet discretization for the transverse Dirichlet condition, which
puts $\vb{B}_{0,t}$ in the ghost cells, and by the zero-gradient Neumann discretization for the
transverse Neumann condition. The normal component is then obtained by discretizing
Eq. {eq}`eq:numerics_bc_divergence_free_normal_b` cell by cell from the boundary outwards: each
ghost value of $B_n$ is computed from the one just set closer to the domain, so that the discrete
divergence of $\vb{B}$ vanishes in every ghost cell. The value of $B_n$ on $\Gamma$ is left unchanged.

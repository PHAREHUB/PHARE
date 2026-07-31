# The MHD models

## Dimensional equations

Here are the equations for deriving the MHD model:

$$
\begin{align}
\pdv{\rho}{t} + \div{\rho \vb{u}} &= 0 \\
\pdv{t}(\rho \vb{u}) + \div(\rho \vb{u} \vb{u}) &= -\grad P + \vb{j} \cross \vb{B} \\
\pdv{t}(\rho e_t) + \div(\rho \vb{u} e_t) &= - \div(P \vb{u}) + \vb{j} \cdot \vb{E} \\
\pdv{\vb{B}}{t} &= - \curl{\vb{E}} \\
\vb{E} &= -\vb{u} \cross \vb{B} + \frac{1}{ne} \vb{j} \cross \vb{B} + \eta \vb{j} - \nu \laplacian \vb{j} \\
\div \vb{B} &= 0 \\
\curl \vb{B} &= \mu_0 \vb{j}
\end{align}
$$

where $e$ is the internal energy, $\vb{u}$ the bulk velocity, $\vb{B}$ the
magnetic field, $\vb{E}$ the electric field, $\vb{j}$ the electric current
density, $P$ the thermal pressure, and $e_t = e + \flatfrac{\vb{u}^2}{2}$ the
total plasma energy. In Ohm's law, $\eta$ is the Joule resistivity and $\nu$
the hyper-resistivity.

With some manipulations, and by defining the total pressure $P* = P + \flatfrac{\vb{B}^2}{2 \mu_0}$,
the total energy per unit volume  $E = \rho e_t + \flatfrac{\vb{B}^2}{2 \mu_0}$,
and the Poynting vector $\vb*{\Pi} = \frac{\vb{E} \cross \vb{B}}{\mu_0}$ ,
we can obtain the following system of equations:

$$
\begin{align}
    & \pdv{\rho}{t} = -\div{\rho \vb{u}} \\
    & \pdv{t}(\rho \vb{u}) = - \div[  \rho \vb{u}\vb{u} + P^* \vb{I} - \frac{\vb{B} \vb{B}}{\mu_0}] \\
    & \pdv{E}{t} =  - \div [ ( \rho e_t + P) \vb{u} + \vb*{\Pi}] = - \div[ \left(E + P^* - \frac{\vb{B}^2}{\mu_0} \right) \vb{u} + \vb*{\Pi} ] \\
    & \pdv{B}{t} = - \curl \vb{E} \\
    & \vb{E} = -\vb{u} \cross \vb{B} + \frac{1}{ne} \vb{j} \cross \vb{B} + \eta \vb{j} - \nu \laplacian \vb{j} \\
    & \div \vb{B} = 0 \\
    & \curl \vb{B} = \mu_0 \vb{j}
\end{align}
$$

The system is closed with a polytropic equation of state:

$$
    \dv{}{t}\left(\frac{P}{\rho^\gamma}\right) = 0
$$ (eq:theory_mhd_dimensional_eos)

## Normalized equations

Using the reference quantities and normalization convention described in
{doc}`normalization`, the MHD equations become:

$$
\begin{aligned}
    & \pdv{\tilde{\rho}}{\tilde{t}} = -\tilde{\grad} \cdot \left( \tilde{\rho} \tilde{\vb{u}} \right) \\
    & \pdv{\tilde{t}}(\tilde{\rho} \tilde{\vb{u}}) = - \tilde{\grad} \cdot \left[  \tilde{\rho} \tilde{\vb{u}}\tilde{\vb{u}} + \tilde{P}^* \vb{I} - \tilde{\vb{B}} \tilde{\vb{B}} \right] \\
    & \pdv{\tilde{E}}{\tilde{t}} =  - \tilde{\grad} \cdot \left[  ( \tilde{\rho} \tilde{e}_t + \tilde{P}) \tilde{\vb{u}} + \tilde{\vb*{\Pi}}\right] = - \tilde{\grad} \cdot  \left[ \left(\tilde{E} + \tilde{P}^* - \tilde{\vb{B}}^2 \right) \tilde{\vb{u}} + \tilde{\vb*{\Pi}} \right] \\
    & \pdv{\tilde{\vb{B}}}{\tilde{t}} = - \tilde{\grad} \cross \vb{\tilde{E}} %= - \curl[ \eta \vb{j} - \tilde{\vb{u}} \cross \tilde{\vb{B}} ]
\end{aligned}
$$ 

with

$$
\begin{aligned}
    & \tilde{\vb*{\Pi}} = \tilde{\vb{E}} \cross \tilde{\vb{B}} \\
    & \tilde{\vb{E}} = \tilde{\eta} \tilde{\vb{j}} - \tilde{\vb{u}} \cross \tilde{\vb{B}} + \frac{1}{\tilde{n}} \tilde{\vb{j}} \cross \tilde{\vb{B}} - \tilde{\nu}  \tilde{\nabla}^2 \tilde{\vb{j}}
\end{aligned}
$$

In the following, equations will be considered under their non-dimensional form, with the tildes dropped for readability.

## Splitting of the magnetic field

To numerically impose an external magnetic field, the magnetic field is decomposed as

$$
\vb{B} = \vb{B}_0 + \vb{B}_1
$$

where $\vb{B}_0$ is the imposed external field and $\vb{B}_1$ is the difference between the total and external magnetic fields.

Inserting this splitting in the definition of the total energy $E$, one can define the "reduced" total energy $E_1$:

$$
E = \underbrace{\rho e_t + \frac{\vb{B}_1^2}{2}}_{E_1} + \frac{\vb{B}_0^2}{2} + \vb{B}_0 \cdot \vb{B}_1
$$

The derivative of $E$ is related to $E_1$'s one following

$$
\begin{aligned}
\pdv{E}{t} & = \pdv{E_1}{t} + \vb{B}_0 \cdot \pdv{\vb{B}_0}{t} + \vb{B}_0 \cdot \pdv{\vb{B}_1}{t} + \pdv{\vb{B}_0}{t} \cdot \vb{B}_1 \\
& = \pdv{E_1}{t} + \vb{B}_0 \cdot \pdv*{\vb{B}}{t} + \pdv{\vb{B}_0}{t} \cdot \vb{B}_1
\end{aligned}
$$

Moreover the term $\vb{B}_0 \cdot \pdv*{\vb{B}}{t}$ can be expanded using the induction equation as:

$$
\begin{aligned}
\vb{B}_0 \cdot \pdv*{\vb{B}}{t} &= - \vb{B}_0 \cdot \left( \curl \vb{E} \right) \\
&= - \div(\vb{E} \cross \vb{B}_0) - \vb{E} \cdot \underbrace{ \left( \curl \vb{B}_0\right)}_{ = \vb{j}_0  = \vb{0}}
\end{aligned} 
$$
The rotational of $\vb{B}_0$ is zero because it is external; the currents generating it are located outside of the computational domain.

Therefore, the equation for energy and magnetic field can be written in terms of $E_1$ and $B_1$ as follows

$$
\begin{aligned}
& \pdv{E_1}{t} =  - \grad \cdot \left[  ( \rho e_t + P) \vb{u} + \vb{E} \cross \vb{B}_1 \right] - \pdv{\vb{B}_0}{t} \cdot \vb{B}_1 % = - \grad \cdot  \left[ \left(E_1 + P_1^* - \vb{B}_1^2 \right) \vb{u} + \vb{E} \cross \vb{B}_1 \right]
\\
& \pdv{\vb{B}_1}{t} = - \grad \cross \vb{E} - \pdv{\vb{B}_0}{t}
\end{aligned}
$$

## Numerical scheme

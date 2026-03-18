# The MHD formalism

## Dimensional equations

Here are the following equations for deriving the MHD model.

- Maxwell-Thomson equation

$$
    \div \vb{B} = 0
$$ (eq:theory_mhd_dimensional_maxwell_thomson)

- Induction or Maxwell-Faraday equation

$$
    \pdv{B}{t} = - \curl{\vb{E}}
$$ (eq:theory_mhd_dimensional_induction_equation)

- Maxwell-Ampere without displacement currents

$$
    \curl \vb{B} = \mu_0 \vb{j}
$$ (eq:theory_mhd_dimensional_maxwell_ampere)

- Ohm's law with resistive, Hall, and hyper-resistive terms

$$
\vb{j} = \sigma \left( \vb{E} + \vb{u} \cross \vb{B} - \frac{1}{ne} \vb{j} \cross \vb{B} + \nu \laplacian \vb{j} \right)
$$ (eq:label:theory_mhd_dimensional_ohm_law)

- Mass conservation for the plasma

$$
\pdv{\rho}{t} + \div{\rho \vb{u}} = 0
$$ (eq:theory_mhd_dimensional_mass_conservation)

- Momentum conservation for the plasma assuming a scalar pressure

$$
\pdv{t}(\rho \vb{u}) + \div(\rho \vb{u} \vb{u}) = -\grad P + \vb{j} \cross \vb{B}
$$ (eq:theory_mhd_dimensional_momentum_conservation)

- Energy conservation for the total energy of the plasma $e_t =  e + \flatfrac{\vb{u}^2}{2} $ 

$$
\pdv{t}(\rho e_t)  + \div(\rho \vb{u} e_t) = - \div(P \vb{u}) + \vb{j} \cdot \vb{E}
$$ (eq:theory_mhd_dimensional_fluid_total_energy_conservation)

With some manipulations, and by defining the total pressure $P* = P + \flatfrac{\vb{B}^2}{2 \mu_0}$,
the total energy per unit volume  $E = \rho e_t + \flatfrac{\vb{B}^2}{2 \mu_0}$,
and the Poynting vector $\vb*{\Pi} = \frac{\vb{E} \cross \vb{B}}{\mu_0}$ ,
we can obtain the following system of equations:

$$
\begin{aligned}
    & \pdv{\rho}{t} = -\div{\rho \vb{u}} \\
    & \pdv{t}(\rho \vb{u}) = - \div[  \rho \vb{u}\vb{u} + P^* \vb{I} - \frac{\vb{B} \vb{B}}{\mu_0}] \\
    & \pdv{E}{t} =  - \div [ ( \rho e_t + P) \vb{u} + \vb*{\Pi}] = - \div[ \left(E + P^* - \frac{\vb{B}^2}{\mu_0} \right) \vb{u} + \vb*{\Pi} ] \\
    & \pdv{B}{t} = - \curl \vb{E} %= - \curl[ \eta \vb{j} - \vb{u} \cross \vb{B} ]
\end{aligned}
$$ (eq:theory_mhd_dimensional_conservative_system)

where the electric field is expressed thanks to Ohm's law as:

$$
\vb{E} = \eta \vb{j} - \vb{u} \cross \vb{B} + \frac{1}{n e} \vb{j} \cross \vb{B} - \nu \laplacian \vb{j}
$$

The system is closed with a polytropic equation of state:

$$
    P \rho^{-\gamma} = \mathrm{constant} \implies \rho e = \frac{P}{\gamma - 1}
$$ (eq:theory_mhd_dimensional_eos)

## Dimensionless equations

Let us introduce the following reference quantities:
- A reference magnetic field $B_0$
- A reference plasma density $n_0$
- The mass of the proton $m_p$
- The Alfven velocity $v_{0}$:

$$
v_{0} = \frac{B_0}{\sqrt{m_p n_0 \mu_0}}
$$
- The proton gyro-frequency $\omega_0$:

$$
\omega_0 = \frac{e B_0}{m_p}
$$

From these reference quantities, we can deduce the other reference quantities
- the proton inertial length $\delta_{0}$

$$
 \delta_0 = \frac{v_0}{\omega_0} = \frac{1}{e} \sqrt{\frac{m_p}{n_0 \mu_0}}
$$

- the reference mass density $\rho_0$:

$$
\rho_0 = m_p n_0
$$

- the reference pressure $P_0$

$$
P_0 = \rho_0 v_0^2 = \frac{B_0^2}{\mu_0}
$$

- the reference electric field $E_0$

$$
E_0 = V_0 B_0 = \frac{B_0^2}{\sqrt{m_p n_0 \mu_0}}
$$

- the reference charge current $j_0$:

$$
j_0 = \frac{B_0}{\delta_0 \mu_0} = e B_0 \sqrt{\frac{n_0 }{m_p \mu_0}}
$$

- the reference resistivity $\eta_0$

$$
 \eta_0 = \mu_0 v_0 \delta_0 = \frac{B_0 }{e n_0}
$$

- the reference hyper-resistivity $\nu_0$

$$
 \nu_0 = \mu v_0 \delta_0^3  = \frac{B_0 m_p}{e^3 n_0^2 \mu_0 }
$$

By defining non-dimensional fields as $\tilde{\phi} = \frac{\phi}{\phi_0}$, we can obtain the following equations:

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
\begin{align*}
\pdv{E}{t} & = \pdv{E_1}{t} + \vb{B}_0 \cdot \pdv{\vb{B}_0}{t} + \vb{B}_0 \cdot \pdv{\vb{B}_1}{t} + \pdv{\vb{B}_0}{t} \cdot \vb{B}_1 \\
& = \pdv{E_1}{t} + \vb{B}_0 \cdot \pdv*{\vb{B}}{t} + \pdv{\vb{B}_0}{t} \cdot \vb{B}_1
\end{align*}
$$

Moreover the term $\vb{B}_0 \cdot \pdv*{\vb{B}}{t}$ can be expanded using the induction equation as:

$$
\begin{align*}
\vb{B}_0 \cdot \pdv*{\vb{B}}{t} &= - \vb{B}_0 \cdot \left( \curl \vb{E} \right) \\
&= - \div(\vb{E} \cross \vb{B}_0) - \vb{E} \cdot \underbrace{ \left( \curl \vb{B}_0\right)}_{ = \vb{j}_0  = \vb{0}}
\end{align*} 
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

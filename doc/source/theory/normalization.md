# Normalization

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
E_0 = v_0 B_0 = \frac{B_0^2}{\sqrt{m_p n_0 \mu_0}}
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
 \nu_0 = \mu_0 v_0 \delta_0^3  = \frac{B_0 m_p}{e^3 n_0^2 \mu_0 }
$$

By defining non-dimensional fields as $\tilde{\phi} = \frac{\phi}{\phi_0}$,
any physical quantity $\phi$ is normalized to its reference value $\phi_0$.

## What this means in an input script

Every number you give PHARE, and every number it writes, is in these units:

| Quantity | Unit | In words |
|---|---|---|
| length (`dl`, domain size, `x`, `y`, `z`) | $\delta_0$ | proton inertial length |
| time (`time_step`, `final_time`) | $1/\omega_0$ | inverse proton gyro-frequency |
| magnetic field (`bx`, `by`, `bz`) | $B_0$ | |
| density (`density`) | $n_0$ | |
| velocity (`vbulk*`, `vth*`) | $v_0$ | Alfvén speed |
| mass, charge of a population | $m_p$, $e$ | proton mass and charge |
| temperature (`Te`) | $m_p v_0^2 = B_0^2/(\mu_0 n_0)$ | |
| resistivity, hyper-resistivity | $\eta_0$, $\nu_0$ | |

So a uniform plasma with `density=1` and `bx=1` has an Alfvén speed of 1, and a proton
gyro-period of $2\pi$.

Two useful consequences:

- In a Maxwellian population, `vth` is the standard deviation of the velocity along one
  direction, so the ion temperature is $T_i = m\,v_{th}^2$ (in units of $m_p v_0^2$).
- Pressure is in units of $B_0^2/\mu_0$, so the plasma beta is
  $\beta = 2nT/B^2$: for example $\beta_e = 2\,n\,T_e/B^2$.


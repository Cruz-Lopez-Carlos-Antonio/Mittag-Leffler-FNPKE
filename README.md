# Analytical solutions of the FNPKE, MATLAB codes
The present repository contains the MATLAB codes that were developed to solve the Fractional Neutron Point Kinetics Equations, using the Laplace transform and the Green Function and which were reported in the paper "Analytical Solution of the Fractional Neutron Point Kinetic Equations using the Mittag-Leffler function", submitted to the journal Computer Physics Communications.

The programs are licensed under a Creative Commons Attribution 4.0 International License: http://creativecommons.org/licenses/by/4.0/

Authors: Carlos-Antonio Cruz-López (cacl.nucl@gmail.com), Gilberto Espinosa-Paredes (gepe@xanum.uam.mx)

# Mathematical description of the problem
The codes solve the following fractional differential equation system:
$$\tau_\alpha D_C^{\alpha+1}\circ\ n\left(t\right)+\frac{dn\left(t\right)}{dt}+\tau_\alpha\left[\frac{P_{NL}\left(1-\rho\right)-1+\beta}{\Lambda}\right]D_C^\alpha\circ n(t)$$
$$=\frac{\rho-\beta}{\Lambda}n\left(t\right)+\tau_\alpha\lambda\ D_C^\alpha\circ\ C\left(t\right)+\lambda\ C\left(t\right),\tag{1}$$
$$\frac{dC\left(t\right)}{dt}=\frac{\beta}{\Lambda}n\left(t\right)-\lambda\ C(t) \tag{2}$$
where: $$n(t)$$ is the neutron density
```math
\left( \sum_{k=1}^n a_k b_k \right)^2 \leq \left( \sum_{k=1}^n a_k^2 \right) \left( \sum_{k=1}^n b_k^2 \right)
```


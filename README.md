# Analytical solutions of the FNPKE, MATLAB codes
The present repository contains the MATLAB codes that were developed to solve the Fractional Neutron Point Kinetics Equations, using the Laplace transform and the Green Function and which were reported in the paper "Analytical Solution of the Fractional Neutron Point Kinetic Equations using the Mittag-Leffler function", submitted to the journal Computer Physics Communications.

The programs are licensed under a Creative Commons Attribution 4.0 International License: http://creativecommons.org/licenses/by/4.0/

Authors: Carlos-Antonio Cruz-López (cacl.nucl@gmail.com), Gilberto Espinosa-Paredes (gepe@xanum.uam.mx)

# Mathematical description of the problem
The following fractional differential equation system is solved in the submitted paper:
$$\tau_\alpha D_C^{\alpha+1}\circ\ n\left(t\right)+\frac{dn\left(t\right)}{dt}+\tau_\alpha\left[\frac{P_{NL}\left(1-\rho\right)-1+\beta}{\Lambda}\right]D_C^\alpha\circ n(t)$$
$$=\frac{\rho-\beta}{\Lambda}n\left(t\right)+\tau_\alpha\lambda\ D_C^\alpha\circ\ C\left(t\right)+\lambda\ C\left(t\right),\tag{1}$$
$$\frac{dC\left(t\right)}{dt}=\frac{\beta}{\Lambda}n\left(t\right)-\lambda\ C(t) \tag{2}$$
where: $$\tau_{\alpha}, n(t),\rho, \alpha, P_{NL} \tag{3}$$ are the fractional relaxation time, the neutron density, the reactivity (assumed constant), the fractional order and the nonleakage probability, respectively. And:
 $$C(t), \beta, \Lambda, \lambda \tag{4}$$
are the concentration, the fraction and the decay of the precursors of the delayed neutrons. On the other hand, the following operations:
$$D_C^\alpha\circ\ n\left(t\right),\ \ D_C^\alpha\circ\ C\left(t\right)\ \tag{5} $$
represent the application of the fractional derivative operator, of the Caputo's type, defined as:
$$D_C^\alpha\circ\ f\left(t\right)=\frac{1}{\Gamma\left(n-1\right)}\int_{0}^{t}{f^{\left(n\right)}\left(\tau\right)\left(t-\tau\right)^{n-\alpha-1}d\tau \tag{6}}$$
# General description of the solution
In the submitted paper is showed that the analytical solution of the neutron density of the system given in Eq. (1) and (2) can be written as:
$$n\left(t\right)=\sum_{j=1}^{5}{b_jG_{4,Y_j,\delta}(t)} \tag{7}$$
and the analytical solution of the precursor of the delayed neutron as:

$$C(t)=\sum_{j=1}^{5}{h_jG_{6,\Omega_j,\delta}(t)}+C_0\exp(-\lambda\ t) \tag{8}$$
where the numbers:
$$b_j,\ h_j,\ 1\le\ j\le5 \tag{9}$$
are constants whose value depends on the initial conditions as well as on nuclear parameters. On the other hand, the sets:

$$Y_j=\set{y_{j,1},y_{j,2,} \ldots,\ y_{y,5}}, \Omega_j=\set{\omega_{j,1}, \omega_{j,2},\ldots,\omega_{j,5}}$$
$$\delta=\set{a_1,a_2,\ldots,a_5}\ \mathrm{and}\ \ \Upsilon=\set{c_1,c_2,\ldots,c_6} \tag{10}$$
are defined in terms of the fractional orders and the nuclear parameter. See the article for a full decription of such coefficients. Finally, the following function:
$$G_{n,\Psi,\delta}\left(t\right)=\frac{1}{a_n}\sum_{m=0}^{\infty}\frac{\left(-1\right)^m}{m!}\sum_{k_0+k_1+\ldots+k_{n-2}=m}\frac{m!}{k_0!k_1!\cdots k_{n-2}!}$$
$$\times\prod_{i=0}^{n-2}{\left(\frac{a_i}{a_n}\right)^{k_i}t^{\left(\phi_n-\phi_{n-1}\right)m+\phi_n+\sum_{j=0}^{n-2}{\left(\phi_{n-1}-\phi_j\right)k_j-1}}}$$
is known as the Green Function [1, p.158] [2,p.225]


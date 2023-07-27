# Analytical solutions of the FNPKE, MATLAB codes
The present repository contains the MATLAB codes that were developed to solve the Fractional Neutron Point Kinetics Equations, using the Laplace transform and the Green Function and which were reported in the paper "Analytical Solution of the Fractional Neutron Point Kinetic Equations using the Mittag-Leffler function", submitted to the journal of [Computer Physics Communications](https://www.sciencedirect.com/journal/computer-physics-communications).

The programs are licensed under a Creative Commons Attribution 4.0 International License: http://creativecommons.org/licenses/by/4.0/

Authors: Carlos-Antonio Cruz-López (cacl.nucl@gmail.com), Gilberto Espinosa-Paredes (gepe@xanum.uam.mx)

Mathematical and algorithmical generalities of the codes are described in the following lines with the purpose to provide some insight of the developed work. Nevertheless, a more detailed and precise discussion is provided in the submitted article.

## 1. Mathematical description of the problem
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
## 2. General description of the solution
As it is showed in the submitted paper, the analytical solution of the neutron density of the system given in Eq. (1) and (2) can be written as:
$$n\left(t\right)=\sum_{j=1}^{5}{b_jG_{4,Y_j,\delta}(t)} \tag{7}$$
and the analytical solution of the precursor of the delayed neutron as:

$$C(t)=\sum_{j=1}^{5}{h_jG_{6,\Omega_j,\delta}(t)}+C_0\exp(-\lambda\ t) \tag{8}$$
where the numbers:
$$b_j,\ h_j,\ 1\le\ j\le5 \tag{9}$$
are constants whose value depends on the initial conditions as well as on nuclear parameters, **see the paper for more details as well as the discussion of the theoretical convergence**. On the other hand, the sets:

$$Y_j=\set{y_{j,1},y_{j,2,} \ldots,\ y_{j,5}}, \Omega_j=\set{\omega_{j,1}, \omega_{j,2},\ldots,\omega_{j,5}}$$
$$\delta=\set{a_1,a_2,\ldots,a_5}\ \mathrm{and}\ \ \Upsilon=\set{c_1,c_2,\ldots,c_6} \tag{10}$$
are defined in terms of the fractional orders and the nuclear parameter. Finally, the following function:
$$G_{n,\Psi,\delta}\left(t\right)=\frac{1}{a_n}\sum_{m=0}^{\infty}\frac{\left(-1\right)^m}{m!}\sum_{k_0+k_1+\ldots+k_{n-2}=m}\frac{m!}{k_0!k_1!\cdots k_{n-2}!}$$
$$\times\prod_{i=0}^{n-2}{\left(\frac{a_i}{a_n}\right)^{k_i}t^{\left(\phi_n-\phi_{n-1}\right)m+\phi_n+S-1}}$$
$$\times\ E_{\phi_n-\phi_{n-1},\phi_n+S}^{\left(m\right)}\left(-\frac{a_{n-1}}{a_n}t^{\phi_n-\phi_{n-1}}\right) \tag{11}$$
is known as the Green Function [1, p.158] [2,p.225], which in turn includes the the m-derivative of the Mittag-Leffler function, defined as [3, p.66]:

$$E_{\alpha,\beta}^{\left(m\right)}\left(z\right)=\sum_{k=0}^{\infty}\frac{\left(m+k\right)!z^k}{k!\Gamma(\alpha\left(m+k\right)+\beta)} \tag{12},$$

where the following set of real numbers is included:

$$\Psi=\set{\phi_n,\phi_{n-1},\ldots\phi_1,\phi_0} \tag{13}$$

which, in turn allows defining the following sum:
$$S_j=\sum_{j=0}^{n-2}{\left(\phi_{n-2}-\phi_j\right)k_j} \tag{14}$$
## 3. Brief description of the developed algorithms
The algorithms that are described in the submitted article have the purpose to compute the equations (7) and (8) using the MATLAB language. Such task can be divided in three parts: 
1. Defining the constants of the sets given in terms of the nuclear reactor parameters.
2. Computing the sets related to the partitions of the integers.
3. Computing the 2-parameter Mittag-Leffler and its derivatives  

The three points are discussed in the submitted paper in a detailed way. Some general aspects of the second and third steps are described in the following lines.
### 3.1 Partitions of integers 
In order to optimize the way in which is computed the following sum of the Eq. (11):

$$\sum_{m=0}^{p}{\frac{\left(-1\right)^m}{m!}\sum_{k_0+k_1+\ldots+k_{n-2}=m}{f(k_0,k_2,\ldots,k_{n-2})}\ \tag{15} }$$
the set of different partitions:
$$K_{n,m}=\set{(k_1,k_2,\ldots,k_{n-2}|k_i\geq0,\ 1\le\ i\le\ n,\ \sum_{j=1}^{n-2}k_j=m}\ \tag{16} $$
is precomputed in the developed algorithm. Therefore, the original sum given in (15) can be valuated, in a fast way, as follows:





<details><summary>EXPAND SECTION 3.1. Partitions of integers</summary>
<p>


```ruby
   puts "Hello World"
```

</p>
</details>




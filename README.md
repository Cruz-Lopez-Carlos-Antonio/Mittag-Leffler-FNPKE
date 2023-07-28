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
$$\sum_{m=0}^{p}{\frac{\left(-1\right)^m}{m!}\sum_{k_0+k_1+\ldots+k_{n-2}=m}{f(k_0,k_2,\ldots,k_{n-2})}\ \tag{15}, }$$

(where the function f represents the different operations that are carried out on the right of such equation, and p the approximation of the infinite sum) the set of different partitions:

$$K_{n,m}=\set{(k_0,k_1,\ldots,k_{n-2}|k_i\geq0,\ 0\le\ i\le\ n,\ \sum_{j=0}^{n-2}k_j=m}\ \tag{16}$$

is precomputed in the developed algorithm. Therefore, the original sum given in (15) can be valuated, in a fast way, as follows:

$$\sum_{K_{n,0}}{f(K_{n,0})}-\sum_{K_{n,1}}\ f\left(K_{n,1}\right)+\ldots+\frac{\left(-1\right)^p}{p!}\sum_{K_{n,m}}{f(K_{n,m}). \tag{17}}$$

The mentioned sets are built using the following function, which is included in all the developed codes:

<details><summary>EXPAND PARTITIONS CODE. Partitions of integers</summary>
<p>


```MATLAB
%------------ Partitions of integers for the neutron density------------
%Section 5.1 of the paper
%Input: a natural number, n
%Output: all the partitios of the number n, considering 3 elements, i.e.
%x_1+x_2+x_3=n. The partitions are stored as arrays [x1, x_2, x_3]
function B = particiones(n)
L=[];

for k_0=0:n
    for k_1=0:n
        for k_2=0:n
            if deltakronecker(k_0+k_1+k_2,n)==1
                L = [L;k_0 k_1 k_2];
            end
        end
    end
end
B =L;
end
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------

```

</p>
</details>

A similar algorithm is used to compute the partitions that used for the precursors of the delayed neutrons. 
### 3.2 Computation of the 2-parameter Mittag-Leffler and its derivatives 
As it is discussed in the submitted paper, one of the most difficult parts when the analytical solutions is computed, is related to evaluate the 2-parameter Mittag-Leffler function. In fact, the definition provided in Eq. (12) is not convenient for large arguments, z, due to several digit precision issues (**see the article for more details**). 

Therefore it is necessary to use more advanced methods, as the one developed by R. Garrappa [4] who used an integral representation of such function in terms of its Laplace transform. R. Garrappa implement his powerful algorithm in MATLAB programming language, which is freely available in the [MathWorks official site](https://www.mathworks.com/matlabcentral/fileexchange/48154-the-mittag-leffler-function) .Garrappa's code can be used to compute the m-derivative (with m a non-negative integer), for which it is necessary to use the following relationship:

$$E_{\alpha,\beta}^{\left(m\right)}\left(z\right)=m!E_{\alpha,\alpha m+\beta}^{m+1}(z) \tag{18} $$

where the function that appears in the right side is known as the 3-parameter Mittag-Leffler function. 
#### 3.2.1 Invoking the Garrappa function. 
Garrappa's function is used in the developed algorithms in the following parts:

<details><summary>EXPAND MITTAG-LEFFLER CALCULATIONS. Mittag-Leffler calculations</summary>
<p>


```MATLAB

      %for m==0, we have the standard Mittag-Leffler function
      if m==0
          M1=ml(-1*(time^alpha_p)/a1,alpha_p,1+s_k012);
          M2=ml(-1*(time^alpha_p)/a1,alpha_p,1+alpha_p+s_k012);
          M3=ml(-1*(time^alpha_p)/a1,alpha_p,2+s_k012);
          M4=ml(-1*(time^alpha_p)/a1,alpha_p,3+s_k012);
          M5=ml(-1*(time^alpha_p)/a1,alpha_p,2+alpha_p+s_k012);
          M6=ml(-1*(time^alpha_p)/a1,alpha_p,3+alpha_p+s_k012);
      %for m>=1, we have derivatives of the Mittag-Leffler function
      %and the Eq. 60) is used.
      else
          M1=factorial(m)*ml(-1*(time^alpha_p)/a1,alpha_p,alpha_p*m+1+s_k012,m+1);
          M2=factorial(m)*ml(-1*(time^alpha_p)/a1,alpha_p,1+alpha_p+alpha_p*m+s_k012,m+1);
          M3=factorial(m)*ml(-1*(time^alpha_p)/a1,alpha_p,2+alpha_p*m+s_k012,m+1);
          M4=factorial(m)*ml(-1*(time^alpha_p)/a1,alpha_p,3+alpha_p*m+s_k012,m+1);
          M5=factorial(m)*ml(-1*(time^alpha_p)/a1,alpha_p,2+alpha_p*m+alpha_p+s_k012,m+1);
 
```
</p>
</details>

<details><summary>EXPAND MITTAG-LEFFLER CALCULATIONS. Mittag-Leffler calculations</summary>
<p>

```MATLAB
%for m==0, we have the standard Mittag-Leffler function
      if m==0
          M1=ml(-1*(time^alpha_p)/c1,alpha_p,2+s_k01234);
          M2=ml(-1*(time^alpha_p)/c1,alpha_p,2+alpha_p+s_k01234);
          M3=ml(-1*(time^alpha_p)/c1,alpha_p,3+s_k01234);
          M4=ml(-1*(time^alpha_p)/c1,alpha_p,4+s_k01234);
          M5=ml(-1*(time^alpha_p)/c1,alpha_p,3+alpha_p+s_k01234);
      %for m>=1, we have derivatives of the Mittag-Leffler function
      %and the Eq. 60) is used.    
      else
          M1=factorial(m)*ml(-1*(time^alpha_p)/c1,alpha_p,alpha_p*m+2+s_k01234,m+1);
          M2=factorial(m)*ml(-1*(time^alpha_p)/c1,alpha_p,alpha_p*m+2+alpha_p+s_k01234,m+1);
          M3=factorial(m)*ml(-1*(time^alpha_p)/c1,alpha_p,alpha_p*m+3+s_k01234,m+1);
          M4=factorial(m)*ml(-1*(time^alpha_p)/c1,alpha_p,alpha_p*m+4+s_k01234,m+1);
          M5=factorial(m)*ml(-1*(time^alpha_p)/c1,alpha_p,alpha_p*m+3+alpha_p+s_k01234,m+1);
      end
```
</p>
</details>

<details><summary>EXPAND MITTAG-LEFFLER CALCULATIONS. Mittag-Leffler calculations</summary>
<p>

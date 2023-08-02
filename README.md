# Analytical solutions of the FNPKE, MATLAB codes
The present repository contains the MATLAB codes that were developed to solve the Fractional Neutron Point Kinetics Equations, using the Laplace transform and the Green Function and which were reported in the paper "Analytical Solution of the Fractional Neutron Point Kinetic Equations using the Mittag-Leffler function", submitted to the journal of [Computer Physics Communications](https://www.sciencedirect.com/journal/computer-physics-communications).

The programs are licensed under a Creative Commons Attribution 4.0 International License: http://creativecommons.org/licenses/by/4.0/

Authors: Carlos-Antonio Cruz-López (cacl.nucl@gmail.com), Gilberto Espinosa-Paredes (gepe@xanum.uam.mx)

Mathematical and algorithmical generalities of the codes are described in the following lines with the purpose to provide some insight of the developed work. Nevertheless, a more detailed and precise discussion is provided in the submitted article.

## Financial Support.
The authors appreciate the financial support received from the Consejo Nacional de Ciencia y Tecnología, CONACYT, under the program “Estancias Posdoctorales por México, 2022”, with the project entitled: “Desarrollo de modelos fenomenológicos energéticos de orden fraccional, para la optimización y simulación en reactores nucleares de potencia”, by which the present development was possible.

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

<details><summary>CLICK HERE to expand the MATLAB code for partitions</summary>
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
#### 3.2.1 Invoking the Garrappa's function. 
Garrappa's function, written as ml(arg1,arg2,arg3,arg4), is used in the developed algorithms in the following parts:

<details><summary>CLICK HERE to expand the Implementation of the Mittag-Leffler in neutron density calculations.</summary>
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

<details><summary>CLICK HERE to expand the Mittag-Leffler implementation in the calculations of the precursors</summary>
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

where the arguments, denoted by arg1, arg2, arg3 and arg4 represent the parameters:
$$z,\ \alpha,\ \beta,\ \mathrm{and}\ m\ \tag{19},$$

respectively. As it can be observed in the CODES, it is possible to left the fourth argument blank for m=0. In order to use the Garrappa's code, it is necessary to download his MATLAB code and saving it in the same folder tha the developed codes. The following diagram resumes such procedure:

<details><summary>CLICK HERE to expand the diagram.</summary>
<p>

![image](https://github.com/Cruz-Lopez-Carlos-Antonio/Mittag-Leffler-FNPKE/assets/139827225/e909ffdd-d1ab-463f-a904-0beca1d44bc2)

</p>
</details>

## 4. Description of the developed codes and some examples of their application. 
### 4.1 FNPK-Insertion code
The FPNK-Insertion code, written in MATLAB language, solves the system given in Eq. (1) and Eq. (2) considering a constant reactivity. In the submitted paper it was concluded that such solution provided accurate results, always that the following inequality be fullfilled:

$$\left(\frac{t}{\tau}\right)^\alpha<80.42 \tag{20}$$

#### 4.1.1 First example of application. 
FNPK-Insertion code was used to reproduce data that was reported by Polo-Labarrios et al. [6]. For such task the following data was used:

<details><summary>CLICK HERE to expand the parameters that were used to the first comparison.</summary>
<p>
 
| Nuclear parameter | Value     | Nuclear parameter | Value           |
| ------------- | ------------- | -------------     | --------------  |
| tau           | 0.000125591s  | lambda            | 0.0769478s^(-1) |
| PNL           | 0.975         | beta              |0.00645          |
| rho           |0.002          |Lambda             |0.00005s         |

where the following notation was used:

$$\mathrm{tau}=\tau,\mathrm{lambda}=\lambda,\mathrm{PNL}=P_{NL},$$
$$\mathrm{beta}=\beta,\mathrm{rho}=\rho,\mathrm{Lambda}=\Lambda$$

</p>
</details>

The last data can be introduced in the code in the following lines:

<details><summary>CLICK HERE to expand the section of code where the data is introduced.</summary>
<p>
 
 ```MATLAB

 %-----------------Nuclear Data from the Model-------------------------
global tau lambda_p beta_p beta_p PNL rho LAMBDA_p;
tau = 3*9.21/220000
lambda_p =0.0769478
beta_p = 0.00645
PNL=0.975
rho = 0.002
LAMBDA_p=0.00005
%---------------------------------------------------------------------

```
</p>
</details>
Finally, the following initial conditions, related to a nuclear reactor in a steady state, were used:

$$n(t=0)=n(0)=n_0, \frac{dn\left(t\right)}{dt}|_{t=0}=0 \tag{21}$$

$$C(t=0)=C(0)=C_0=\frac{\beta}{\Lambda\lambda}n_0$$

which are introduced as follows:

<details><summary>CLICK HERE to expand the section of code where the initial conditions are introduced.</summary>
<p>
 
 ```MATLAB

%--------------------Initial conditions-------------------------------
global n_0 C_0 dn_0 q0;
n_0=1
C_0=n_0*beta_p/(LAMBDA_p*lambda_p)
%---------------------------------------------------------------------
```
</p>
</details>

**Example of an input and output #1**

For small times that fulfill the condition Eq. (20), it is possible to compute the neutron and the density of precursor specifiying the fractional order and the numbers of terms that will be used to approximate the sum. For such case the time step, denoted by **h**, will be equal to the final time. In order words, it will be not necessary to use sub-intervals. In the following lines is an example of calculations, which coincides with the value that is reported in Table 7 of the paper:

<details><summary>CLICK HERE to expand the first example of Input and output.</summary>
<p>

The wanted calculation: n(0.011) and C(0.0011)

INPUT:

```MATLAB

%---------------------Calling the function----------------------------
%first parameter = Target time
%second parameter = step
%third parameter = fractional order
%fourth parameter = number of the terms used to approx the infinite sum
Insertion(0.0011,0.0011,0.999999,10)
%----------------------------------------------------------------------

```

OUTPUT (console)

```
vect_sol =

   1.0e+03 *

                   0   0.001000000000000   1.676461185375020
   0.000001100000000   0.001039367880451   1.676464093008690
```
In this case we have three columns: the first one contains the time, the second one the neutron density and the final one the precursor's density. In the first row, we have the results for t=0, which are:

$$n\left(t=0 \right)=1,\ C\left(t=0 \right)=1676.461185375$$

In the second raw, in the other hand, we have:
$$n(0.0011s )=1.039367,\ \ C(0.0011s)=1676.464093$$

A similar output is provided in a .xls file. 
</p>
</details>

**Example of an input and output #2**
For large times or times that does not fulfill the Eq. (20), it is necessary to use a different approach that consists of dividing the time interval, and use an iterative procedure. For example, to find the neutron density at t=0.5 and alpha = 0.7, it is possible to use a step of **h=0.01** as follows:

<details><summary>CLICK HERE to expand the second example of Input and output.</summary>
<p>
The wanted calculation: n(0.5) and C(0.5)

INPUT:

```MATLAB


Insertion(0.5,0.01,0.7,10)


```

OUTPUT (console)

```
vect_sol =

   1.0e+03 *

                   0   0.001000000000000   1.676461185375020
   0.000010000000000   0.001293448128059   1.676772641112274
   0.000020000000000   0.001395244404846   1.677348012599482
   0.000030000000000   0.001430821255074   1.678014771682379
   0.000040000000000   0.001443518123640   1.678713306472076
   0.000050000000000   0.001448309667173   1.679423019937833
   0.000060000000000   0.001450369933096   1.680136795473135
   0.000070000000000   0.001451486659728   1.680852174231273
   0.000080000000000   0.001452277546963   1.681568306715986
   0.000090000000000   0.001452956022707   1.682284899494246
   0.000100000000000   0.001453595830699   1.683001851264176
   0.000110000000000   0.001454222450971   1.683719127109141
   0.000120000000000   0.001454844687462   1.684436715049993
   0.000130000000000   0.001455465582187   1.685154611032981
   0.000140000000000   0.001456086186297   1.685872813742782
   0.000150000000000   0.001456706863058   1.686591322810286
   0.000160000000000   0.001457327738062   1.687310138193358
   0.000170000000000   0.001457948854775   1.688029259962877
   0.000180000000000   0.001458570228292   1.688748688228805
   0.000190000000000   0.001459191863900   1.689468423114644
   0.000200000000000   0.001459813763500   1.690188464748613
   0.000210000000000   0.001460435927825   1.690908813260595
   0.000220000000000   0.001461058357199   1.691629468781087
   0.000230000000000   0.001461681051810   1.692350431440832
   0.000240000000000   0.001462304011796   1.693071701370696
   0.000250000000000   0.001462927237280   1.693793278701624
   0.000260000000000   0.001463550728377   1.694515163564624
   0.000270000000000   0.001464174485202   1.695237356090762
   0.000280000000000   0.001464798507868   1.695959856411163
   0.000290000000000   0.001465422796489   1.696682664657006
   0.000300000000000   0.001466047351178   1.697405780959529
   0.000310000000000   0.001466672172049   1.698129205450021
   0.000320000000000   0.001467297259215   1.698852938259832
   0.000330000000000   0.001467922612789   1.699576979520365
   0.000340000000000   0.001468548232886   1.700301329363081
   0.000350000000000   0.001469174119618   1.701025987919495
   0.000360000000000   0.001469800273099   1.701750955321178
   0.000370000000000   0.001470426693443   1.702476231699761
   0.000380000000000   0.001471053380764   1.703201817186925
   0.000390000000000   0.001471680335176   1.703927711914412
   0.000400000000000   0.001472307556792   1.704653916014018
   0.000410000000000   0.001472935045726   1.705380429617596
   0.000420000000000   0.001473562802092   1.706107252857054
   0.000430000000000   0.001474190826004   1.706834385864358
   0.000440000000000   0.001474819117576   1.707561828771528
   0.000450000000000   0.001475447676923   1.708289581710643
   0.000460000000000   0.001476076504158   1.709017644813836
   0.000470000000000   0.001476705599396   1.709746018213297
   0.000480000000000   0.001477334962750   1.710474702041273
   0.000490000000000   0.001477964594335   1.711203696430066
   0.000500000000000   0.001478594494266   1.711933001512036
```

The last value correspond to t=0.5, where it follows that:
$$n\left(0.5s\right)=1.47859\ $$
$$C\left(0.5s\right)=1711.93$$

which is very similar to the result given in Table 11, part 1 of the paper, where it was computed n(0.501) and C(0.501). In fact:

$$n\left(0.501s\right)=1.4785090\ $$
$$C\left(0.501s\right)=1705.2439$$


</p>
</details>

### 4.2 FNPK-ramp code
Analytical solutions can be used to simulate cases where reactivity is function of time. In such case it is necessary to divide the time interval in smalls sub intervals, assuming a constant value of the reactivity in each of them, given by:

$$\bar{\rho}=\frac{\ \rho\left(t_n\right)+\rho(t_{n-1})}{2} \tag{22}$$

where the limit and upper times are defined as:

$$t_{n}=\Delta t \cdot n = h \cdot n$$

$$t_{n-1}=\Delta t \cdot (n-1) = h \cdot (n-1) \tag{23}$$

Therefore, it follows that Eq. (22) is reduced to:

$$\bar{\rho}=\frac{\rho\left(nh\right)+\rho(\left(n-1\right)h)} $$

For the particular case of a ramp, it follows that:
$$\rho(t)=kt \tag{25}$$







 

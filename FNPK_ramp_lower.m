%Code used to simulate a ramp reactivity, Section 6.3.1
%Authors: Cruz-López, C., Espinosa-Paredes, G. 

%------------------------------------------------------------------
%Input: nuclear data, target time, initial conditions, step and ramp value
%Output: a .xls file with neutron and delayed precursors densities
%------------------------------------------------------------------

%--------------------------Required file----------------------------
%This code computes the Mittag-Leffler function (and its derivatives)
%using the MATLAB script that was developed by Roberto Garrappa in [1], 
%which can be freely download from the site of MATLAB [2]. Such function 
% is named as "ml" and returns the Mittag-Leffler function as well as its
%derivatives. It is only neceesary to save the two files (i.e. the present
%one and the developed by R. Garrappa) in the same Folder.
%---------------------------------------------------------------------

format long

%-----------------Nuclear Data from the Model-------------------------
global tau lambda_p beta_p beta_p PNL rho LAMBDA_p ramp;

tau = 3*9.21/220000
lambda_p =0.0787
beta_p = 0.00755
PNL=0.975
ramp = 0.0005
LAMBDA_p=0.003
%---------------------------------------------------------------------

%--------------------Initial conditions for the first time step--------
global n_0 C_0 dn_0 q0;
n_0=1
C_0=n_0*beta_p/(LAMBDA_p*lambda_p)
%---------------------------------------------------------------------

%---------------------Calling the function----------------------------
%first parameter = Target time
%second parameter = step
%third parameter = fractional order
%fourth parameter = number of the terms used to approx the infinite sum

ramp_reactivity(15,0.01,0.999999,10)
%---------------------------------------------------------------------


%---------Definition of the time discretization-------------------------
%Section 6.3.2 of the paper.
%The variable "final_time" is the target time,
%the "step" is the lenght of the time in which was
%divided the interval. "order" is the fractional order of the derivative
%and "approx" is the number of terms used in equations (45) and (50)

function D = ramp_reactivity(final_time,step,order, approx)
global vect_sol;
global tau lambda_p beta_p beta_p PNL LAMBDA_p ramp;
global n_0 C_0; 

vect_sol = [];
malla = ceil(final_time/step);
paso=step;
for i=0:malla
    i
    n_i = n_0;
    c_i=C_0;
   
    %Lower approax where the first point includes a time step
    
    rho = (ramp*i*paso+ramp*(i-1)*paso)/2;

    
    n_f = solution_neutrons(paso,n_0,C_0,order,approx,rho);
    c_f = solution_precursors(paso,n_0,C_0,order,approx,i,rho);
    n_0=n_f;
    C_0=c_f;
    
    vect_sol = [vect_sol; i*step n_f c_f ];
    size(vect_sol);
    i
end
vect_sol

filename = 'Densities_output_ramp_results_f_1.xlsx';
xlswrite(filename,vect_sol);

end


%-------------Analytical solution of the neutron density-------------
%Section 4.5, Eq.(45)of the paper.
%Input = time, initial conditions, order, approax
%Output = neutron density
%In this case the parameter alpha_p is related to the fractional order.

function C = solution_neutrons(time,N0,C0,alpha_p,approx,rho)
global tau lambda_p beta_p beta_p PNL LAMBDA_p;
a1=tau^alpha_p;
a2=1;
a3=a1*(lambda_p+(PNL*(1-rho)+beta_p-1)/LAMBDA_p);
a4=lambda_p-(rho-beta_p)/LAMBDA_p;
a5=(lambda_p*a1/LAMBDA_p)*((PNL*(1-rho)-1));
b1=N0*a1;
b2=N0;
b3= N0*a1*lambda_p+(N0*a1*(PNL*(1-rho)+beta_p-1))/LAMBDA_p;
b4=-lambda_p*C0*a1+(N0*a1*lambda_p*(PNL*(1-rho)+beta_p-1)/LAMBDA_p);
b5=lambda_p*N0+C0*lambda_p;
sol = 0;
for m=0:approx
    factor1 = (-1)^m;
    D = [ ]; 
    partial_sum=0;
    D = particiones(m);
    for u=1:size(D,1)
      s1 = factorial(m)/(factorial(D(u,1))*factorial(D(u,2))*factorial(D(u,3)));
      s2 = ((a5/a1)^(D(u,1)))*((a4/a1)^(D(u,2)))*((a3/a1)^(D(u,3)));
      s_k012 = (2-alpha_p)*D(u,1)+D(u,2)+(1-alpha_p)*D(u,3);
      s3 = (time^(alpha_p*m+s_k012))/factorial(m);
      if m==0
          M1=ml(-1*(time^alpha_p)/a1,alpha_p,1+s_k012);
          M2=ml(-1*(time^alpha_p)/a1,alpha_p,1+alpha_p+s_k012);
          M3=ml(-1*(time^alpha_p)/a1,alpha_p,2+s_k012);
          M4=ml(-1*(time^alpha_p)/a1,alpha_p,3+s_k012);
          M5=ml(-1*(time^alpha_p)/a1,alpha_p,2+alpha_p+s_k012);
          M6=ml(-1*(time^alpha_p)/a1,alpha_p,3+alpha_p+s_k012);
      else
          M1=factorial(m)*ml(-1*(time^alpha_p)/a1,alpha_p,alpha_p*m+1+s_k012,m+1);
          M2=factorial(m)*ml(-1*(time^alpha_p)/a1,alpha_p,1+alpha_p+alpha_p*m+s_k012,m+1);
          M3=factorial(m)*ml(-1*(time^alpha_p)/a1,alpha_p,2+alpha_p*m+s_k012,m+1);
          M4=factorial(m)*ml(-1*(time^alpha_p)/a1,alpha_p,3+alpha_p*m+s_k012,m+1);
          M5=factorial(m)*ml(-1*(time^alpha_p)/a1,alpha_p,2+alpha_p*m+alpha_p+s_k012,m+1);
          
      end
      SUMA_G=b1*M1+b2*(time^alpha_p)*M2+b3*time*M3+b4*(time^2)*M4+b5*(time^(1+alpha_p))*M5;
      partial_sum=partial_sum+s1*s2*s3*SUMA_G;
    end
    sol=sol+factor1*partial_sum;
end
solution=sol/a1;
C=solution;
end
%----------------------------------------------------------------------

%----------------Analytical solution of the precursors of delayed
%neutrons
function C = solution_precursors(time,N0,C0,alpha_p,approx,i,rho)
global tau lambda_p beta_p beta_p PNL LAMBDA_p;

a1=tau^alpha_p;
a2=1;
a3=a1*(lambda_p+(PNL*(1-rho)+beta_p-1)/LAMBDA_p);
a4=lambda_p-(rho-beta_p)/LAMBDA_p;
a5=(lambda_p*a1/LAMBDA_p)*((PNL*(1-rho)-1));
b1=N0*a1;
b2=N0;
b3= N0*a1*lambda_p-lambda_p*C0*a1+lambda_p*C0*a1+(N0*a1*(PNL*(1-rho)+beta_p-1))/LAMBDA_p;
b4=-lambda_p*C0*a1+(N0*a1*lambda_p*(PNL*(1-rho)+beta_p-1)/LAMBDA_p);
b5=lambda_p*N0+C0*lambda_p;
c1=a1;
c2=a2;
c3=a3+lambda_p*a1;
c4=a4+lambda_p*a2;
c5=a5+lambda_p*a3;
c6=lambda_p*a4;
c7=lambda_p*a5;
sol = 0;
i
for m=0:approx
    factor1 = (-1)^m;
    D = [ ]; 
    partial_sum=0;
    D = particiones_c(m);
    
    for u=1:size(D,1)
      s1 = 1/(factorial(D(u,1))*factorial(D(u,2))*factorial(D(u,3))*factorial(D(u,4)*factorial(D(u,5))));
      s2 = ((c7)^(D(u,1)))*((c6)^(D(u,2)))*((c5)^(D(u,3))*((c4)^(D(u,4)))*((c3)^(D(u,5))));
      s_k01234 = (3-alpha_p)*D(u,1)+2*D(u,2)+(2-alpha_p)*D(u,3)+D(u,4)+(1-alpha_p)*D(u,5);
      s3 = (time^(alpha_p*m+1+s_k01234));
      if m==0
          M1=ml(-1*(time^alpha_p)/c1,alpha_p,2+s_k01234);
          M2=ml(-1*(time^alpha_p)/c1,alpha_p,2+alpha_p+s_k01234);
          M3=ml(-1*(time^alpha_p)/c1,alpha_p,3+s_k01234);
          M4=ml(-1*(time^alpha_p)/c1,alpha_p,4+s_k01234);
          M5=ml(-1*(time^alpha_p)/c1,alpha_p,3+alpha_p+s_k01234);
          
      else
          M1=factorial(m)*ml(-1*(time^alpha_p)/c1,alpha_p,alpha_p*m+2+s_k01234,m+1);
          M2=factorial(m)*ml(-1*(time^alpha_p)/c1,alpha_p,alpha_p*m+2+alpha_p+s_k01234,m+1);
          M3=factorial(m)*ml(-1*(time^alpha_p)/c1,alpha_p,alpha_p*m+3+s_k01234,m+1);
          M4=factorial(m)*ml(-1*(time^alpha_p)/c1,alpha_p,alpha_p*m+4+s_k01234,m+1);
          M5=factorial(m)*ml(-1*(time^alpha_p)/c1,alpha_p,alpha_p*m+3+alpha_p+s_k01234,m+1);
      end
      SUMA_G=b1*M1+b2*(time^alpha_p)*M2+b3*time*M3+b4*(time^2)*M4+b5*(time^(1+alpha_p))*M5;
      partial_sum=partial_sum+s1*s2*s3*SUMA_G/c1^m;
    end
    sol=sol+factor1*partial_sum;
end
o1=sol/c1;
l=(beta_p/LAMBDA_p)*sol/c1;
l1=C0*exp(-1*lambda_p*time);
solution=(beta_p/LAMBDA_p)*sol/c1+C0*exp(-1*lambda_p*time);
C=solution;
end
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------



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

%---------------------Kronecker's delta defined in Eq. (56)-----------
function A = deltakronecker(n,m)
if n==m
    A=1;
else
    A=0;
end
end
%------------------------------------------------------------------------


%------ Partitions of integers for the Precursors of delayed neutrons-----
%Section 5.1 of the paper
%Input: a natural number, n
%Output: all the partitios of the number n, considering 5 elements, i.e.
%x_1+x_2+x_3+x_4+x_5=n, The partitions are stored as arrays 
%[x1, x_2, ...,x_5]
function B = particiones_c(n)
L=[];

for k_0=0:n
    for k_1=0:n
        for k_2=0:n
            for k_3=0:n
                for k_4=0:n
                    
                    if deltakronecker(k_0+k_1+k_2+k_3+k_4,n)==1
                        L = [L;k_0 k_1 k_2 k_3 k_4];
                    end
                end
            end
        end
    end
end
B =L;
end
%------------------------------------------------------------------------

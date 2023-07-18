%MATLAB code corresponding to the Eq. (45) of the Paper. 
% 
%Data from the Model
global tau lambda_p beta_p beta_p PNL rho LAMBDA_p;

tau = 3*9.21/220000
lambda_p =0.0787
beta_p = 0.00755
PNL=0.975
rho = 0.004
LAMBDA_p=0.003

format long
%Constants in terms of a's
global alpha_p a1 a2 a3 a4 a5 a6;
alpha_p=0.99999999
a1=tau^alpha_p
a2=1
a3=a1*(lambda_p+(PNL*(1-rho)+beta_p-1)/LAMBDA_p)
a4=lambda_p-(rho-beta_p)/LAMBDA_p
a5=(lambda_p*a1/LAMBDA_p)*((PNL*(1-rho)-1))

%Fractional and Model data
global n_0 C_0 
n_0=1
C_0=beta_p*n_0/(LAMBDA_p*lambda_p)

%Constants in terms of b
global b1 b2 b3 b4 b5 
b1=n_0*a1
b2=n_0
b3= n_0*a1*lambda_p+(n_0*a1*(PNL*(1-rho)+beta_p-1))/LAMBDA_p
b4=-lambda_p*C_0*a1+(n_0*a1*lambda_p*(PNL*(1-rho)+beta_p-1)/LAMBDA_p)
b5=lambda_p*n_0+C_0*lambda_p

solution(10,10)

function C = solution(time,approx)
global tau lambda_p beta_p beta_p PNL rho LAMBDA_p;
global alpha_p a1 a2 a3 a4 a5 a6; 
global n_0 C_0 dn_0 q;
global b1 b2 b3 b4 b5 b6;

sol = 0
for m=0:approx
    factor1 = (-1)^m;
    D = [ ]; 
    partial_sum=0;
    D = particiones(m);
    for u=1:size(D,1)
      s1 = 1/(factorial(D(u,1))*factorial(D(u,2))*factorial(D(u,3)));
      s2 = ((a5)^(D(u,1)))*((a4)^(D(u,2)))*((a3)^(D(u,3)));
      s_k012 = (2-alpha_p)*D(u,1)+D(u,2)+(1-alpha_p)*D(u,3);
      s3 = (time^(alpha_p*m+s_k012));
      if m==0
          M1=ml(-1*(time^alpha_p)/a1,alpha_p,1+s_k012);
          M2=ml(-1*(time^alpha_p)/a1,alpha_p,1+alpha_p+s_k012);
          M3=ml(-1*(time^alpha_p)/a1,alpha_p,2+s_k012);
          M4=ml(-1*(time^alpha_p)/a1,alpha_p,3+s_k012);
          M5=ml(-1*(time^alpha_p)/a1,alpha_p,2+alpha_p+s_k012);
          M6=ml(-1*(time^alpha_p)/a1,alpha_p,3+alpha_p+s_k012);
      else
          M1=ml(-1*(time^alpha_p)/a1,alpha_p,alpha_p*m+1+s_k012,m+1);
          M2=ml(-1*(time^alpha_p)/a1,alpha_p,1+alpha_p+alpha_p*m+s_k012,m+1);
          M3=ml(-1*(time^alpha_p)/a1,alpha_p,2+alpha_p*m+s_k012,m+1);
          M4=ml(-1*(time^alpha_p)/a1,alpha_p,3+alpha_p*m+s_k012,m+1);
          M5=ml(-1*(time^alpha_p)/a1,alpha_p,2+alpha_p*m+alpha_p+s_k012,m+1);
          
      end
      SUMA_G=b1*M1+b2*((time)^(alpha_p))*M2+b3*time*M3+b4*(time^2)*M4+b5*(time^(1+alpha_p))*M5;
      partial_sum=partial_sum+s1*s2*s3*factorial(m)*SUMA_G;
    end
    sol=sol+factor1*((1/a1)^(m))*partial_sum;
    m
    sol/a1
end
solution=sol/a1
end
           


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


function A = deltakronecker(n,m)
if n==m
    A=1;
else
    A=0;
end
end

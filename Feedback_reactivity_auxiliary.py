#Auxiliary code to solve the Point Kinetic Equations for a feedback reactivity
#considering a single group of precursors, mentioned in Section 6.3.3

import math

#------Nuclear parameters proposed by the authors----
LAMBDA_d = 5E-5
lambda_d=0.0787
Beta_d=0.00755


#----------Function that computes the roots-------
#INPUT: Lambda, lambda, Beta, Ramp
#OUTPUT: a list with two roots that do not depend
#on the initial conditions.
def roots(rho):
    a1 = (lambda_d*LAMBDA_d-rho+Beta_d)/(2*LAMBDA_d)
    b1 = lambda_d*rho/LAMBDA_d
    r1 = math.sqrt(a1**2+b1)-a1
    r2 = -math.sqrt(a1**2+b1)-a1
    return [r1,r2]


#------------Analytical solution---------------------
#See the Equations 29-31 of the GitHub's repository
#----------------------------------------------------
#Input parameters are the reactivity rho, the initial
#conditions of the neutron density and of the precursors
#respectively, and the time t.

def solution_t(rho,n0,c0,t):
    global LAMBDA_d,lambda_d,Beta_d
    sol = 0
    for i in range(2):
        num=n0*(R[i]+lambda_d)+lambda_d*c0
        a1 = (lambda_d*LAMBDA_d-rho+Beta_d)/LAMBDA_d
        den = 2*R[i]+a1
        sol = sol+(num/den)*math.exp(R[i]*t)
    return sol
#----------------------------------------------------
#----------------------------------------------------


#-----------Precursors of the delayed neutrons-----------
#See the Equations 29-31 of the GitHub's repository
#--------------------------------------------------------
#Input parameters are the same that the ones mentioned in
# neutron density
def solution_pre(rho,n0,c0,t):
    global LAMBDA_d,lambda_d,Beta_d
    sol = 0
    for i in range(2):
        num=(Beta_d/LAMBDA_d)*n0*(R[i]+lambda_d)+lambda_d*c0
        a1 = (lambda_d*LAMBDA_d-rho+Beta_d)/LAMBDA_d
        den = (2*R[i]+a1)*(R[i]+lambda_d)
        sol = sol+(num/den)*(math.exp(R[i]*t)-math.exp(-lambda_d*t))
    sol = sol+c0*math.exp(-lambda_d*t)
    return sol


#-----------------Solver---------------------------------
#--------------------------------------------------------
#-----------------------INPUT---------------------------- 
#Target time, step, parameters "a" and "b" from Table 17
#Nuclear parameters
#-----------------------OUTPUT---------------------------
#Three vectors that contains: 1) the time, 2)the neutron density
# and the precursors of delayed neutrons.
# Additionally the maximum is localted using the function
#max.

#Three vectos are used to save the results
Vector_solucion_n = [ ]
Vector_solucion_C = [ ]
Vector_tiempo = [ ]
#The step that is used for calculations
step = 0.001
#The final or disered time
Target = 2
#The initial condition for the neutron density
n0=1
#The initial condition for the precursors of the delayed neutron

c0=Beta_d/(lambda_d*LAMBDA_d)

for k in range(int(Target/step)+1):
    if k==0:
        rho = 0
    else:
        #Parameters given in Table 17 of the paper
        #a = 0.01, b=1E1-3
        rho = rho+step*(0.01-10**(-13)*n0)
    R = roots(rho)
    nt = solution_t(rho,n0,c0,step)
    Ct = solution_pre(rho,n0,c0,step)
    Vector_solucion_n.append(nt)
    Vector_solucion_C.append(Ct)
    Vector_tiempo.append(k*step)
    print(step*k,nt,Ct)
    n0=nt
    c0=Ct
    for w in R[:]:
        R.remove(w)

#Computation of the max values
print(max(Vector_solucion_n),Vector_tiempo[Vector_solucion_n.index(max(Vector_solucion_n))])
#The results are stored in the vectors

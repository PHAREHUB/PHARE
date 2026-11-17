#!/usr/bin/env python3
# Documentation check (not part of the test suite): rebuilds the SSPRK4_5 Shu-Osher recursion from the coefficients in
# src/amr/solvers/time_integrator/ssprk4_5_integrator.hpp:130-144, prints the amplification polynomial, the measured order
# on a nonlinear ODE, and the effective weights on L(U^n..U^4) (compared with the weights the code accumulates, :56-106).
import numpy as np
w0=0.391752226571890; w10=0.444370493651235; w11=0.555629506348765; w12=0.368410593050371
w20=0.620101851488403; w21=0.379898148511597; w22=0.251891774271694
w30=0.178079954393132; w31=0.821920045606868; w32=0.544974750228521
w40=0.517231671970585; w41=0.096059710526147; w42=0.063692468666290; w43=0.386708617503268; w44=0.226007483236906
print("row sums", w10+w11, w20+w21, w30+w31, w40+w41+w43)  # consistency (coefficients on U's)
def step(F,u,dt):
    U1=u+w0*dt*F(u)
    U2=w10*u+w11*U1+w12*dt*F(U1)
    U3=w20*u+w21*U2+w22*dt*F(U2)
    U4=w30*u+w31*U3+w32*dt*F(U3)
    return w40*U2+w41*U3+w42*dt*F(U3)+w43*U4+w44*dt*F(U4)
# Taylor coefficients for y'=lam*y
lam=1.0
# stability polynomial via series: apply to y'=z y with symbolic poly by sampling
import numpy.polynomial.polynomial as P
zs=np.linspace(-0.5,0.5,11)
vals=[step(lambda y:z*y,1.0,1.0) for z in zs]
c=np.polyfit(zs,vals,6)[::-1]
print("stability poly coeffs", c[:6], "exact", [1,1,.5,1/6,1/24,1/120])
# measured order on a nonlinear 2D ODE
def F(y): return np.array([y[1], -np.sin(y[0])])
def solve(n):
    y=np.array([1.0,0.0]); dt=2.0/n
    for _ in range(n): y=step(F,y,dt)
    return y
ref=solve(4096)
errs=[np.linalg.norm(solve(n)-ref) for n in (8,16,32,64)]
print(errs, [np.log2(errs[i]/errs[i+1]) for i in range(3)])
# effective weights on F(Un),F(U1),F(U2),F(U3),F(U4) from my recursion vs code's accumulateButcherFluxes_ expressions
import itertools
# represent state as (coeff on Un, [coeffs on dt*F_k])
def lin(a,s,b,t): return (a*s[0]+b*t[0], [a*x+b*y for x,y in zip(s[1],t[1])])
Un=(1.0,[0]*5)
def addF(s,k,c): f=list(s[1]); f[k]+=c; return (s[0],f)
U1=addF(Un,0,w0)
U2=addF(lin(w10,Un,w11,U1),1,w12)
U3=addF(lin(w20,Un,w21,U2),2,w22)
U4=addF(lin(w30,Un,w31,U3),3,w32)
fin=lin(w40,U2,w41,U3); fin=addF(fin,3,w42); fin=lin(1,fin,w43,U4); fin=addF(fin,4,w44)
mine=fin[1]
code=[w0*w11*w21*w31*w43+w0*w11*w21*w41+w0*w11*w40,
      w12*w21*w31*w43+w12*w21*w41+w12*w40,
      w22*w31*w43+w22*w41,
      w32*w43+w42,
      w44]
print("mine",mine); print("code",code); print("max diff",max(abs(a-b) for a,b in zip(mine,code)), "coef on Un", fin[0])

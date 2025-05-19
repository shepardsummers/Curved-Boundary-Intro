function f_eq_d2q9=eqm_d2q9(Rho,U)
%% f_eq_d2q9=eqm_d2q9(Rho,U) compute the D2Q9 equilibrium PDF based on the given physical density Rho and velocity U at any location.
%% Rho is the density and must be a scalar
%% U is the velocity and must be a column vector
%% f_eq_d2q9 is the PDF and must be a column vector

vLatt= [0     1     0    -1     0     1    -1    -1     1;...
              0     0     1     0    -1     1     1    -1    -1];

w=[4/9 1/9 1/9 1/9 1/9 1/36 1/36 1/36 1/36];

f_eq_d2q9=Rho*(1+3*vLatt'*U+9/2*((vLatt'*U).*(vLatt'*U))-3/2*(U')*U).*w';
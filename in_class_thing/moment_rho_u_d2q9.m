function [Rho,U]=moment_rho_u_d2q9(f)
%% function [Rho,U] = moment_rho_u_d2q9(f) compute the physical density and velocity based on the given PDF at any locaition using the D2Q9 lattice
%% f is the PDF and must be a column vector
%% Rho is the density and must be a scalar
%% U is the velocity and must be a column vector

vLatt= [0     1     0    -1     0     1    -1    -1     1;...
        0     0     1     0    -1     1     1    -1    -1];
Rho=sum(f);
U=(vLatt*f)/Rho;
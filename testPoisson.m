clc; clear; close all
Nx = 140; Ny = 80; Nz = 20;
x = linspace(0,1,Nx);
y = linspace(0,1,Ny);
z = linspace(0,1,Nz);
M = Nx*Ny*Nz;
hx = x(2) - x(1);
hy = y(2) - y(1);
hz = z(2) - z(1);

[X, Y, Z] = meshgrid(y,x,z);

epsilon = .1;
A = fem.createPoissonMatrix(Nx, Ny, Nz, hx, hy, hz);
RHS = ones(M,1)/epsilon;
wnodes = fem.bnodes('w', X, Y, Z);
enodes = fem.bnodes('e', X, Y, Z);
bnodes = union(wnodes,enodes);
inodes = setdiff(1:M,bnodes);
u = zeros(M,1); u(wnodes) = 0; u(enodes) = 0;

u = fem.solvePoisson(A, RHS, inodes, bnodes, u);

%show slice
figure
u1_slice = reshape(u,Nx,Ny,[]);
slice(X,Y,Z,u1_slice,[0.2 0.5 0.8],0.5,0.5)
shading interp
colorbar
title('u1')
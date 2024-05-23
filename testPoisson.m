clc; clear; close all
Nx = 20; Ny = 20; Nz = 20;
x = linspace(0,1,Nx);
y = linspace(0,1,Ny);
z = linspace(0,1,Nz);
M = Nx*Ny*Nz;
hx = x(2) - x(1);
hy = y(2) - y(1);
hz = z(2) - z(1);

[X, Y, Z] = meshgrid(y,x,z);

epsilon = 1;
A = fem.createPoissonMatrix(Nx, Ny, Nz, hx, hy, hz);
f = X*0;
f(4211) = 1;
RHS = f(:)/epsilon;
uex = sin(X.*Y) + exp(X) .*Z;
wnodes = fem.bnodes('w', X, Y, Z);
enodes = fem.bnodes('e', X, Y, Z);
nnodes = fem.bnodes('n', X, Y, Z);
snodes = fem.bnodes('s', X, Y, Z);
fnodes = fem.bnodes('f', X, Y, Z);
bcnodes = fem.bnodes('b', X, Y, Z);
bnodes = union(wnodes,enodes);
bnodes = union(bnodes,nnodes);
bnodes = union(bnodes,snodes);
bnodes = union(bnodes,fnodes);
bnodes = union(bnodes,bcnodes);
inodes = setdiff(1:M,bnodes);
u(bnodes) = 0;

u = fem.solvePoisson(A, RHS, inodes, bnodes, u(:));
%% 
err = norm(u(:)-uex(:),inf);

%show slice
figure
u1_slice = reshape(u,Nx,Ny,[]);
slice(X,Y,Z,u1_slice,[],0.5,[])
shading interp
colorbar
title('u')

%show slice
figure
slice(X,Y,Z,uex,[],0.5,[])
shading interp
colorbar
title('uex')


%% show slice
figure
slice(X,Y,Z,u1_slice-uex,[0.2 0.5 0.8],0.5,0.5)
shading interp
colorbar
title('uex')


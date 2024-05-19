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

q = [-1 1 -1 1];

posizioni = [0.16455696,0.54676259,0.47368421;   
             0.06329114,0.32374101,0.21052632; 
             0.22784810,0.75539568,0.21052632; 
             0.72151899,0.08633094,1];   

G = 0.7;
u = @(r) (G/sqrt(pi))^3*exp(-G^2*r.^2);

rho_lr = zeros(size(X));

for k = 1:length(q)
    r = sqrt((X - posizioni(k, 1)).^2 + (Y - posizioni(k, 2)).^2 + (Z - posizioni(k, 3)).^2);
    rho_lr = rho_lr + q(k)*u(r);
end

figure;
slice(X,Y,Z,rho_lr,[0.2 0.5 0.8],0.5,0.5)
shading interp
hold on;
scatter3(posizioni(:, 1), posizioni(:, 2), posizioni(:, 3), 100, q, 'filled');
text(posizioni(:, 1), posizioni(:, 2), posizioni(:, 3),num2str(q'), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
title('Densità di carica nello spazio');
xlabel('X');
ylabel('Y');
zlabel('Z');
grid on;
colorbar;
axis equal;
hold off;


epsilon = .1;
A = fem.createPoissonMatrix(Nx, Ny, Nz, hx, hy, hz);
RHS = reshape(rho_lr,M,1)/epsilon;
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








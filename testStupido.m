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

c = [-1 1]; %cariche
m = [10 10];
epsilon = .1;

%init
q0 = [0.16455696,0.54676259,0.47368421;   
             0.06329114,0.32374101,0.21052632];
             % 0.22784810,0.75539568,0.21052632; 
             % 0.72151899,0.08633094,1];

p0 = zeros(size(q0));


NP = size(q0,1);
ND = size(q0,2);

K = @(p)  sum(vecnorm(p').^2/2./m);

dKdp = @(p) p./m';

t = 0:0.1:50;

%% matrix
A = fem.createPoissonMatrix(Nx, Ny, Nz, hx, hy, hz);
wnodes = fem.bnodes('w', X, Y, Z);
enodes = fem.bnodes('e', X, Y, Z);
bnodes = union(wnodes,enodes);
inodes = setdiff(1:M,bnodes);
phi = zeros(M,1); phi(wnodes) = 0; phi(enodes) = 0;


%% Smooth long range 
G = 0.7;
u = @(r) (G/sqrt(pi))^3*exp(-G^2*r.^2);

F=@(q) force_long_range(q, X, Y, Z, Nx, Ny, Nz, hx, hy, hz, M, epsilon, u, c, A, inodes, bnodes, phi);


[q p] = int.velVerlet(q0,p0,F,dKdp,t);

figure
for i = 1:length(t)
    scatter3(q(:,1,i), q(:,2,i), q(:,3,i), 100, c, 'filled');
    text(q(:,1,i), q(:,2,i), q(:,3,i),num2str(c'), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
    xlim([0 1])
    ylim([0 1])
    zlim([0 1])
    drawnow

    clf
end




function F = force_long_range(q, X ,Y, Z, Nx, Ny, Nz, hx, hy, hz, M, epsilon, u, c, A, inodes, bnodes, phi)

rho_lr = zeros(size(X));

for k = 1:length(c)
    r = sqrt((X - q(k, 1)).^2 + (Y - q(k, 2)).^2 + (Z - q(k, 3)).^2);
    rho_lr = rho_lr + c(k)*u(r);
end

RHS = reshape(rho_lr,M,1)/epsilon;
phi = fem.solvePoisson(A, RHS, inodes, bnodes, phi);

phi_slice = reshape(phi,Nx,Ny,[]);

[dphi_x, dphi_y, dphi_z] = gradient(phi_slice, hx, hy, hz);

phix = c'.*interp3(X,Y,Z,dphi_x,q(:,1),q(:,2),q(:,3));
phiy = c'.*interp3(X,Y,Z,dphi_y,q(:,1),q(:,2),q(:,3));
phiz = c'.*interp3(X,Y,Z,dphi_z,q(:,1),q(:,2),q(:,3));

Fx = 1/2*phix;
Fy = 1/2*phiy;
Fz = 1/2*phiz;

F = [Fx,Fy,Fz];
end 


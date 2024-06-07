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

c = [-.5 .5]; %cariche
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
rCut = 3;
H = 3/pi/rCut^2;
u = @(r) H*(1-r/rCut).*(r<=rCut);

F=@(q) force_long_range(q,c, X ,Y, Z, Nx, Ny, Nz, hx, hy, hz, M, epsilon, u, A, inodes, bnodes, phi);


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

function F = force_long_range(q,c, X ,Y, Z, Nx, Ny, Nz, hx, hy, hz, M, epsilon, u, A, inodes, bnodes, phi)

    ptcls.x = q';
    ptcls.q = c;
    rho_lr = zeros(size(X));

    for k = 1:length(ptcls.q)
        r = sqrt((X - ptcls.x(1, k)).^2 + (Y - ptcls.x(2, k)).^2 + (Z - ptcls.x(3, k)).^2);
        rho_lr = rho_lr + ptcls.q(k)*u(r);
    end

    RHS = rho_lr(:)/epsilon;
    phi = fem.solvePoisson(A, RHS, inodes, bnodes, phi);

    phirec = reshape(phi,Nx,Ny,[]);

    [dphi_x, dphi_y, dphi_z] = gradient(phirec, hx, hy, hz);

    phix = ptcls.q.*interp3(X,Y,Z,dphi_x,ptcls.x(1,:),ptcls.x(2,:),ptcls.x(3,:),"nearest");
    phiy = ptcls.q.*interp3(X,Y,Z,dphi_y,ptcls.x(1,:),ptcls.x(2,:),ptcls.x(3,:),"nearest");
    phiz = ptcls.q.*interp3(X,Y,Z,dphi_z,ptcls.x(1,:),ptcls.x(2,:),ptcls.x(3,:),"nearest");

    Fx = - 1/2*phix;
    Fy = - 1/2*phiy;
    Fz = - 1/2*phiz;

    Fl = [Fx;Fy;Fz]';

    qiqj = -1;
    rc = 3;

    dx = ptcls.x(1,:)-ptcls.x(1,:)';
    dy = ptcls.x(2,:)-ptcls.x(2,:)';
    dz = ptcls.x(3,:)-ptcls.x(3,:)';

    r2 = dx.^2+dy.^2+dz.^2+eye(2);

    Fmat = -qiqj/epsilon.*(1./(r2.^(3/2))+ 4/rc^2 - 3*sqrt(r2)/rc^3);
    Fx = sum(Fmat.*dx,1);
    Fy = sum(Fmat.*dy,1);
    Fz = sum(Fmat.*dz,1);

    Fs = [Fx;Fy;Fz]';
    Fs = 0;

    F = Fs + Fl;
end 

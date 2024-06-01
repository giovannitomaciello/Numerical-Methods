clc; clear all; close all

%% generate particles
rng("default")
clc; clear
%% parameters of the simulation
L1 = 32; % L of the domain
L2 = 32; % H of the domain
L3 = 32; % D of the domain
epsilon = 10;
epsi = .1;
sigma = .5;
rCut = 4*sigma;
m = 10;
dt = 0.00025;

%% generate two rotating circles colliding
scale = 1;
delta = sigma*2^(1/6);

H1 = 10; W1 = 10; D1 = 10;

H1_l = (H1-1)*delta; W1_l = (W1-1)*delta; D1_l = (D1-1)*delta;

grd.ncy = L2/rCut;
grd.ncx = L1/rCut;
grd.ncz = L3/rCut;
grd.x = linspace (0, L1, grd.ncx+1);
grd.y = linspace (0, L2, grd.ncy+1);
grd.z = linspace (0, L3, grd.ncz+1);

[X, Y, Z] = meshgrid(linspace(-H1_l/2,H1_l/2,H1),linspace(-W1_l/2,W1_l/2,W1),linspace(-D1_l/2,D1_l/2,D1));

% cut particles outside the circle
X = X(:); Y = Y(:); Z = Z(:);
Xc = X(X.^2 + Y.^2 + Z.^2 <= (H1_l/2)^2);
Yc = Y(X.^2 + Y.^2 + Z.^2 <= (H1_l/2)^2);
Zc = Z(X.^2 + Y.^2 + Z.^2 <= (H1_l/2)^2);

Px = Yc*2; 
Py = -Xc*2;
Pz = zeros(size(Zc));

% copy to first circle
X1 = Xc + 16; Y1 = Yc-2.3 + 16; Z1 = Zc + 16;
%Px1 = Px/3; Py1 = Py/3+15; Pz1 = Pz + 0;
Px1 = Px*0; Py1 = 0*Px; Pz1 = 0*Px;

% copy to second circle
X2 = Xc + 16; Y2 = Yc+2.3 + 16; Z2 = Zc + 16;
%Px2 = Px/3; Py2 = Py/3-15; Pz2 = Pz + 0;
Px2 = 0*Px; Py2 = 0*Px2; Pz2 = 0*Px2;

ptcls.x = [[X1(:);X2(:)], [Y1(:);Y2(:)], [Z1(:);Z2(:)]]';
ptcls.p = [[Px1(:);Px2(:)], [Py1(:);Py2(:)], [Pz1(:);Pz2(:)]]';

carica1 = ones(1,numel(X1));
carica2 = -1*carica1;

ptcls.q = [carica1 carica2]*3;

grd_to_ptcl = sint.init_ptcl_mesh (grd, ptcls);
[nLy, nLx, nLz] = size(grd_to_ptcl);
counter = 0;
grd.removeIndex = [1:nLx*nLy, nLx*nLy*nLz-nLx*nLy+1:nLx*nLy*nLz];
for i = 1:nLz
    grd.removeIndex = unique([grd.removeIndex, 1+counter:nLy+counter, 1+counter:nLy:nLx*nLy+counter, ...
        nLy+counter:nLy:nLx*nLy+counter, nLy*nLx-nLy+1+counter:nLx*nLy+counter]);
        counter = counter + nLx*nLy;
end

d = cellfun (@numel, grd_to_ptcl, 'UniformOutput', true);

% number of time steps
tFinal = 2;
nTime = round(tFinal/dt);


Nx = grd.ncx*2; Ny = grd.ncy*2; Nz = grd.ncz*2;
x = linspace (0, L2, Nx);
y = linspace (0, L1, Ny);
z = linspace (0, L3, Nz);
M = Nx*Ny*Nz;
hx = x(2) - x(1);
hy = y(2) - y(1);
hz = z(2) - z(1);

[X, Y, Z] = meshgrid(y,x,z);


A = fem.createPoissonMatrix(Nx, Ny, Nz, hx, hy, hz);
A.Lx = L1; A.Ly = L2; A.Lz = L3;
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

remnodes = union(enodes,fnodes);
remnodes = union(remnodes,nnodes);

mat = A.A;
mat(wnodes,:) = mat(wnodes,:) + mat(enodes,:);
mat(snodes,:) = mat(snodes,:) + mat(nnodes,:);
mat(bcnodes,:) = mat(bcnodes,:) + mat(fnodes,:);

mat(:,wnodes) = mat(:,wnodes) + mat(:,enodes);
mat(:,snodes) = mat(:,snodes) + mat(:,nnodes);
mat(:,bcnodes) = mat(:,bcnodes) + mat(:,fnodes);

mat(:,remnodes) = [];
mat(remnodes,:) = [];
A.A = mat;
phi = 0*X;

H = 3/pi/rCut^2;
u = @(r) H*(1-r/rCut).*(r.^2<=rCut^2);

%% functions to pass to the integrator
force = @(dx, dy, dz, r2, ptcls, fc,indexPtclLocal,indexPtclAd) lennardJonesForce(dx, dy, dz, r2, ptcls, fc, indexPtclLocal,indexPtclAd, ...
    sigma, epsilon,epsi,rCut); % this is -Force (negative)
boundaryConditions = @(ptcls) updateBoundaryConditions(ptcls, L1, L2, L3, rCut, rCut, rCut);
ghost = @(ptcls, NP) updateGhost(ptcls, NP, L1, L2, L3, rCut, rCut, rCut);
dKdp = @(p) p/m;
forcelr=@(ptcls) force_long_range(ptcls, X ,Y, Z, Nx, Ny, Nz, hx, hy, hz, M, epsilon, u, A, remnodes, phi);


%% run the simulation
savingStep = 50;
[q, p] = sint.cellVelVerlet(force,forcelr, rCut^2, dKdp,dt,nTime,grd,ptcls,grd_to_ptcl,boundaryConditions,ghost,savingStep,1);

%% plot the results
figure
% save the video
v = VideoWriter('smallCollq3D.avi');
open(v)
for i = 1:size(q,3)
    scatter3(q(1,1:numel(X1),i), q(2,1:numel(X1),i), q(3,1:numel(X1),i), 10,"red" ,'filled')
    hold on
    scatter3(q(1,numel(X1)+1:end,i), q(2,numel(X1)+1:end,i), q(3,numel(X1)+1:end,i), 10,"blue" ,'filled')
    axis([rCut L1-rCut rCut L2-rCut rCut L3-rCut])
    view(0,90)
    drawnow
    frame = getframe(gcf);
    writeVideo(v,frame);
    axis equal
    hold off
end
close(v)

%% plot the results p
set(0,'DefaultFigureColor',[1 1 1]);
figure
% save the video
v = VideoWriter('smallCollp3D.avi');
open(v)
for i = 1:size(q,3)
    scatter3(q(1,:,i), q(2,:,i),  q(3,:,i), 10,vecnorm(p(:,:,i)) ,'filled')
    hold on
    axis([rCut L1-rCut rCut L2-rCut rCut L3-rCut])
    view(0,0)
    drawnow
    frame = getframe(gcf);
    writeVideo(v,frame);
    hold off
end
close(v)


%% FUNCTIONS

function [Fx, Fy, Fz] = lennardJonesForce(dx, dy, dz, r2, ptcls, fc, indexPtclLocal, indexPtclAd, sigmaij, epsij, epsilon, rc)
    qiqj = ptcls.q(indexPtclLocal)'.*ptcls.q(indexPtclAd);
    qiqj = qiqj(fc);

    % calculate sigma2 - 6
    sigma2 = (sigmaij^2)./r2;
    sigma6 = sigma2.*sigma2.*sigma2;

    % calculate force
    Fmat = 48*epsij*sigma6.*(sigma6 - 0.5)./r2 + qiqj/epsilon.*(1./(r2.^(3/2)) + 4/rc^3 - 3.*sqrt(r2)/rc^3);
    Fx = Fmat.*dx;
    Fy = Fmat.*dy;
    Fz = Fmat.*dz;
end

function F = force_long_range(ptcls, X ,Y, Z, Nx, Ny, Nz, hx, hy, hz, M, epsilon, u, A, remnodes, phi)

    % higly optimized with eigen3
    rho_lr = sint.ptclsToMeshInterp(X, Y, Z, ptcls.q, ptcls.x);

    RHS = reshape(rho_lr,[],1)/epsilon;
    phi = fem.solvePoissonPeriodicFFT(A, RHS);

    phirec = reshape(phi,Nx,Ny,[]);

    [dphi_x, dphi_y, dphi_z] = gradient(phirec, hx, hy, hz);

    phix = ptcls.q.*interp3(X,Y,Z,dphi_x,ptcls.x(1,:),ptcls.x(2,:),ptcls.x(3,:),"linear");
    phiy = ptcls.q.*interp3(X,Y,Z,dphi_y,ptcls.x(1,:),ptcls.x(2,:),ptcls.x(3,:),"linear");
    phiz = ptcls.q.*interp3(X,Y,Z,dphi_z,ptcls.x(1,:),ptcls.x(2,:),ptcls.x(3,:),"linear");

    F = [phix;phiy;phiz];
end 
function x = updateBoundaryConditions(x, L1, L2, L3, hx, hy, hz)
    if any(x(1,:) < hx) || any(x(1,:) > L1 - hx) || any(x(2,:) < hy) ...
        || any(x(2,:) > L2 - hy) || any(x(3,:) < hz) || any(x(3,:) > L3 - hz)
        % periodic boundary conditions using ghost cells
        % left boundary
        left = x(1,:) < hx;
        x(1,left) = x(1,left) + L1 - 2*hx;
        % right boundary
        right = x(1,:) > L1 - hx;
        x(1,right) = x(1,right) - L1 + 2*hx;
        % bottom boundary
        bottom = x(2,:) < hy;
        x(2,bottom) = x(2,bottom) + L2 - 2*hy;
        % top boundary
        top = x(2,:) > L2 - hy;
        x(2,top) = x(2,top) - L2 + 2*hy;
        % front boundary
        front = x(3,:) < hz;
        x(3,front) = x(3,front) + L3 - 2*hz;
        % back boundary
        back = x(3,:) > L3 - hz;
        x(3,back) = x(3,back) - L3 + 2*hz;
    end
end

function ptcls = updateGhost(ptcls, NP, L1, L2, L3, hx, hy, hz)
    ptcls.x = ptcls.x(:,1:NP);
    ptcls.q = ptcls.q(:,1:NP);
    if any(ptcls.x(1,:) < 2*hx) || any(ptcls.x(1,:) > L1-2*hx) ...
            || any(ptcls.x(2,:) < 2*hy) || any(ptcls.x(2,:) > L2-2*hx)...
            || any(ptcls.x(3,:) < 2*hz) || any(ptcls.x(3,:) > L3-2*hz)
        % periodic boundary conditions using ghost cells
        % left boundary
        left = ptcls.x(1,1:NP) < 2*hx;
        left = find(left);
        ptcls.x = [ptcls.x, [ptcls.x(1,left)+L1-2*hx; ptcls.x(2,left); ptcls.x(3,left)]];
        ptcls.q = [ptcls.q,ptcls.q(left)];

        % right boundary
        right = ptcls.x(1,1:NP) > L1-2*hx;
        right = find(right);
        ptcls.x = [ptcls.x, [ptcls.x(1,right)-L1+2*hx; ptcls.x(2,right); ptcls.x(3,right)]];
        ptcls.q = [ptcls.q,ptcls.q(right)];

        % bottom boundary
        bottom = ptcls.x(2,1:NP) < 2*hy;
        bottom = find(bottom);
        ptcls.x = [ptcls.x, [ptcls.x(1,bottom); ptcls.x(2,bottom)+L2-2*hy; ptcls.x(3,bottom)]];
        ptcls.q = [ptcls.q,ptcls.q(bottom)];

        % top boundary
        top = ptcls.x(2,1:NP) > L2-2*hy;
        top = find(top);
        ptcls.x = [ptcls.x, [ptcls.x(1,top); ptcls.x(2,top)-L2+2*hy; ptcls.x(3,top)]];
        ptcls.q = [ptcls.q,ptcls.q(top)];

        % front boundary
        front = ptcls.x(3,1:NP) < 2*hz;
        front = find(front);
        ptcls.x = [ptcls.x, [ptcls.x(1,front); ptcls.x(2,front); ptcls.x(3,front)+L3-2*hz]];
        ptcls.q = [ptcls.q,ptcls.q(front)];

        % back boundary
        back = ptcls.x(3,1:NP) > L3-2*hz;
        back = find(back);
        ptcls.x = [ptcls.x, [ptcls.x(1,back); ptcls.x(2,back); ptcls.x(3,back)-L3+2*hz]];
        ptcls.q = [ptcls.q,ptcls.q(back)];
    end
end


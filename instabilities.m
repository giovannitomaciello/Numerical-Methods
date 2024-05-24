clear all; clc; close all

scale = 1;

%% parameters of the simulation
rCut = 3;
L1 = 60 + 84*scale + 2*rCut;
L2 = 60 + 2*rCut;
L3 = 1 + 2*rCut;
epsilon = 1;
epsi = 1e-3;
sigma = 1;
m = 23;
noise = 0.1;

H1 = 19;
H2 = 38;
H3 = 19;
W = 180*L1/144;

grd.ncy = L2/rCut;
grd.ncx = L1/rCut;
grd.ncz = 3;
grd.x = linspace (0, L1 + 2*rCut, grd.ncx+1);
grd.y = linspace (0, L2 + 2*rCut, grd.ncy+1);
grd.z = linspace (0, L3 + 2*rCut, grd.ncz);

delta = 0.4;

[X1, Y1] = meshgrid (linspace(delta,L1-delta,W) + rCut, linspace(delta,L2/4 - delta, H1) + rCut);
[X2, Y2] = meshgrid (linspace(delta,L1-delta,W) + rCut, linspace(L2/4 + delta, 3*L2/4 - delta, H2) + rCut);
[X3, Y3] = meshgrid (linspace(delta,L1-delta,W) + rCut, linspace(3*L2/4 + delta ,L2 - delta, H3) + rCut);

N1 = numel(X1);
N2 = numel(X2);
N3 = numel(X3);
Nparticelle = N1 + N2 + N3;

ptcls.x = [[X1(:);X2(:);X3(:)], [Y1(:);Y2(:);Y3(:)],  4.7500+[Y1(:);Y2(:);Y3(:)]*0]';
ptcls.q(1:N1) = -0.5.*ones(1,N1);
ptcls.q((N1+1):(N1+N2)) = 0.5.*ones(1,N2);
ptcls.q((N1+N2+1):(N1+N2+N3)) = -0.5.*ones(1,N3);
ptcls.p = randn (size (ptcls.x)) * noise;

grd_to_ptcl = sint.init_ptcl_mesh (grd, ptcls);
d = cellfun (@numel, grd_to_ptcl, 'UniformOutput', true);

[nLy, nLx, nLz] = size(grd_to_ptcl);
counter = 0;
grd.removeIndex = [1:nLx*nLy, nLx*nLy*nLz-nLx*nLy+1:nLx*nLy*nLz];
for i = 1:nLz
    grd.removeIndex = unique([grd.removeIndex, 1+counter:nLy+counter, 1+counter:nLy:nLx*nLy+counter, ...
        nLy+counter:nLy:nLx*nLy+counter, nLy*nLx-nLy+1+counter:nLx*nLy+counter]);
        counter = counter + nLx*nLy;
end

%% functions to pass to the integrator
force = @(dx, dy, dz, r2, ptcls, fc,indexPtclLocal,indexPtclAd) lennardJonesForce(dx, dy, dz, r2, ptcls, fc, indexPtclLocal,indexPtclAd, ...
    sigma, epsilon, epsi, rCut); % this is -Force (negative)
boundaryConditions = @(ptcls) updateBoundaryConditions(ptcls, L1, L2, L3, rCut, rCut, rCut);
ghost = @(ptcls, NP) updateGhost(ptcls, NP, L1, L2, L3, rCut, rCut, rCut);
dKdp = @(p) p/m;

%% run the simulation
% number of time steps
dt = 0.001;
tFinal = 10;
nTime = round(tFinal/dt);
savingStep = 100;
[q, p] = sint.cellVelVerlet(force, rCut^2, dKdp,dt,nTime,grd,ptcls,grd_to_ptcl,boundaryConditions,ghost,savingStep);
%[q, p] = sint.cellForest(force, rCut^2, dKdp,dt,nTime,grd,ptcls,grd_to_ptcl,boundaryConditions,ghost,savingStep);

%% plot the results
figure
% save the video
v = VideoWriter('collq.avi');
open(v)
for i = 1:size(q,3)
    scatter(q(1,:,i), q(2,:,i), 5,ptcls.q ,'filled')
    axis([0 L1 0 L2])
    % axis equal
    % annotate the time
    text(0.8*L1,0.9*L2,sprintf("t = %.2f",i*tFinal/size(q,3)))
    drawnow
    frame = getframe(gcf);
    writeVideo(v,frame);
    hold off
end
close(v)

%% plot the results p
figure
% save the video
v = VideoWriter('collp.avi');
open(v)
for i = 1:size(q,3)
    scatter(q(1,:,i), q(2,:,i), 5,vecnorm(p(:,:,i)) ,"filled")
    hold on
    axis([0 L1 0 L2])
    axis equal
    colorbar
    colormap("turbo")
    clim([0 max(p,[],"all")])
    % annotate the time
    text(0.8*L1,0.9*L2,sprintf("t = %.2f",i*tFinal/size(q,3)))
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
    Fmat = 48*epsij*sigma6.*(sigma6 - 0.5)./r2 + qiqj/epsilon.*(1./(r2.^(3/2))+ 4/rc^2 - 3*sqrt(r2)/rc^3);
    Fx = Fmat.*dx;
    Fy = Fmat.*dy;
    Fz = Fmat.*dz;
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
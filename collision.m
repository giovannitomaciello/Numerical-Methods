clc; clear all; close all

%% parameters of the simulation
L1 = 250; % length of the domain
L2 = 200; % width of the domain
epsilon = 5;
sigma = 1;
rCut = 2.5*sigma;
m = 1;
dt = 0.00005;

%% generate particles
scale = 1;
delta = sigma*2^(1/6);

H1 = 40; H2 = 40;
W1 = 40; W2 = 160;

H1_l = (H1-1)*delta; W1_l = (W1-1)*delta;
H2_l = (H2-1)*delta; W2_l = (W2-1)*delta;

grd.ncy = L2/rCut;
grd.ncx = L1/rCut;
grd.x = linspace (0, L1, grd.ncx+1);
grd.y = linspace (0, L2, grd.ncy+1);

Dist = 5;
yc_1 = 120-20*(1-(1/scale));
yc_2 = 80-Dist+20*(1-(1/scale));

[X1, Y1] = meshgrid (linspace ((L1-W1_l)/2, (L1+W1_l)/2, W1),...
		     linspace ((yc_1)-H1_l/2, (yc_1)+H1_l/2, H1)); 
[X2, Y2] = meshgrid (linspace ((L1-W2_l)/2, (L1+W2_l)/2, W2),...
		     linspace ((yc_2)-H2_l/2, (yc_2)+H2_l/2, H2)); 

N1 = numel(X1);
N2 = numel(X2);
Nparticelle = N1 + N2;

ptcls.x = [[X1(:);X2(:)], [Y1(:);Y2(:)]]';
ptcls.p = randn (size (ptcls.x)) * .1;
ptcls.p(2, 1:N1) = ptcls.p(2, 1:N1) -10;
grd_to_ptcl = sint.init_ptcl_mesh (grd, ptcls);
[nLy, nLx] = size(grd_to_ptcl);
grd.removeIndex = unique([1:nLy, 1:nLy:nLx*nLy, ...
        nLy:nLy:nLx*nLy, nLy*nLx-nLy+1:nLx*nLy]);
d = cellfun (@numel, grd_to_ptcl, 'UniformOutput', true);

% number of time steps
tFinal = 20;
nTime = round(tFinal/dt);

%% functions to pass to the integrator
force = @(dx, dy, r2) lennardJonesForce(dx, dy, r2, sigma, epsilon); % this is -Force (negative)
boundaryConditions = @(ptcls) updateBoundaryConditions(ptcls, L1, L2, rCut, rCut);
ghost = @(ptcls, NP) updateGhost(ptcls, NP, L1, L2, rCut, rCut);
dKdp = @(p) p/m;

%% run the simulation
savingStep = 1000;
[q, p] = sint.cellVelVerlet(force, rCut^2, dKdp,dt,nTime,grd,ptcls,grd_to_ptcl,boundaryConditions,ghost,savingStep);
%[q, p] = sint.cellForest(force, rCut^2, dKdp,dt,nTime,grd,ptcls,grd_to_ptcl,boundaryConditions,ghost,savingStep);

%% plot the results
figure
% save the video
v = VideoWriter('collq.avi');
open(v)
for i = 1:size(q,3)
    scatter(q(1,1:numel(X1),i), q(2,1:numel(X1),i), 5,"red" ,'filled')
    hold on
    scatter(q(1,numel(X1)+1:end,i), q(2,numel(X1)+1:end,i), 5,"blue" ,'filled')
    axis([0 L1 0 L2])
    axis equal
    xline(L1-rCut)
    xline(rCut)
    yline(L2-rCut)
    yline(rCut)
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
    xline(L1-rCut)
    xline(rCut)
    yline(L2-rCut)
    yline(rCut)
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
function [Fx, Fy] = lennardJonesForce(dx, dy, r2, sigmaij, epsij)
    % calculate sigma2 - 6
    sigma2 = (sigmaij^2)./r2;
    sigma6 = sigma2.*sigma2.*sigma2;

    % calculate force
    Fmat = 48*epsij*sigma6.*(sigma6 - 0.5)./r2;
    Fx = Fmat.*dx;
    Fy = Fmat.*dy;
end


function x = updateBoundaryConditions(x, L1, L2, hx, hy)
    if any(x(1,:) < hx) || any(x(1,:) > L1 - hx) || any(x(2,:) < hy) || any(x(2,:) > L2 - hy)
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
    end
end

function ptcls = updateGhost(ptcls, NP, L1, L2, hx, hy)
    ptcls.x = ptcls.x(:,1:NP);
    if any(ptcls.x(1,:) < 2*hx) || any(ptcls.x(1,:) > L1-2*hx) ...
            || any(ptcls.x(2,:) < 2*hy) || any(ptcls.x(2,:) > L2-2*hx)
        % periodic boundary conditions using ghost cells
        % left boundary
        left = ptcls.x(1,1:NP) < 2*hx;
        left = find(left);
        ptcls.x = [ptcls.x, [ptcls.x(1,left)+L1-2*hx; ptcls.x(2,left)]];

        % right boundary
        right = ptcls.x(1,1:NP) > L1-2*hx;
        right = find(right);
        ptcls.x = [ptcls.x, [ptcls.x(1,right)-L1+2*hx; ptcls.x(2,right)]];

        % bottom boundary
        bottom = ptcls.x(2,1:NP) < 2*hy;
        bottom = find(bottom);
        ptcls.x = [ptcls.x, [ptcls.x(1,bottom); ptcls.x(2,bottom)+L2-2*hy]];

        % top boundary
        top = ptcls.x(2,1:NP) > L2-2*hy;
        top = find(top);
        ptcls.x = [ptcls.x, [ptcls.x(1,top); ptcls.x(2,top)-L2+2*hy]];
    end
end

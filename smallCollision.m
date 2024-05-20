clc; clear all; close all

%% generate particles
rng("default")
clc; clear
%% parameters of the simulation
L1 = 70; % L of the domain
L2 = 100; % H of the domain
epsilon = 5;
sigma = 1;
rCut = 5*sigma;
m = 1;
dt = 0.0005;
1;

%% generate two rotating circles colliding
scale = 1;
delta = sigma*2^(1/6);

H1 = 30; W1 = 30;

H1_l = (H1-1)*delta; W1_l = (W1-1)*delta;

grd.ncy = L2/rCut;
grd.ncx = L1/rCut;
grd.x = linspace (0, L1, grd.ncx+1);
grd.y = linspace (0, L2, grd.ncy+1);

[X, Y] = meshgrid(linspace(-H1_l/2,H1_l/2,H1),linspace(-W1_l/2,W1_l/2,W1));

% cut particles outside the circle
X = X(:); Y = Y(:);
Xc = X(X.^2 + Y.^2 <= (H1_l/2)^2);
Yc = Y(X.^2 + Y.^2 <= (H1_l/2)^2);

Px = Yc*2; 
Py = -Xc*2;

% copy to first circle
X1 = Xc + 32; Y1 = Yc-20 + 50;
Px1 = Px; Py1 = Py+10;

% copy to second circle
X2 = Xc + 38; Y2 = Yc+20 + 50;
Px2 = Px; Py2 = Py-10;

ptcls.x = [[X1(:);X2(:)], [Y1(:);Y2(:)]]';
ptcls.p = [[Px1(:);Px2(:)], [Py1(:);Py2(:)]]';
grd_to_ptcl = sint.init_ptcl_mesh (grd, ptcls);
[nLy, nLx] = size(grd_to_ptcl);
grd.removeIndex = unique([1:nLy, 1:nLy:nLx*nLy, ...
        nLy:nLy:nLx*nLy, nLy*nLx-nLy+1:nLx*nLy]);
d = cellfun (@numel, grd_to_ptcl, 'UniformOutput', true);

% number of time steps
tFinal = 4;
nTime = round(tFinal/dt);

%% functions to pass to the integrator
force = @(dx, dy, r2) lennardJonesForce(dx, dy, r2, sigma, epsilon); % this is -Force (negative)
boundaryConditions = @(ptcls) updateBoundaryConditions(ptcls, L1, L2, rCut, rCut);
ghost = @(ptcls, NP) updateGhost(ptcls, NP, L1, L2, rCut, rCut);
dKdp = @(p) p/m;

%% run the simulation
savingStep = 10;
[q, p] = sint.cellVelVerlet(force, rCut^2, dKdp,dt,nTime,grd,ptcls,grd_to_ptcl,boundaryConditions,ghost,savingStep);
%[q, p] = sint.cellForest(force, rCut^2, dKdp,dt,nTime,grd,ptcls,grd_to_ptcl,boundaryConditions,ghost,savingStep);

%% plot the results
figure
% save the video
v = VideoWriter('smallCollq.avi');
open(v)
for i = 1:5:size(q,3)
    scatter(q(1,1:numel(X1),i), q(2,1:numel(X1),i), 10,"red" ,'filled')
    hold on
    scatter(q(1,numel(X1)+1:end,i), q(2,numel(X1)+1:end,i), 10,"blue" ,'filled')
    axis([0 L1 0 L2])
    xline(L1-rCut)
    xline(rCut)
    yline(L2-rCut)
    yline(rCut)
    drawnow
    frame = getframe(gcf);
    writeVideo(v,frame);
    hold off
end
close(v)

%% plot the results p
figure
% save the video
v = VideoWriter('smallCollp.avi');
open(v)
for i = 1:5:size(q,3)
    scatter(q(1,:,i), q(2,:,i), 10,vecnorm(p(:,:,i)) ,'filled')
    hold on
    axis([0 L1 0 L2])
    xline(L1-rCut)
    xline(rCut)
    yline(L2-rCut)
    yline(rCut)
    drawnow
    frame = getframe(gcf);
    writeVideo(v,frame);
    hold off
end
close(v)

%% calculate Energy
kb = 1.380658e-23;
for i = 1:size(q,3)
    KU(i,:) = Energy(q(:,:,i), p(:,:,i), m, sigma, epsilon)*(1.66e-27);
end

%%
figure
plot(linspace(1,tFinal,size(q,3)),KU(:,1),"LineWidth",1.4)
hold on
plot(linspace(1,tFinal,size(q,3)),KU(:,2),"LineWidth",1.4)
plot(linspace(1,tFinal,size(q,3)),KU(:,1) + KU(:,2),"LineWidth",1.4)
xlabel("Time [-]","FontSize",11,"FontWeight","bold")
ylabel("Energy","FontSize",11,"FontWeight","bold")
legend(["K" "U" "Etot"])

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

function E = Energy(q, p, m, sigmaij, epsij)
    n = size(q,2);

    % calculate kinetic energy
    K = sum(m.*( sum( (p./m).^2 ,2) )) /2;

    % calculate distance
    dx = q(1,:) - q(1,:)';
    dy = q(2,:) - q(2,:)';
    r2 = dx.^2 + dy.^2  + eye(n);

    % calculate sigma2 - 6
    sigma2 = (sigmaij^2)./r2;
    sigma6 = sigma2.*sigma2.*sigma2;

    %set diagonal to 0
    sigma6(1:n+1:end) = 0;

    % calculate potential energy
    U = 4*epsij*sigma6.*(sigma6 - 1);

    % sum potential energy
    U = sum(sum(U))/2;

    E = [K U];
end
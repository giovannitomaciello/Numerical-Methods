clc; clear all; close all

%% generate particles
rng("default")
clc; clear
%% parameters of the simulation
L1 = 40; % L of the domain
L2 = 40; % H of the domain
L3 = 40; % D of the domain
epsilon = 20;
sigma = .5;
rCut = 4*sigma;
m = 1;
dt = 0.0005;

%% generate two rotating circles colliding
scale = 1;
delta = sigma*2^(1/6);

H1 = 15; W1 = 15; D1 = 15;

H1_l = (H1-1)*delta; W1_l = (W1-1)*delta; D1_l = (D1-1)*delta;

grd.ncy = L2/rCut;
grd.ncx = L1/rCut;
grd.ncz = L3/rCut;
grd.x = linspace (0, L1, grd.ncx+1);
grd.y = linspace (0, L2, grd.ncy+1);
grd.z = linspace (0, L3, grd.ncz+1);
nLx = grd.ncx;
nLy = grd.ncy;
nLz = grd.ncz;
counter = 0;
grd.removeIndex = [1:nLx*nLy, nLx*nLy*nLz-nLx*nLy+1:nLx*nLy*nLz];
for i = 1:nLz
    grd.removeIndex = unique([grd.removeIndex, 1+counter:nLy+counter, 1+counter:nLy:nLx*nLy+counter, ...
        nLy+counter:nLy:nLx*nLy+counter, nLy*nLx-nLy+1+counter:nLx*nLy+counter]);
        counter = counter + nLx*nLy;
end

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
X1 = Xc + 18; Y1 = Yc-5 + 20; Z1 = Zc + 15;
Px1 = Px/3; Py1 = Py/3+15; Pz1 = Pz + 0;

% copy to second circle
X2 = Xc + 18; Y2 = Yc+5 + 20; Z2 = Zc + 15;
Px2 = Px/3; Py2 = Py/3-15; Pz2 = Pz + 0;

ptcls.x = [[X1(:);X2(:)], [Y1(:);Y2(:)], [Z1(:);Z2(:)]]';
ptcls.p = [[Px1(:);Px2(:)], [Py1(:);Py2(:)], [Pz1(:);Pz2(:)]]';

carica1 = ones(1,numel(X1));
carica2 = -1*carica1;

ptcls.q = [[carica1 carica2];[carica1 carica2];[carica1 carica2]];
epsi = .1;
 

grd_to_ptcl = sint.init_ptcl_mesh (grd, ptcls);
d = cellfun (@numel, grd_to_ptcl, 'UniformOutput', true);

% number of time steps
tFinal = 4;
nTime = round(tFinal/dt);

%% functions to pass to the integrator
force = @(dx, dy, dz, r2, ptcls, fc,indexPtclLocal,indexPtclAd) lennardJonesForce(dx, dy, dz, r2, ptcls, fc, indexPtclLocal,indexPtclAd, ...
    sigma, epsilon,epsi,rCut); % this is -Force (negative)
boundaryConditions = @(ptcls) updateBoundaryConditions(ptcls, L1, L2, L3, rCut, rCut, rCut);
ghost = @(ptcls, NP) updateGhost(ptcls, NP, L1, L2, L3, rCut, rCut, rCut);
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
for i = 1:size(q,3)
    scatter3(q(1,1:numel(X1),i), q(2,1:numel(X1),i), q(3,1:numel(X1),i), 10,"red" ,'filled')
    hold on
    scatter3(q(1,numel(X1)+1:end,i), q(2,numel(X1)+1:end,i), q(3,numel(X1)+1:end,i), 10,"blue" ,'filled')
    axis([0 L1 0 L2 0 L3])
    drawnow
    frame = getframe(gcf);
    writeVideo(v,frame);
    axis equal
    hold off
end
close(v)

%% plot the results p
figure
% save the video
v = VideoWriter('smallCollp.avi');
open(v)
for i = 1:size(q,3)
    scatter3(q(1,:,i), q(2,:,i),  q(2,:,i), 10,vecnorm(p(:,:,i)) ,'filled')
    hold on
    axis([0 L1 0 L2 0 L3])
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

function [Fx, Fy, Fz] = lennardJonesForce(dx, dy, dz, r2, ptcls, fc, indexPtclLocal,indexPtclAd, sigmaij, epsij, epsilon, rc)
    qiqj = ptcls.q(indexPtclLocal)*ptcls.q(indexPtclAd)';
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

function x = updateGhost(x, NP, L1, L2, L3, hx, hy, hz)
    if any(x(1,:) < 2*hx) || any(x(1,:) > L1-2*hx) ...
            || any(x(2,:) < 2*hy) || any(x(2,:) > L2-2*hx)...
            || any(x(3,:) < 2*hz) || any(x(3,:) > L3-2*hz)
        % periodic boundary conditions using ghost cells
        % left boundary
        left = x(1,1:NP) < 2*hx;
        left = find(left);
        x = [x, [x(1,left)+L1-2*hx; x(2,left); x(3,left)]];

        % right boundary
        right = x(1,1:NP) > L1-2*hx;
        right = find(right);
        x = [x, [x(1,right)-L1+2*hx; x(2,right); x(3,right)]];

        % bottom boundary
        bottom = x(2,1:NP) < 2*hy;
        bottom = find(bottom);
        x = [x, [x(1,bottom); x(2,bottom)+L2-2*hy; x(3,bottom)]];

        % top boundary
        top = x(2,1:NP) > L2-2*hy;
        top = find(top);
        x = [x, [x(1,top); x(2,top)-L2+2*hy; x(3,top)]];

        % front boundary
        front = x(3,1:NP) < 2*hz;
        front = find(front);
        x = [x, [x(1,front); x(2,front); x(3,front)+L3-2*hz]];

        % back boundary
        back = x(3,1:NP) > L3-2*hz;
        back = find(back);
        x = [x, [x(1,back); x(2,back); x(3,back)-L3+2*hz]];
    end
end

function E = Energy(q, p, m, sigmaij, epsij)
    n = size(q,2);

    % calculate kinetic energy
    K = sum(m.*( sum( (p./m).^2 ,2) )) /2;

    % calculate distance
    dx = q(1,:) - q(1,:)';
    dy = q(2,:) - q(2,:)';
    dz = q(3,:) - q(3,:)';
    r2 = dx.^2 + dy.^2  +  dz.^2 + eye(n);

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
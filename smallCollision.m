clc; clear all; close all

%% generate particles
rng("default")
clc; clear
%% parameters of the simulation
L1 = 70; % L of the domain
L2 = 100; % H of the domain
epsilon = 20;
sigma = .5;
rCut = 1.25*sigma;
m = 1;
dt = 0.001;
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
d = cellfun (@numel, grd_to_ptcl, 'UniformOutput', true);

% number of time steps
tFinal = 4;
nTime = round(tFinal/dt);

%% functions to pass to the integrator
force = @(ptcls, grd_to_ptcl) forceCells(ptcls, grd_to_ptcl, epsilon, sigma, rCut^2); % this is -Force (negative)
boundaryConditions = @(ptcls) updateBoundaryConditions(ptcls, L1, L2, rCut, rCut);
ghost = @(ptcls, NP) updateGhost(ptcls, NP, L1, L2, rCut, rCut);
dKdp = @(p) p/m;

%% run the simulation
savingStep = 10;
[q, p] = sint.cellVelVerlet(force,dKdp,dt,nTime,grd,ptcls,grd_to_ptcl,boundaryConditions,ghost,savingStep);
%[q, p] = sint.cellForest(force,dKdp,dt,nTime,grd,ptcls,grd_to_ptcl,boundaryConditions,ghost,savingStep);

%% plot the results
figure
for i = 1:size(q,3)
    scatter(q(1,1:numel(X1),i), q(2,1:numel(X1),i), 10,"red" ,'filled')
    hold on
    scatter(q(1,numel(X1)+1:end,i), q(2,numel(X1)+1:end,i), 10,"blue" ,'filled')
    axis([0 L1 0 L2])
    xline(L1-rCut)
    xline(rCut)
    yline(L2-rCut)
    yline(rCut)
    drawnow
    hold off
end

%% plot the results p
figure
for i = 1:size(q,3)
    scatter(q(1,:,i), q(2,:,i), 10,vecnorm(p(:,:,i)) ,'filled')
    hold on
    axis([0 L1 0 L2])
    xline(L1-rCut)
    xline(rCut)
    yline(L2-rCut)
    yline(rCut)
    drawnow
    hold off
end

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

function F = forceCells(ptcls, grd_to_ptcl, epsij, sigmaij, rcut2)
    index = cellfun(@numel, grd_to_ptcl, 'UniformOutput', true);
    Fvectx = zeros(numel(ptcls.x)/2,1);
    Fvecty = zeros(numel(ptcls.x)/2,1);

    [nLy, nLx] = size(grd_to_ptcl);
    removeIndex = unique([1:nLy, 1:nLy:nLx*nLy, nLy:nLy:nLx*nLy, nLy*nLx-nLy+1:nLx*nLy]);
    nonEmpty = find(index>=1);
    nonEmpty = setdiff(nonEmpty,removeIndex);

    for i = nonEmpty(:)'
        % w-n-e-s-sw-nw-ne-su adiacent cells
        adCells = [i, i-nLy, i-1, i+nLy, i+1, i-nLy+1, i-nLy-1, i+nLy-1, i+nLy+1];
        indexPtclAd = grd_to_ptcl(adCells); indexPtclAd = vertcat(indexPtclAd{:});
        indexPtclLocal = grd_to_ptcl(i); indexPtclLocal = vertcat(indexPtclLocal{:});

        if isempty(indexPtclAd) == 1
            continue
        end
    
        % calculate pairwise distance between local and adiacent particles
        dx = ptcls.x(1,indexPtclAd) - ptcls.x(1,indexPtclLocal)';
        dy = ptcls.x(2,indexPtclAd) - ptcls.x(2,indexPtclLocal)';
        distanceMat2 = dx.^2 + dy.^2;

        % get pairwise inside/outside cut radius
        fc = find(distanceMat2 < rcut2 & distanceMat2 > eps);
        [fc1, ~] = ind2sub(size(distanceMat2),fc);
    
        if isempty(fc) == 1
            continue
        end

        r2Local = distanceMat2(fc);
        dx = dx(fc);
        dy = dy(fc);

        % calculate force between local and adiacent particles
        fc1Unique = unique(fc1);
        Fx = zeros(length(indexPtclLocal),1); Fy = Fx;
        [Fxv, Fyv] = lennardJonesForce(dx, dy, r2Local, sigmaij, epsij);

        % calc forces on localPtcls
        fc1Offsetted = discretize(fc1,length(fc1Unique));
        Fx(fc1Unique) = accumarray(fc1Offsetted(:),Fxv(:));
        Fy(fc1Unique) = accumarray(fc1Offsetted(:),Fyv(:));

        % sum forces
        Fvectx(indexPtclLocal) = Fvectx(indexPtclLocal) + Fx;
        Fvecty(indexPtclLocal) = Fvecty(indexPtclLocal) + Fy;
    end

    F = [Fvectx'; Fvecty'];

    %% internal functions
    function [Fx, Fy] = lennardJonesForce(dx, dy, r2, sigmaij, epsij)
        % calculate sigma2 - 6
        sigma2 = (sigmaij^2)./r2;
        sigma6 = sigma2.*sigma2.*sigma2;

        % calculate force
        Fmat = 48*epsij*sigma6.*(sigma6 - 0.5)./r2;
        Fx = Fmat.*dx;
        Fy = Fmat.*dy;
    end
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

function x = updateGhost(x, NP, L1, L2, hx, hy)
    if any(x(1,:) < 2*hx) || any(x(1,:) > L1-2*hx) ...
            || any(x(2,:) < 2*hy) || any(x(2,:) > L2-2*hx)
        % periodic boundary conditions using ghost cells
        % left boundary
        left = x(1,1:NP) < 2*hx;
        left = find(left);
        x = [x, [x(1,left)+L1-2*hx; x(2,left)]];

        % right boundary
        right = x(1,1:NP) > L1-2*hx;
        right = find(right);
        x = [x, [x(1,right)-L1+2*hx; x(2,right)]];

        % bottom boundary
        bottom = x(2,1:NP) < 2*hy;
        bottom = find(bottom);
        x = [x, [x(1,bottom); x(2,bottom)+L2-2*hy]];

        % top boundary
        top = x(2,1:NP) > L2-2*hy;
        top = find(top);
        x = [x, [x(1,top); x(2,top)-L2+2*hy]];
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
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
d = cellfun (@numel, grd_to_ptcl, 'UniformOutput', true);

% number of time steps
tFinal = 20;
nTime = round(tFinal/dt);

%% functions to pass to the integrator
force = @(ptcls, grd_to_ptcl) forceCells(ptcls, grd_to_ptcl, epsilon, sigma, rCut^2); % this is -Force (negative)
boundaryConditions = @(ptcls) updateBoundaryConditions(ptcls, L1, L2, rCut, rCut);
ghost = @(ptcls, NP) updateGhost(ptcls, NP, L1, L2, rCut, rCut);
dKdp = @(p) p/m;

%% run the simulation
savingStep = 1000;
[q, p] = sint.cellVelVerlet(force,dKdp,dt,nTime,grd,ptcls,grd_to_ptcl,boundaryConditions,ghost,savingStep);
%[q, p] = sint.cellForest(force,dKdp,dt,nTime,grd,ptcls,grd_to_ptcl,boundaryConditions,ghost,savingStep);

%% plot the results
figure
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
    hold off
end

%% plot the results p
figure
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
    hold off
end
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

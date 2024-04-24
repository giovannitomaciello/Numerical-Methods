clc; clear all; close all

%% distribute particles in a 2D domain as two colliding circles

% number of particles
N = 1000;

% radius of the circles
r = 20;

% center of the circles
x = [0 0];
y = [-50 50];

%% generate particles
grd.ncy = 100;
grd.ncx = 40;
grd.x = linspace (0, 250, grd.ncx+1);
grd.y = linspace (0, 200, grd.ncy+1);
[X1, Y1] = meshgrid (linspace (125-39/2*2^(1/6), 125+39/2*2^(1/6), 40),...
                 linspace (120-39/2*2^(1/6), 120+39/2*2^(1/6), 40));
[X2, Y2] = meshgrid (linspace (125-159/2*2^(1/6), 125+159/2*2^(1/6), 40),...
                 linspace (50-39/2*2^(1/6), 50+39/2*2^(1/6), 40));
ptcls.x = [[X1(:);X2(:)], [Y1(:);Y2(:)]]';
ptcls.v = randn (size (ptcls.x)) * .1;
ptcls.v(2, 1:numel(X1)) = -10;
grd_to_ptcl = sint.init_ptcl_mesh (grd, ptcls);

%% parameters of the simulation
L1 = 250; % length of the domain
L2 = 200; % width of the domain
epsilon = 5;
sigma = 1;
rCut = 2.5*sigma;
m = 1;
dt = 0.001;

% number of time steps
tFinal = 10;
nTime = round(tFinal/dt);

%% FUNCTIONS

function [Fvectx, Fvecty] = forceCells(Fvectx, Fvecty, grd, ptcls, grd_to_ptcl, epsij, sigmaij)
    index = cellfun(@numel, grd_to_ptcl, 'UniformOutput', true);

    nLy = numel(grd.y);
    removeIndex = unique([1:nLy, 1:nLy:nLx*nLy, nLy:nLy:nLx*nLy, nLy*nLx-nLy+1:nLx*nLy]);
    nonEmpty = find(index>=1);
    nonEmpty = setdiff(nonEmpty,removeIndex);

    for i = nonEmpty(:)'
    % w-n-e-s-sw-nw-ne-su adiacent cells
    adCells = [i-nLy, i-1, i+nLy, i+1, i-nLy+1, i-nLy-1, i+nLy-1, i+nLy+1];
    indexPtclAd = grd_to_ptcl(adCells); indexPtclAd = horzcat(indexPtclAd{:});
    indexPtclLocal = grd_to_ptcl(i); indexPtclLocal = horzcat(indexPtclLocal{:});

    if isempty(indexPtclAd) == 1
        continue
    end

    % calculate pairwise distance between local particles
    dx = ptcls.x(1,indexPtclLocal) - ptcls.x(1,indexPtclLocal)';
    dy = ptcls.x(2,indexPtclLocal) - ptcls.x(2,indexPtclLocal)';
    distanceMat2 = dx.^2 + dy.^2 + eye(numel(indexPtclLocal));

    % calculate pairwise force between local particles
    [Fx, Fy] = LennardJonesForceMatrix(dx, dy, distanceMat2, sigmaij, epsij);
    Fvectx(indexPtclLocal) = Fvectx(indexPtclLocal) + Fx;
    Fvecty(indexPtclLocal) = Fvecty(indexPtclLocal) + Fy;

    if isempty(indexPtclAd) == 1
        continue
    end
    
    % calculate pairwise distance between local and adiacent particles
    dx = ptcls.x(1,indexPtclAd) - ptcls.x(1,indexPtclLocal)';
    dy = ptcls.x(2,indexPtclAd) - ptcls.x(2,indexPtclLocal)';
    distanceMat2 = dx.^2 + dy.^2;

    % get pairwise inside/outside cut radius
    fc = find(distanceMat2 < rcut2);
    
    if isempty(fc) == 1
        continue
    end

    r2Local = distanceMat2(fc);
    dx = dx(fc);
    dy = dy(fc);

    % calculate force between local and adiacent particles
    [Fx, Fy] = lennardJonesForce(dx, dy, r2Local, sigmaij, epsij);
    Fmatx = zeros(numel(indexPtclLocal),numel(indexPtclAd));
    Fmaty = zeros(numel(indexPtclLocal),numel(indexPtclAd));

    % indexing forces
    Fmatx(fc) = Fx; Fx = sum(Fmatx,2);
    Fmaty(fc) = Fy; Fy = sum(Fmaty,2);

    % sum forces
    Fvectx(indexPtclLocal) = Fvectx(indexPtclLocal) + Fx;
    Fvecty(indexPtclLocal) = Fvecty(indexPtclLocal) + Fy;
    end

    %% internal functions
    function [Fx, Fy] = lennardJonesForce(dx, dy, r2, sigmaij, epsij)
        % calculate sigma2 - 6
        sigma2 = (sigmaij^2)./r2;
        sigma6 = sigma2.*sigma2.*sigma2;

        % calculate force
        Fmat = 48*epsij*sigma6.*(sigma6 - 0.5)./r2;
        Fx = -Fmat.*dx;
        Fy = -Fmat.*dy;
    end

    function [Fx, Fy] = LennardJonesForceMatrix(dx, dy ,r2 , sigmaij, epsij)
        % calculate sigma2 - 6
        sigma2 = (sigmaij^2)./r2;
        sigma6 = sigma2.*sigma2.*sigma2;

        % calculate force
        Fmat = 48*epsij*sigma6.*(sigma6 - 0.5)./r2;
        Fx = -Fmat.*dx;
        Fy = -Fmat.*dy;
    
        % sum forces
        Fx = sum(Fx,2);
        Fy = sum(Fy,2);
    end
end

function ptcls = updateBoundaryConditions(ptcls, L1, L2)
    % update boundary conditions
    ptcls.x(1, ptcls.x(1,:) > L1) = ptcls.x(1, ptcls.x(1,:) > L1) - L1;
    ptcls.x(1, ptcls.x(1,:) < 0) = ptcls.x(1, ptcls.x(1,:) < 0) + L1;
    ptcls.x(2, ptcls.x(2,:) > L2) = ptcls.x(2, ptcls.x(2,:) > L2) - L2;
    ptcls.x(2, ptcls.x(2,:) < 0) = ptcls.x(2, ptcls.x(2,:) < 0) + L2;
end
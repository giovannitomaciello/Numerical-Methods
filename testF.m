%% generate particles
rng("default")
clc; clear
%% parameters of the simulation
L1 = 20; % length of the domain
L2 = 20; % width of the domain
epsilon = 5;
sigma = 1;
rCut = 2.5*sigma;
m = 1;
dt = 0.00001;

%% generate particles
scale = 1;
delta = sigma*2^(1/6);

H1 = 10; W1 = 10;

H1_l = (H1-1)*delta; W1_l = (W1-1)*delta;

grd.ncy = L2/rCut;
grd.ncx = L1/rCut;
grd.x = linspace (0, L1, grd.ncx+1);
grd.y = linspace (0, L2, grd.ncy+1);

yc_1 = 10-20*(1-(1/scale));

[X1, Y1] = meshgrid (linspace ((L1-W1_l)/2, (L1+W1_l)/2, W1),...
		     linspace ((yc_1)-H1_l/2, (yc_1)+H1_l/2, H1)); 

ptcls.x = [[X1(:) + rand(size(X1(:)))*0.1], [Y1(:) + rand(size(X1(:)))*0.1]]';
ptcls.p = randn (size (ptcls.x)) * .1;
grd_to_ptcl = sint.init_ptcl_mesh (grd, ptcls);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

epsij = 5;
sigmaij = 1;

%% plot
% scatter(ptcls.x(1,:), ptcls.x(2,:),80,1:size(ptcls.x,2),"filled")
% colormap(hsv(size(ptcls.x,2)))
% colorbar

mesh (grd.x, grd.y, zeros (numel(grd.y), numel(grd.x)))
view (2)
hold on
for jj = 1 : grd.ncx
  for ii = 1 : grd.ncy
    plot (ptcls.x(1, grd_to_ptcl{ii,jj}),...
          ptcls.x(2, grd_to_ptcl{ii,jj}), 'x')
    axis ([grd.x(1) grd.x(end) grd.y(1) grd.y(end)])
    pause (.1)
  end
end

%% start test 1
index = cellfun(@numel, grd_to_ptcl, 'UniformOutput', true);
Fvectx = zeros(numel(ptcls.x)/2,1);
Fvecty = zeros(numel(ptcls.x)/2,1);

[nLx, nLy] = size(grd_to_ptcl);
removeIndex = unique([1:nLy, 1:nLy:nLx*nLy, nLy:nLy:nLx*nLy, nLy*nLx-nLy+1:nLx*nLy]);
nonEmpty = find(index>=1);
nonEmpty = setdiff(nonEmpty,removeIndex);
rcut2 = min(grd.x(2) - grd.x(1),grd.y(2) - grd.y(1))^2;

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

Ff = [Fvectx'; Fvecty'];

%% start test Matrix
Fm = LennardJonesForceM(ptcls.x', sigmaij, epsij);

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

%% check Function

function F = LennardJonesForceM(q, sigmaij, epsij)
    n = length(q);
    F = zeros(2,n);

    % calculate distance
    dx = q(:,1) - q(:,1)';
    dy = q(:,2) - q(:,2)';
    r2 = dx.^2 + dy.^2 + eye(n);

    % calculate sigma2 - 6
    sigma2 = (sigmaij^2)./r2;
    sigma6 = sigma2.*sigma2.*sigma2;

    % calculate force
    Fmat = 48*epsij*sigma6.*(sigma6 - 0.5)./r2;
    Fx = -Fmat.*dx;
    Fy = -Fmat.*dy;
    
    % sum forces
    F(1,:) = sum(Fx,2)';
    F(2,:) = sum(Fy,2)';
end



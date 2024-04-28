clc; clear all; close all

%% generate particles
rng("default")
clc; clear
%% parameters of the simulation
L1 = 40; % length of the domain
L2 = 70; % width of the domain
epsilon = 5;
sigma = 1;
rCut = 2.5*sigma;
m = 1;
dt = 0.0001;

%% generate two squares colliding
scale = 1;
delta = sigma*2^(1/6);

H1 = 4; W1 = 4;
H2 = 4; W2 = 4;

H1_l = (H1-1)*delta; W1_l = (W1-1)*delta;
H2_l = (H2-1)*delta; W2_l = (W2-1)*delta;

grd.ncy = L2/rCut;
grd.ncx = L1/rCut;
grd.x = linspace (0, L1, grd.ncx+1);
grd.y = linspace (0, L2, grd.ncy+1);

yc_1 = 50-20*(1-(1/scale));
yc_2 = 20-20*(1-(1/scale));

[X1, Y1] = meshgrid (linspace ((L1-W1_l)/2, (L1+W1_l)/2, W1),...
             linspace ((yc_1)-H1_l/2, (yc_1)+H1_l/2, H1));
[X2, Y2] = meshgrid (linspace ((L1-W2_l)/2, (L1+W2_l)/2, W2),...
                linspace ((yc_2)-H2_l/2, (yc_2)+H2_l/2, H2));

ptcls.x = [[X1(:);X2(:)], [Y1(:);Y2(:)]]';
ptcls.p = randn (size (ptcls.x)) * .01;
grd_to_ptcl = sint.init_ptcl_mesh (grd, ptcls);
d = cellfun (@numel, grd_to_ptcl, 'UniformOutput', true);
% initialize the momentum
ptcls.p(:,1:numel(X1)) = [0; -5]*ones(1,numel(X1));
ptcls.p(:,numel(X1)+1:end) = [0;5]*ones(1,numel(X2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

epsij = 5;
sigmaij = 1;

%% plot
% scatter(ptcls.x(1,:), ptcls.x(2,:),80,1:size(ptcls.x,2),"filled")
% colormap(hsv(size(ptcls.x,2)))
% colorbar

% mesh (grd.x, grd.y, zeros (numel(grd.y), numel(grd.x)))
% view (2)
% hold on
% for jj = 1 : grd.ncx
%   for ii = 1 : grd.ncy
%     plot (ptcls.x(1, grd_to_ptcl{ii,jj}),...
%           ptcls.x(2, grd_to_ptcl{ii,jj}), 'Marker','square')
%     axis ([grd.x(1) grd.x(end) grd.y(1) grd.y(end)])
%     pause (.1)
%   end
% end

%% start test 1
index = cellfun(@numel, grd_to_ptcl, 'UniformOutput', true);
Fvectx = zeros(numel(ptcls.x)/2,1);
Fvecty = zeros(numel(ptcls.x)/2,1);

[nLy, nLx] = size(grd_to_ptcl);
removeIndex = unique([1:nLy, 1:nLy:nLx*nLy, nLy:nLy:nLx*nLy, nLy*nLx-nLy+1:nLx*nLy]);
nonEmpty = find(index>=1);
nonEmpty = setdiff(nonEmpty,removeIndex);
rcut2 = min(grd.x(2) - grd.x(1),grd.y(2) - grd.y(1))^2;

for i = 1:20000
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
    fc1 = rem(fc-1,size(distanceMat2,1))+1;

    if isempty(fc) == 1
        continue
    end

    r2Local = distanceMat2(fc);
    dx = dx(fc);
    dy = dy(fc);

    % calculate force between local and adiacent particles
    Fx = zeros(length(indexPtclLocal),1); Fy = Fx;
    [Fxv, Fyv] = lennardJonesForce(dx, dy, r2Local, sigmaij, epsij);

    % calc forces on localPtcls
    offset = min(fc1) - 1;
    maxFc1 = max(fc1);
    fc1Offsetted = fc1 - offset;
    Fx(offset+1:maxFc1) = accumarray(fc1Offsetted(:),Fxv(:));
    Fy(offset+1:maxFc1) = accumarray(fc1Offsetted(:),Fyv(:));

    % sum forces
    Fvectx(indexPtclLocal) = Fvectx(indexPtclLocal) + Fx;
    Fvecty(indexPtclLocal) = Fvecty(indexPtclLocal) + Fy;
end

Ff = [Fvectx'; Fvecty'];
end

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



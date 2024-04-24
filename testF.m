%% generate particles
rand('seed',1)
NP = 100;
grd.ncy = 100;
grd.ncx = 100;
hx = 1/(grd.ncx);
hy = 1/(grd.ncy);
grd.x = linspace (0, 1, grd.ncx);
grd.y = linspace (0, 1, grd.ncy);
X = rand(1,NP);
Y = rand(1,NP);
ptcls.x = [X(:), Y(:)]';
ptcls.v = randn (size (ptcls.x)) * .1;
grd_to_ptcl = sint.init_ptcl_mesh(grd, ptcls);

index = cellfun(@numel, grd_to_ptcl, 'UniformOutput', true);

nLx = numel(grd.x);
nLy = numel(grd.y);
removeIndex = unique([1:nLy, 1:nLy:nLx*nLy, nLy:nLy:nLx*nLy, nLy*nLx-nLy+1:nLx*nLy]);
nonEmpty = find(index>=1);
nonEmpty = setdiff(nonEmpty,removeIndex);
rcut2 = min(hx,hy)^2; 

epsij = 0.1;
sigmaij = 0.1;

Fvectx = zeros(NP,1);
Fvecty = zeros(NP,1);

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




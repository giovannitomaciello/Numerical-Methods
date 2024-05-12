function F = forceCells(forceCalculator, ptcls, grd_to_ptcl, rcut2)
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
        fc1 = rem(fc-1,size(distanceMat2,1))+1;
    
        if isempty(fc) == 1
            continue
        end

        r2Local = distanceMat2(fc);
        dx = dx(fc);
        dy = dy(fc);

        % calculate force between local and adiacent particles
        Fx = zeros(length(indexPtclLocal),1); Fy = Fx;
        [Fxv, Fyv] = forceCalculator(dx, dy, r2Local);

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
    F = [Fvectx'; Fvecty'];
end
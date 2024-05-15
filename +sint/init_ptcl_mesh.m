function grd_to_ptcl = init_ptcl_mesh (grd, ptcls)
    dim = size(ptcls.x,1);
    switch dim
        case 2
            h = [grd.x(2) - grd.x(1); grd.y(2) - grd.y(1)];
            idx = floor(ptcls.x./h) + 1;
            indexPtcls = 1:size(ptcls.x,2);
            grd_to_ptcl = accumarray(idx',indexPtcls(:),[grd.ncx grd.ncy],@(x) {x});
        case 3
            h = [grd.x(2) - grd.x(1); grd.y(2) - grd.y(1);  grd.z(2) - grd.z(1)];
            idx = floor(ptcls.x./h) + 1;
            indexPtcls = 1:size(ptcls.x,2);
            grd_to_ptcl = accumarray(idx',indexPtcls(:),[grd.ncx grd.ncy grd.ncz],@(x) {x});
    end
end
function grd_to_ptcl = init_ptcl_mesh (grd, ptcls)
  h = [grd.x(2) - grd.x(1); grd.y(2) - grd.y(1)];
  idx = floor(ptcls.x./h) + 1;
  indexPtcls = 1:size(ptcls.x,2);

  %grd_to_ptcl = cell (grd.ncy,grd.ncx);
  %for ii = 1 : numel (ptcl_to_grd)
    %grd_to_ptcl{ptcl_to_grd(ii)}(end+1) = ii;
  %end
  grd_to_ptcl = accumarray(idx',indexPtcls(:),[grd.ncx grd.ncy],@(x) {x});
end
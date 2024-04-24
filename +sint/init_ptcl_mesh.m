function grd_to_ptcl = init_ptcl_mesh (grd, ptcls)
  hx = grd.x(2) - grd.x(1); hy = grd.y(2) - grd.y(1); 
  ridx = floor(ptcls.x(2, :)/hy) + 1;
  cidx = floor(ptcls.x(1, :)/hx) + 1;

  ptcl_to_grd = sub2ind ([grd.ncy, grd.ncx], ridx, cidx);
  grd_to_ptcl = cell (grd.ncy,grd.ncx);
  for ii = 1 : numel (ptcl_to_grd)
    grd_to_ptcl{ptcl_to_grd(ii)}(end+1) = ii;
  end
end
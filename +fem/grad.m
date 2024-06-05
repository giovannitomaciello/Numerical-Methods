function [gradx,grady,gradz] = grad(field3,hx,hy,hz)
    gradx = diff(field3,[],1)*hx;
    grady = diff(field3,[],2)*hy;
    gradz = diff(field3,[],3)*hz;
end
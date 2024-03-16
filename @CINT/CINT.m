classdef CINT
    methods (Static)
        [q, p] = rattle(q0,p0,dTdq,dKdp,C,t);
    end

    methods (Static)
        q = Shake(q0,p0,dTdq,dKdp,G,C,m,t);
    end
end
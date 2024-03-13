classdef INT
    methods (Static)
        [q, p] = euleroavanti(q0,p0,dTdq,dKdp,t)
    end
    methods (Static)
        [q, p] = euleroindietro(q0,p0,dTdq,dKdp,t)
    end
    methods (Static)
        [q, p] = crankNick(q0,p0,dTdq,dKdp,t)
    end
    methods (Static)
        [q,p] = velVerlet(q0,p0,dTdq,dKdp,t)
    end
    methods (Static)
        [q,p] = posVerlet(q0,p0,dTdq,dKdp,m,t)
    end
    methods (Static)
        [q, p] = symplecticEuler(q0,p0,dTdq,dKdp,t)
    end
end
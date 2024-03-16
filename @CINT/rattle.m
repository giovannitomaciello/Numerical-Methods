function [q,p] = rattle(q0,p0,dTdq,dKdp,C,t)

[NP, ND] = size(q0);

NT = numel(t);

q = zeros(NP,ND,NT);
p = zeros(NP,ND,NT);

q(:,:,1) = q0;
p(:,:,1) = p0;

[Index] = find(C);

%Ftmp ha dimensioni [NP,ND] e Ftmp(i,j) rappresenta la componente j sulla particella i  
dt = diff(t);
sparseM = sparse(NP, NP);
Ftmp = -dTdq(q0);

for i = 2:NT
    p0 = p(:,:,i-1);
    q0 = q(:,:,i-1);

    UNK = fsolve(@(UNK) SYS1(UNK, q0, p0, Ftmp, dKdp, C, Index, sparseM, NP, ND, ... 
    dt(i-1)), [q0(:); p0(:), zeros(numel(lambda),1)]);
    
    q(:,:,i) = reshape(UNK(1:NP*ND), NP, []);
    ptmp = reshape(UNK(NP*ND+1:2*NP*ND), NP, []);
    Ftmp = -dTdq(q(:,:,i));

    UNK = fsolve(@(UNK) SYS2(UNK, q0, p0, dTdq, dKdp, C, Index, sparseM, NP, ND, ...
    dt(i-1)), [p0(:), zeros(numel(lambda),1)]);

    p(:,:,i) = reshape(UNK(1:NP*ND), NP, []);
end
end

function toZero = SYS1(UNK, q0, p0, F, dKdp, C, Index, sparseM, NP, ND, dt)
    q = reshape(UNK(1:NP*ND), NP, []);
    p = reshape(UNK(NP*ND+1:2*NP*ND), NP, []);
    lambda = UNK(2*NP*ND+1:end);

    [Gout, qCM] = G(q0, lambda, Index, NP, sparseM);

    toZeroq = - q + q0 + dt*p;
    toZerop = - p + p0 - dt/2*(-Ftmp + Gout.*lambda);
    toZerolambda = - nonzeros(C(qCM)) + nonzeros(C).^2;

    toZero = [toZeroq(:); toZerop(:); toZerolambda(:)];
end

function toZero = SYS2(UNK, q0, p0, dTdq, dKdp, C, Index, sparseM, NP, ND, dt)
    p = reshape(UNK(1:NP*ND), NP, []);
    lambda = UNK(NP*ND+1:end);

    [Gout, qCM] = G(q0, lambda, Index, NP, sparseM);

    toZerop = - p + p0 - dt/2*(-Ftmp + Gout.*lambda);
    toZerolambda = ???;

end

function [Gout, qCM] = G(q, lambda, Index, NP, sparseM)
    ???
end
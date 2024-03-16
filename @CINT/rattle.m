function [q,p] = rattle(q0,p0,dTdq,dKdp,G,C,t)

[NP, ND] = size(q0);

NT = numel(t);
NC = numel(nonzeros(C));

q = zeros(NP,ND,NT);
p = zeros(NP,ND,NT);

q(:,:,1) = q0;
p(:,:,1) = p0;

[Index] = find(C);
[Indexr, ~] = ind2sub(NP,Index);
C2 = C.^2;

%Ftmp ha dimensioni [NP,ND] e Ftmp(i,j) rappresenta la componente j sulla particella i  
dt = diff(t);
sparseM = sparse(NP, NP);
Ftmp = -dTdq(q0);

for i = 2:NT
    p0 = p(:,:,i-1);
    q0 = q(:,:,i-1);

    UNK = fsolve(@(UNK) SYS1(UNK, q0, p0, Ftmp, dKdp, C2, G, Index, Indexr, sparseM, NP, ND, ... 
    dt(i-1)), [q0(:); p0(:); zeros(NC,1)]*(1 + 1e-6));
    
    q(:,:,i) = reshape(UNK(1:NP*ND), NP, []);
    ptmp = reshape(UNK(NP*ND+1:2*NP*ND), NP, []);
    Ftmp = -dTdq(q(:,:,i));

    UNK = fsolve(@(UNK) SYS2(UNK, q0, ptmp, Ftmp, dKdp, C2, G, Index, Indexr, sparseM, NP, ND, ...
    dt(i-1)), [p0(:); zeros(NC,1)]*(1 + 1e-6));

    p(:,:,i) = reshape(UNK(1:NP*ND), NP, []);
end
end

function toZero = SYS1(UNK, q0, p0, Ftmp, dKdp, C2, G, Index, Indexr, sparseM, NP, ND, dt)
    % preprocessing
    q = reshape(UNK(1:NP*ND), NP, []);
    p = reshape(UNK(NP*ND+1:2*NP*ND), NP, []);
    lambda = UNK(2*NP*ND+1:end);
    lambdaV = zeros(NP,1); lambdaV(Indexr) = lambda;

    % evaluate prop
    [Gq,r2] = G(q,lambdaV);

    % to zero vectors
    toZeroq = - q + q0 + dt*dKdp(p);
    toZerop = - p + p0 - dt/2*(-Ftmp + Gq.*lambdaV);
    toZerolambda = r2(Index) - nonzeros(triu(C2));
    
    % all to zero
    toZero = [toZeroq(:); toZerop(:); toZerolambda(:)]
end

function toZero = SYS2(UNK, q0, p0, Ftmp, dKdp, C2, G, Index, Indexr, sparseM, NP, ND, dt)
    % preprocessing
    p = reshape(UNK(1:NP*ND), NP, []);
    lambda = UNK(NP*ND+1:end);
    lambdaV = zeros(NP,1); lambdaV(Indexr) = lambda;

    % evaluate prop
    [Gq,r2] = G(q0,lambdaV);

    % to zero vectors
    toZerop = - p + p0 - dt/2*(-Ftmp + Gq.*lambdaV);
    toZerolambda = (r2(Index) - nonzeros(triu(C2)))*norm(dKdp(p(Indexr,:)),2);

    % all to zero
    toZero = [toZerop(:); toZerolambda(:)];
end
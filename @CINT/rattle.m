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

for i = 2:NT
    p0 = p(:,:,i-1);
    q0 = q(:,:,i-1);

    
end
end

function toZero = SYS1(UNK, q0, p0, F, dKdp, C, Index, sparseM, NP, ND, dt)

end

function toZero = SYS2(UNK, q0, p0, q, F, dKdp, C, Index, sparseM, NP, ND, dt)

end

function [Gout, qCM] = G(q, lambda, Index, NP, sparseM)

end
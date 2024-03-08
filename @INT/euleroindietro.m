function [q, p] = euleroindietro(q0,p0,dTdq,dKdp,t)

[NP, ND] = size(q0);

NT = numel(t);

q = zeros(NP,ND,NT);
p = zeros(NP,ND,NT);

q(:,:,1) = q0;
p(:,:,1) = p0;

for i=2:NT
    % time step
    dt = t(i) - t(i-1);
    % init val
    q0 = q(:,:,i-1);  
    p0 = p(:,:,i-1);
    
    opt = optimoptions("fsolve","Display","none");
    unk = fsolve(@(unk) sysEB(unk, q0, p0, dt, dKdp, dTdq),....
        [q0(:); p0(:)], opt);

    q(:,:,i) = reshape(unk(1:NP*ND,:),[],3);
    p(:,:,i) = reshape(unk(NP*ND+1:end,:),[],3);

end
end

function toZero = sysEB(unk, q0, p0, dt, dKdp, dTdq)

    [NP, ND] = size(q0);
    q = reshape(unk(1:NP*ND,:),[],3);
    p = reshape(unk(NP*ND+1:end,:),[],3);
    
    % position
    zeroq = dKdp(p) - (q - q0)/dt;

    %momentum
    zerop = - dTdq(q) - (p - p0)/dt;

    toZero = [zeroq(:); zerop(:)];
end
function [q, p] = crankNick(q0,p0,dTdq,dKdp,t)

[NP, ND] = size(q0);

NT = numel(t);

q = zeros(NP,ND,NT);
p = zeros(NP,ND,NT);

q(:,:,1) = q0;
p(:,:,1) = p0;

unk = [q0(:); p0(:)];
dt = diff(t);

for i=2:NT
    % init val
    q0 = q(:,:,i-1);  
    p0 = p(:,:,i-1);
    
    opt = optimoptions("fsolve","Display","none","Algorithm","levenberg-marquardt");
    unk = fsolve(@(unk) sysEB(unk, q0, p0, dt(i-1), dKdp, dTdq),...
        unk, opt);

    q(:,:,i) = reshape(unk(1:NP*ND,:),[],3);
    p(:,:,i) = reshape(unk(NP*ND+1:end,:),[],3);
end
end

function toZero = sysEB(unk, q0, p0, dt, dKdp, dTdq)

    [NP, ND] = size(q0);
    q = reshape(unk(1:NP*ND,:),[],3);
    p = reshape(unk(NP*ND+1:end,:),[],3);

    % position
    zeroq = (dKdp(p) + dKdp(p0))*dt/2 - q + q0;

    %momentum
    zerop = (dTdq(q) + dTdq(q0))*dt/2 + p - p0;

    toZero = [zeroq(:); zerop(:)];
end
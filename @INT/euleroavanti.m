function [q, p] = euleroavanti(q0,p0,dTdq,dKdp,t)

[NP, ND] = size(q0);

NT = numel(t);

q = zeros(NP,ND,NT);
p = zeros(NP,ND,NT);

q(:,:,1) = q0;
p(:,:,1) = p0;

for i=2:length(t)
    % time step
    dt = t(i) - t(i-1);
    
    % position
    q(:,:,i) = dKdp(p(:,:,i-1))*dt + q(:,:,i-1);

    %momentum
    p(:,:,i) = - dTdq(q(:,:,i-1))*dt + p(:,:,i-1);

end
end
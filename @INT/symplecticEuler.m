function [q, p] = symplecticEuler(q0,p0,dTdq,dKdp,t)

[NP, ND] = size(q0);

NT = numel(t);
t = [t, (t(end) - t(end-1) + t(end))];

q = zeros(NP,ND,NT);
p = zeros(NP,ND,NT);

q(:,:,1) = q0;
p(:,:,1) = p0;
dt = diff(t);

for i=2:NT
    % position
    q(:,:,i) = q(:,:,i-1) + dKdp(p(:,:,i-1))*dt(i-1);

    %momentum
    p(:,:,i) = p(:,:,i-1) - dTdq(q(:,:,i))*dt(i);

end
end
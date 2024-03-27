function [q,p] = posVerlet(q0,p0,dTdq,dKdp,m,t)

[NP, ND] = size(q0);
NT = numel(t);
t = [t, (t(end) - t(end-1) + t(end))];
q = zeros(NP,ND,NT+1);
p = zeros(NP,ND,NT);

q(:,:,1) = q0;
p(:,:,1) = p0;

%Force tmp

%Ftmp ha dimensioni [NP,ND] e Ftmp(i,j) rappresenta la componente j sulla particella i
Ftmp = - dTdq(q0);
i = 2;
dt = t(i) - t(i-1);

%momentum tmp
ptmp = Ftmp*dt/2 + p(:,:,i-1);

% position
q(:,:,i) = dKdp(ptmp)*dt + q(:,:,i-1);


%Force tmp
Ftmp = - dTdq(q(:,:,i));

%momentum
p(:,:,i) = Ftmp*dt/2 + ptmp;

for i = 3:NT+1
    % time step
    dt = t(i) - t(i-1);

    % position
    q(:,:,i) = 2*q(:,:,i-1) - q(:,:,i-2) - dt^2 * dTdq(q(:,:,i-1))./m;

    %momentum
    p(:,:,i-1) = m.*(q(:,:,i) - q(:,:,i-2))/2/dt;
end

q(:,:,end) = [];
end
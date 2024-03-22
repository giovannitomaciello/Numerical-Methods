function q = shake(q0,p0,dTdq,dKdp,G,C,m,t)

[NP, ND] = size(q0);
NT = numel(t);

ind_constraints = find(triu(C));

q = zeros(NP,ND,NT);
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

unk = C(ind_constraints);
for i = 3:NT
    % time step
    dt = t(i) - t(i-1);

    % position
    q_tilde = 2*q(:,:,i-1) - q(:,:,i-2) - dt^2 * dTdq(q(:,:,i-1))./m;
    q_prec = q(:,:,i-1);

    opt = optimoptions("fsolve","Display","none","OptimalityTolerance",1e-18,...
        "FunctionTolerance",1e-18);
    [unk,iszeros] = fsolve(@(unk) sysEB(unk,NP,q_tilde,q_prec,G,C,m,ind_constraints,dt),unk, opt);

    lambda = zeros(NP,NP);
    lambda(ind_constraints) = unk;
    lambda =lambda+lambda';

    [Gq,r2] = G(q_prec,lambda);

    q(:,:,i)=q_tilde-dt^2*Gq./m;

end
end

function toZero = sysEB(unk,NP,q_tilde,q_prec,G,C,m,ind_constraints,dt)
             
        lambda = zeros(NP,NP);
        lambda(ind_constraints) = unk;
        lambda =lambda+lambda';

        [Gq,~] = G(q_prec,lambda);
        q = q_tilde-dt^2*Gq./m;

        dx = q(:,1) - q(:,1)';
        dy = q(:,2) - q(:,2)';
        dz = q(:,3) - q(:,3)';
        r = dx.^2 + dy.^2 + dz.^2; 
       
        toZero = (r(ind_constraints) - C(ind_constraints).^2);

end
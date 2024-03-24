function [q,p] = rattle2(q0,p0,dTdq,dKdp,G,C,m,t)

[NP, ND] = size(q0);
NT = numel(t);

ind_constraints = find(triu(C));


q = zeros(NP,ND,NT);
p = zeros(NP,ND,NT);

q(:,:,1) = q0;
p(:,:,1) = p0;


unk1 = C(ind_constraints);
unk2 = C(ind_constraints);
for i = 2:NT
    % time step
    dt = t(i) - t(i-1);
    
    p_prec = p(:,:,i-1);
    q_prec = q(:,:,i-1);
    opt = optimoptions("fsolve","Display","none","OptimalityTolerance",1e-20,...
        "FunctionTolerance",1e-20);
    unk1 = fsolve(@(unk) sys1(unk,NP,p_prec,q_prec,G,C,m,ind_constraints,dt,dTdq),unk1,opt);
    
    lambda1 = zeros(NP,NP);
    lambda1(ind_constraints) = unk1;
    lambda1 =lambda1+lambda1';
    [G1,~] = G(q_prec,lambda1);
    ptmp = p(:,:,i-1) - dt/2*(dTdq(q(:,:,i-1)) + G1 );
    q(:,:,i) = q_prec + dt*ptmp./m;
 

    unk2 = fsolve(@(unk) sys2(unk,NP,ptmp,q(:,:,i),G,G1,ind_constraints,dt,dKdp,dTdq),unk2,opt);
    lambda2 = zeros(NP,NP);
    lambda2(ind_constraints) = unk2;
    lambda2 =lambda2+lambda2';
    [G2,~] = G(q(:,:,i),lambda2);
    p(:,:,i) = ptmp - dt/2 * (dTdq(q(:,:,i)) + G1 + G2);

end
end

function toZero = sys1(unk,NP,p_prec,q_prec,G,C,m,ind_constraints,dt,dTdq)
             
        lambda = zeros(NP,NP);
        lambda(ind_constraints) = unk;
        lambda =lambda+lambda';

        [Gq,~] = G(q_prec,lambda);
        ptmp = p_prec - dt/2*(dTdq(q_prec) + Gq);
        q = q_prec + dt*ptmp./m;

        dx = q(:,1) - q(:,1)';
        dy = q(:,2) - q(:,2)';
        dz = q(:,3) - q(:,3)';
        r = dx.^2 + dy.^2 + dz.^2; 
       
        toZero = r(ind_constraints) - C(ind_constraints).^2;

end
function toZero = sys2(unk,NP,ptmp,q,G,G1,ind_constraints,dt,dKdp,dTdq)
             
        lambda = zeros(NP,NP);
        lambda(ind_constraints) = unk;
        lambda =lambda+lambda';

        [G2,~] = G(q,lambda);
        p = ptmp - dt/2 * (dTdq(q) + G1 + G2);

        dqdt = dKdp(p);
        dx = 2*(dqdt(:,1) - dqdt(:,1)').*(q(:,1) - q(:,1)'); 
        dy = 2*(dqdt(:,2) - dqdt(:,2)').*(q(:,2) - q(:,2)'); 
        dz = 2*(dqdt(:,3) - dqdt(:,3)').*(q(:,3) - q(:,3)');
               
        toZero = dx(ind_constraints) + dy(ind_constraints) + dz(ind_constraints);
end
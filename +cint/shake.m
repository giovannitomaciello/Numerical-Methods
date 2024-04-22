function [q,p] = shake(q0,p0,dTdq,dKdp,G,S,t,m)

[NP, ND] = size(q0);
NT = numel(t);

q = zeros(NP,ND,NT);
p = zeros(NP,ND,NT);

q(:,:,1) = q0;
p(:,:,1) = p0;

%Second step velocity verlet
Ftmp = - dTdq(q0);
dt = t(2) - t(1);
ptmp = Ftmp*dt/2 + p(:,:,1);
q(:,:,2) = dKdp(ptmp)*dt + q(:,:,1);
p(:,:,2) =  ptmp - dTdq(q(:,:,2))*dt/2;

%
S0 = S(q0); NC = numel(S0);
unk = zeros(NC,1); 

opt = optimoptions("fsolve","Display","none","OptimalityTolerance",1e-20,...
        "FunctionTolerance",1e-20,"FiniteDifferenceType","central","StepTolerance",1e-20);

for i = 3:NT
    % time step
    dt = t(i) - t(i-1);

    % position
    q_tilde = 2*q(:,:,i-1) - q(:,:,i-2) - dt^2 * dTdq(q(:,:,i-1))./m;
    q_prec = q(:,:,i-1);

    [lambda,iszeros1] = fsolve(@(lambda) sysEB(lambda,G(q_prec),q_tilde,m,dt,S),unk, opt);
    if any(abs(iszeros1) > 1e-8)
        warning(strcat("CONSTRAINTS NOT RESPECTED, value of sum(zeros) in SYS1 is:",num2str(sum(abs(iszeros1)))," time =",num2str(t(i))))
    end

    Gq = tensorprod(G(q_prec),lambda,2,1);
    q(:,:,i) = q_tilde + dt^2 * Gq./m;
    p(:,:,i) =  p(:,:,i-1) + (-dTdq(q(:,:,i-1))+Gq)*dt;  

end
end

function toZero = sysEB(lambda,Gq0,q_tilde,m,dt,S)
             
        Gq = tensorprod(Gq0,lambda,2,1);
        q = q_tilde + dt^2 * Gq./m;
        toZero = S(q);

end
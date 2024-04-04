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
p(:,:,2) =  p(:,:,1)- dTdq(q(:,:,2))*dt/2;

%
S0 = S(q0); NC = numel(S0);
unk = zeros(NC,1); 

opt = optimoptions("fsolve","Display","none","OptimalityTolerance",1e-20,...
        "FunctionTolerance",1e-20,"FiniteDifferenceType","central","StepTolerance",1e-18);
for i = 3:NT
    % time step
    dt = t(i) - t(i-1);

    % position
    q_tilde = 2*q(:,:,i-1) - q(:,:,i-2) - dt^2 * dTdq(q(:,:,i-1))./m;
    q_prec = q(:,:,i-1);

    [lambda,iszeros] = fsolve(@(lambda) sysEB(lambda,G(q_prec),q_tilde,m,dt,S),unk, opt);
    if any(abs(iszeros) > 1e-8)
        warning(strcat("CONSTRAINTS NOT RESPECTED, value of sum(zeros) in SYS1 is:",strcat(num2str(sum(abs(iszeros))))))
    end

    Gq = reshape(G(q_prec)*lambda,[],3);
    q = q_tilde - dt^2 * Gq./m;

    q(:,:,i) = q_tilde - dt^2 * Gq./m;
    p(:,:,i) =  p(:,:,i-1)- dTdq(q(:,:,i))*dt/2;

end
end

function toZero = sysEB(lambda,Gq0,q_tilde,m,dt,S)
             
        Gq = reshape(Gq0*lambda,[],3);
        q = q_tilde - dt^2 * Gq./m;
        toZero = S(q);
end
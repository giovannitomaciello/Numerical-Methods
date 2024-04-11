function [q,p] = rattle(q0,p0,dTdq,dKdp,G,S,t,dGdt)

q0 = q0(:);
p0 = p0(:);

N = numel(q0);
NT = numel(t);

q = zeros(N,NT);
p = zeros(N,NT);

q(:,1) = q0;
p(:,1) = p0;

S0 = S(q0); NC = numel(S0);
unk2 = zeros(NC,1);

opt = optimoptions("fsolve","Display","none","OptimalityTolerance",1e-20,...
        "FunctionTolerance",1e-20,"FiniteDifferenceType","central","StepTolerance",1e-18);

for i = 2:NT
    % time step
    dt = t(i) - t(i-1);
    
    p0 = p(:,i-1);
    q0 = q(:,i-1);
    dTdq0 = dTdq(q0);

    % explicit halved first try
    F = - dTdq0;
    ptmp = p0 + dt/2*F;
    qtry = dKdp(ptmp)*dt/2 + q(:,i-1);

    [lambda_q,iszeros1] = fsolve(@(X) sys1(X,p0,q0,G,S,dt,dTdq0,dKdp),[qtry; zeros(NC,1)],opt);
    if any(abs(iszeros1) > 1e-8)
        warning(strcat("CONSTRAINTS NOT RESPECTED, value of sum(zeros) in SYS1 is:",num2str(sum(abs(iszeros1)))," time =",num2str(t(i))))
    end
    
    q(:,i) = lambda_q(1:numel(q0));
    lambda = lambda_q(numel(q0)+1:end);
    Gq = G(q(:,i))*lambda;
    F = Gq - dTdq0; 
    ptmp = p0 + dt/2*F;
 
    % @q
    F =  G(q(:,i))*lambda - dTdq(q(:,i)); 

    [p_mu,iszeros2] = fsolve(@(X) sys2(X,ptmp,G(q(:,i)),F,dt,dKdp,dGdt),[ptmp;0*lambda],opt);
    if any(abs(iszeros2) > 1e-8)
        warning(strcat("CONSTRAINTS NOT RESPECTED, value of sum(zeros) in SYS2 is:",strcat(num2str(sum(abs(iszeros2))))))
    end

    p(:, i) = p_mu(1:N);
    %Frv = G(q(:,i))*mu + F;   
    %p(:,i) = ptmp + dt/2 * Frv;
end
end

function toZero = sys1(lambda_q,p0,q0,G,S,dt,dTdq0,dKdp)
        
        q = lambda_q(1:numel(q0));
        lambda = lambda_q(numel(q0)+1:end);

        Gq = G(q)*lambda;
        F = Gq - dTdq0; 
        ptmp = p0 + dt/2*F;
        toZero = q0 + dt*dKdp(ptmp)-q;
 
        toZero = [toZero;S(q)];
        
end

function toZero = sys2(p_mu,ptmp,Gq,F,dt,dKdp,dGdt)

        p = p_mu(1:numel(ptmp));
        mu = p_mu(numel(ptmp)+1:end);
        Frv = Gq*mu + F;   
        toZero = p - ptmp - dt/2 * Frv;
        toZero = [toZero; Gq'*dKdp(p)];
        % toZero = dGdt(Gq,dKdp(p));
end
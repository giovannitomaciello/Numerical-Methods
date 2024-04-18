function [q,p] = rattleNew(q0,p0,dTdq,dKdp,G,S,t,dGdt)

NP = size(q0,1);
NT = numel(t);

q = zeros(NP,3,NT);
p = zeros(NP,3,NT);

q(:,:,1) = q0;
p(:,:,1) = p0;

opt = optimoptions("fsolve","Display","none","OptimalityTolerance",1e-20,...
        "FunctionTolerance",1e-20,"FiniteDifferenceType","central","StepTolerance",1e-18);
NC = numel(S(q0));

for i = 2:NT
    % time step
    dt = t(i) - t(i-1);
    
    p0 = p(:,:,i-1);
    q0 = q(:,:,i-1);
    dTdq0 = dTdq(q0);

    % explicit halved first try
    % F = - dTdq0;
    % ptmp = p0 + dt/2*F;
    % qtry = dKdp(ptmp)*dt/2 + q(:,i-1);

    [lambda,iszeros1] = fsolve(@(X) sys1(X,p0,q0,G(q0),S,dt,dTdq0,dKdp),[zeros(NC,1)],opt);
    if any(abs(iszeros1) > 1e-8)
        warning(strcat("CONSTRAINTS NOT RESPECTED, value of sum(zeros) in SYS1 is:",num2str(sum(abs(iszeros1)))," time =",num2str(t(i))))
    end
    
    Gq = G(q0)*lambda; Gq = reshape(Gq,[],3);
    F = Gq - dTdq0; 
    ptmp = p0 + dt/2*F;
    q(:,:,i) = q0 + dt*dKdp(ptmp);
 
    % @q
    %F =  G(q(:,i))*lambda - dTdq(q(:,i)); 
    %Force tmp
    Ftmp = - dTdq(q(:,:,i));

    %momentum
    p(:,:,i) = Ftmp*dt/2 + ptmp;

    % [p_mu,iszeros2] = fsolve(@(X) sys2(X,ptmp,G(q(:,i)),F,dt,dKdp,dGdt),[ptmp;0*lambda],opt);
    % if any(abs(iszeros2) > 1e-8)
    %     warning(strcat("CONSTRAINTS NOT RESPECTED, value of sum(zeros) in SYS2 is:",strcat(num2str(sum(abs(iszeros2))))))
    % end
    % 
    % p(:, i) = p_mu(1:N);
    %Frv = G(q(:,i))*mu + F;   
    %p(:,i) = ptmp + dt/2 * Frv;
end
end

function toZero = sys1(lambda,p0,q0,G0,S,dt,dTdq0,dKdp)
        Gq = G0*lambda; Gq = reshape(Gq,[],3);
        F = Gq - dTdq0; 
        ptmp = p0 + dt/2*F;
        q = q0 + dt*dKdp(ptmp);
 
        toZero = S(q); 
end

function toZero = sys2(p_mu,ptmp,Gq,F,dt,dKdp,dGdt)

        p = p_mu(1:numel(ptmp));
        mu = p_mu(numel(ptmp)+1:end);
        Frv = Gq*mu + F;   
        toZero = p - ptmp - dt/2 * Frv;
        toZero = [toZero; Gq'*dKdp(p)];
        % toZero = dGdt(Gq,dKdp(p));
end
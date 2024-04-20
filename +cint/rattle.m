function [q,p] = rattle(q0,p0,dTdq,dKdp,G,S,t)

NP = size(q0,1);
NT = numel(t);

q = zeros(NP,3,NT);
p = zeros(NP,3,NT);

q(:,:,1) = q0;
p(:,:,1) = p0;

NC = numel(S(q0));
opt = optimoptions("fsolve","Display","none","OptimalityTolerance",1e-40,...
        "FunctionTolerance",1e-40,"FiniteDifferenceType","central","StepTolerance",1e-25);

for i = 2:NT
    % time step
    dt = t(i) - t(i-1);
    
    p0 = p(:,:,i-1);
    q0 = q(:,:,i-1);
    dTdq0 = dTdq(q0);

    [lambda,iszeros1] = fsolve(@(X) sys1(X,p0,q0,G(q0),S,dt,dTdq0,dKdp),zeros(NC,1), opt);
    if any(abs(iszeros1) > 1e-8)
        warning(strcat("CONSTRAINTS NOT RESPECTED, value of sum(zeros) in SYS1 is:",num2str(sum(abs(iszeros1)))," time =",num2str(t(i))))
    end
    
    Gq = G(q0)*lambda; Gq = reshape(Gq,[],3);
    F = Gq - dTdq0; 
    ptmp = p0 + dt/2*F;
    q(:,:,i) = q0 + dt*dKdp(ptmp);
 
    % second part
    Ftmp = - dTdq(q(:,:,i));

    Gqp = G(q(:,:,i));
    Gqpx = Gqp(1:size(Gqp,1)/3,:)*(ptmp(:,1)*2/dt + Ftmp(:,1));
    Gqpy = Gqp(size(Gqp,1)/3+1:2*size(Gqp,1)/3,:)*(ptmp(:,2)*2/dt + Ftmp(:,2));
    Gqpz = Gqp(2*size(Gqp,1)/3+1:end,:)*(ptmp(:,3)*2/dt + Ftmp(:,3));
    b = Gqpx + Gqpy + Gqpz;
    A = Gqp'*Gqp;

    oneOnMu = A\b; mu = 1./oneOnMu;

    p(:,:,i) = ptmp + dt/2*Ftmp - dt/2*reshape(Gqp*mu,[],3);
end
end

function toZero = sys1(lambda,p0,q0,G0,S,dt,dTdq0,dKdp)
        Gq = G0*lambda; Gq = reshape(Gq,[],3);
        F = Gq - dTdq0; 
        ptmp = p0 + dt/2*F;
        q = q0 + dt*dKdp(ptmp);
 
        toZero = S(q); 
end
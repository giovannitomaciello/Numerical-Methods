function [q,p] = rattle(q0,p0,dTdq,dKdp,G,S,t)

[NP, ND] = size(q0);
NT = numel(t);

q = zeros(NP,ND,NT);
p = zeros(NP,ND,NT);

q(:,:,1) = q0;
p(:,:,1) = p0;

S0 = S(q0); NC = numel(S0);
unk1 = zeros(NC,1); unk2 = unk1;

opt = optimoptions("fsolve","Display","none","OptimalityTolerance",1e-20,...
        "FunctionTolerance",1e-20,"FiniteDifferenceType","central","StepTolerance",1e-18);

for i = 2:NT
    % time step
    dt = t(i) - t(i-1);
    
    p0 = p(:,:,i-1);
    q0 = q(:,:,i-1);
    [lambda,iszeros1] = fsolve(@(X) sys1(X,p0,q0,G(q0),S,dt,dTdq(q0),dKdp),unk1,opt);
    if any(abs(iszeros1) > 1e-8)
        warning(strcat("CONSTRAINTS NOT RESPECTED, value of sum(zeros) in SYS1 is:",strcat(num2str(sum(abs(iszeros1))))))
    end

    % @q0
    Gq0 = reshape(G(q0)*lambda,[],3);
    F = Gq0 - dTdq(q0); 
    ptmp = p0 + dt/2*F;
    q(:,:,i) = q0 + dt*dKdp(ptmp);
 
    % @q
    F = Gq0 - dTdq(q(:,:,i)); 

    [mu,iszeros2] = fsolve(@(X) sys2(X,ptmp,G(q(:,:,i)),F,dt,dKdp),unk2,opt);
    if any(abs(iszeros2) > 1e-8)
        warning(strcat("CONSTRAINTS NOT RESPECTED, value of sum(zeros) in SYS2 is:",strcat(num2str(sum(abs(iszeros2))))))
    end

    Frv = reshape(G(q(:,:,i))*mu,[],3) + F;   
    p(:,:,i) = ptmp + dt/2 * Frv;
end
end

function toZero = sys1(lambda,p0,q0,Gq0,S,dt,dTdq0,dKdp)
        Gq = reshape(Gq0*lambda,[],3);
        F = Gq - dTdq0; 
        ptmp = p0 + dt/2*F;
        q = q0 + dt*dKdp(ptmp);

        toZero = S(q);
end

function toZero = sys2(mu,ptmp,Gq,F,dt,dKdp)
        Frv = reshape(Gq*mu,[],3) + F;   
        p = ptmp + dt/2 * Frv;

        
               
        toZero = sum(sum(reshape(Gq,[],size(p,1),3).*reshape(dKdp(p),size(p,1),1,3) - reshape(dKdp(p)',1,size(p,1),3),2),3);
end
function [qx qy px py Energy]=crankNick(q_x_0,q_y_0,p_x_0,p_y_0,T,K,dTdqx,dTdqy,dK_dp,t,NP,dsdt,s0,theta)

qx=zeros(NP,length(t));
qy=zeros(NP,length(t));

px=zeros(NP,length(t));
py=zeros(NP,length(t));

s=zeros(NP,length(t));

qx(:,1)=q_x_0';
qy(:,1)=q_y_0';

px(:,1)=p_x_0';
py(:,1)=p_y_0';

Energy= zeros(numel(t),1);
Energy(1)=T(qx(:,1),qy(:,1)) + K(px(:,1),py(:,1))+theta*sum(s0);

for i=2:length(t)
    % time step
    dt = t(i) - t(i-1);
    % init val
    qx0 = qx(:,i-1);
    qy0 = qy(:,i-1);
    px0 = px(:,i-1);
    py0 = py(:,i-1);
    
    opt = optimoptions("fsolve","Display","none");
    unk = fsolve(@(unk) sysEB(unk, qx0, qy0, px0, py0, dt, dK_dp, dTdqx, dTdqy),....
        [px0; py0; qx0; qy0], opt);

    px(:,i) = unk(1:NP);
    py(:,i) = unk(NP+1:2*NP);
    qx(:,i) = unk(2*NP+1:3*NP);
    qy(:,i) = unk(3*NP+1:4*NP);

    s(:,i)=(dsdt(px(:,i),py(:,i)) + dsdt(px(:,i-1),py(:,i-1)))*dt/2  +s(:,i-1);

    Energy(i) = T(qx(:,i),qy(:,i)) + K(px(:,i),py(:,i))+theta*sum(s(:,i));
end
end

function toZero = sysEB(unk, qx0, qy0, px0, py0, dt, dK_dp, dTdqx, dTdqy)

    NP = length(unk)/4;
    px = unk(1:NP);
    py = unk(NP+1:2*NP);
    qx = unk(2*NP+1:3*NP);
    qy = unk(3*NP+1:4*NP);

    % position
    zeroqx = (dK_dp(px) + dK_dp(px0))*dt/2 - qx + qx0;
    zeroqy = (dK_dp(py) + dK_dp(py0))*dt/2 - qy + qy0;

    %momentum
    zeropx = (dTdqx(qx,qy,px) + dTdqx(qx0,qy0,px0))*dt/2 + px - px0;
    zeropy = (dTdqy(qx,qy,py) + dTdqy(qx0,qy0,py0))*dt/2 + py - py0;

    toZero = [zeroqx; zeroqy; zeropx; zeropy];
end
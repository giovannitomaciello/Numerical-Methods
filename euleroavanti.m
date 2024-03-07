function [qx qy px py Energy]=euleroavanti(q_x_0,q_y_0,p_x_0,p_y_0,T,K,dTdqx,dTdqy,dK_dp,t,NP,dsdt,s0,theta)

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
    
    % position
    qx(:,i)=dK_dp(px(:,i-1))*dt+qx(:,i-1);
    qy(:,i)=dK_dp(py(:,i-1))*dt+qy(:,i-1);

    %momentum
    px(:,i) = - dTdqx(qx(:,i-1),qy(:,i-1),px(:,i-1))*dt + px(:,i-1);
    py(:,i) = - dTdqy(qx(:,i-1),qy(:,i-1),py(:,i-1))*dt + py(:,i-1);

    s(:,i) = dsdt(px(:,i-1),py(:,i-1))*dt + s(:,i-1);

    % Energy
    Energy(i) = T(qx(:,i),qy(:,i)) + K(px(:,i),py(:,i))+theta*sum(s(:,i));

end
end
function [q,p]=velVerlet(q0,p0,dTdq,dKdp,t)

[NP, ND] = size(q0);

NT=numel(t);

q=zeros(NP,ND,NT);
p=zeros(NP,ND,NT);

q(:,:,1)=q0;
p(:,:,1)=p0;

%Force tmp
 
%Ftmp ha dimensioni [NP,ND] e Ftmp(i,j) rappresenta la componente j sulla particella i  
Ftmp = - dTdq(q0);

for i=2:NT
    % time step
    dt = t(i) - t(i-1);

    %momentum tmp
    ptmp = Ftmp*dt/2 + p(:,:,i-1);
   
    % position
    q(:,:,i) = dKdp(ptmp)*dt + q(:,:,i-1);
    

    %Force tmp
    Ftmp = - dTdq(q(:,:,i));
   
    %momentum
    p(:,:,i) = Ftmp*dt/2 + ptmp;
    
end
end
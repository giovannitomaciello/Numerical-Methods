function [q,p] = velVerlet(q0,p0,dTdq,dKdp,t)

[NP, ND] = size(q0);

NT = numel(t);

q = zeros(NP,ND,NT);
p = zeros(NP,ND,NT);

q(:,:,1) = q0;
p(:,:,1) = p0;

%Force tmp
 
%Ftmp ha dimensioni [NP,ND] e Ftmp(i,j) rappresenta la componente j sulla particella i  
Ftmp = - dTdq(q0);
dt = diff(t);

for i = 2:NT
    %momentum tmp
    ptmp = Ftmp*dt(i-1)/2 + p(:,:,i-1);
   
    % position
    q(:,:,i) = dKdp(ptmp)*dt(i-1) + q(:,:,i-1);
    

    %Force tmp
    Ftmp = - dTdq(q(:,:,i));
   
    %momentum
    p(:,:,i) = Ftmp*dt(i-1)/2 + ptmp;
end
end
clear 
clc
close all


m = [2];
k = [10];
l = [1];
NP = length(m);
ND = 3;
g = [0 -9.81 0]';

K = @(p)  sum(vecnorm(p').^2/2./m);
T = @(q) -m'*q*g'+k/2*((vecnorm(diff([zeros(1,ND);q])')-l).^2)';

dKdp = @(p) p(:)./m(:);
dTdq = @(q) - m*g(:)+k*q(:);
G = @(q) constraints(q);
dGdt = @(G,dKdp) derivative_constraints(G,dKdp);
S = @(q) Sfunc(q);
%init
q0 = [-1 0 0];
p0 = [0 0 0];
t = 0:0.01:10;

%!!! CHOOSE A METHOD !!!%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%[q p] = int.velVerlet(q0,p0,dTdq,dKdp,t);
%[q p] = int.euleroindietro(q0,p0,dTdq,dKdp,t);
%[q p] = int.euleroavanti(q0,p0,dTdq,dKdp,t);
%[q p] = int.crankNick(q0,p0,dTdq,dKdp,t);
% [q,p] = int.posVerlet(q0,p0,dTdq,dKdp,m,t);
%[q,p] = int.symplecticEuler(q0,p0,dTdq,dKdp,t);
[q p] = cint.rattle(q0,p0,dTdq,dKdp,G,S,t);
q = reshape(q,1,3,[]);
p = reshape(p,1,3,[]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% for i=1:length(t)
%     Energy(i) = T(q(:,:,i)) + K(p(:,:,i));
% end

%%
for i = 1:length(t)
    % plot
    hold on
    plot(q(:,1,i),q(:,2,i),'ro','MarkerSize',10,'MarkerFaceColor','r')
    plot([0 q(:,1,i)],[0 q(:,2,i)],'r-');
    xlim([-5 5])
    ylim([-5 5])
    drawnow
    clf
end
%%
%plot(t,Energy)
%%
for i = 1:length(t)
    qxx(i) = q(1,1,i);
end
figure(2)
plot(t,qxx), grid, title('Position Verlet'), xlabel('time t'), ylabel('q_x')
ylim([-1 1])


function G = constraints(q)
    q =reshape(q,[],3);
    n = size(q,1);
 
    G = zeros(n,3);
   
    % sum
    G(:,1) = 2*q(:,1);
    G(:,2) = 2*q(:,2);
    G(:,3) = 2*q(:,3);

    G = G(:);
end
function dGdt = derivative_constraints(G,dKdp)
    G = reshape(G,[],3);
    dKdp = reshape(dKdp,[],3);

    n = size(dKdp,1);
    dGdt = zeros(n,3,n);
    dx = dKdp(:,1) - dKdp(:,1)'; 
    dy = dKdp(:,2) - dKdp(:,2)'; 
    dz = dKdp(:,3) - dKdp(:,3)'; 

    dGdt(:,1,:) = dx;
    dGdt(:,2,:) = dy;
    dGdt(:,3,:) = dz;

    dGdt = reshape(dGdt,n*3,n);
    dGdt = dGdt.*G';
    dGdt = dGdt(:);
end
function out = Sfunc(q)
    out = sum(q.^2)-1;
    %out(1:size(q,1)) = zeros();

end
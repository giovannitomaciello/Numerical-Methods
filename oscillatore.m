clear 
clc
close all

int = INT;

m = [2];
k = [10];
l = [0];
NP = length(m);
ND = 3;
g = [0 -9.81 0].*ones(NP,ND);

K = @(p)  sum(vecnorm(p').^2/2./m);
T = @(q) k/2*((vecnorm(diff([zeros(1,ND);q])')-l).^2)';

dKdp = @(p) p./m';
dTdq = @(q) k*q;

%init
q0 = [-1 1 0];
p0 = [0 0 0];
t = 0:0.01:10;

%!!! CHOOSE A METHOD !!!%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%[q p] = int.velVerlet(q0,p0,dTdq,dKdp,t);
%[q p] = int.euleroindietro(q0,p0,dTdq,dKdp,t);
%[q p] = int.euleroavanti(q0,p0,dTdq,dKdp,t);
%[q p] = int.crankNick(q0,p0,dTdq,dKdp,t);
% [q,p] = int.posVerlet(q0,p0,dTdq,dKdp,m,t);
[q,p] = int.symplecticEuler(q0,p0,dTdq,dKdp,t);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:length(t)
    Energy(i) = T(q(:,:,i)) + K(p(:,:,i));
end

%%
% for i = 1:length(t)
%     % plot
%     hold on
%     plot(q(:,1,i),q(:,2,i),'ro','MarkerSize',10,'MarkerFaceColor','r')
%     plot([0 q(:,1,i)],[0 q(:,2,i)],'r-');
%     xlim([-5 5])
%     ylim([-5 5])
%     drawnow
%     clf
% end
%%
plot(t,Energy)
%%
for i = 1:length(t)
    qxx(i) = q(1,1,i);
end
figure(2)
plot(t,qxx), grid, title('Position Verlet'), xlabel('time t'), ylabel('q_x')
ylim([-1 1])
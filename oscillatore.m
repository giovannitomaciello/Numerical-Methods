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

[q p] = int.velVerlet(q0,p0,dTdq,dKdp,t);
%[q p] = int.euleroindietro(q0,p0,dTdq,dKdp,t);
%[q p] = int.euleroavanti(q0,p0,dTdq,dKdp,t);
%[q p] = int.crankNick(q0,p0,dTdq,dKdp,t);

for i=1:length(t)
 Energy(i) = T(q(:,:,i)) + K(p(:,:,i));
end

%%
for i = 1:length(t)
    % plot
    hold on
    plot(q(:,1,i),q(:,2,i),'ro','MarkerSize',10,'MarkerFaceColor','r')
    plot([0 q(:,1,i)],[0 q(:,2,i)],'r-');
    % plot([qx(1,i) qx(2,i)],[qy(1,i) qy(2,i)],'g-');
    % plot(qx(2,i),qy(2,i),'go','MarkerSize',5,'MarkerFaceColor','g')
    % plot([qx(2,i) qx(3,i)],[qy(2,i) qy(3,i)],'b-');
    % plot(qx(3,i),qy(3,i),'go','MarkerSize',5,'MarkerFaceColor','b')
    % text(0,5,"Timer: "+num2str(t(i),2))

    xlim([-5 5])
    ylim([-5 5])
    drawnow
%     pause(.1)
    clf
end
%%
plot(t,Energy)

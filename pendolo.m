clear 
clc
close all

int = INT;

m = [1 .2 .1];
NP = length(m);
k = [100 100 100];
l = [1 1 1];
ND = 3;
g = [0 -9.81 0];


K = @(p)  sum(vecnorm(p').^2/2./m);
T = @(q) -m*q*g' + k/2*((vecnorm(diff([zeros(1,ND);q])')-l).^2)';

dKdp = @(p) p./m';
v = ones(1,NP);
matrice = diag(v) - diag(v(2:end),1);
dTdqi = @(q,i) - m'*g(i) - matrice*(k.*(diff([0;q(:,i)]))'.*(l./vecnorm(diff([0 0 0;q])')-1))';
dTdq = @(q) [dTdqi(q,1) dTdqi(q,2) dTdqi(q,3)];


%init
q0 = [1 0 0; 2 0 0; 3 0 0];
p0 = [0 0 0;0 0 0; 0 0 0];
t = 0:0.005:10;

%[q p] = int.velVerlet(q0,p0,dTdq,dKdp,t);
%[q p] = int.euleroindietro(q0,p0,dTdq,dKdp,t);
%[q p] = int.euleroavanti(q0,p0,dTdq,dKdp,t);
%[q p] = int.crankNick(q0,p0,dTdq,dKdp,t);

for i = 1:length(t)
    Tpl(i) = T(q(:,:,i));
    Kpl(i) = K(p(:,:,i));
    Energy(i) = T(q(:,:,i)) + K(p(:,:,i));
end

%%
for i = 1:length(t)
    % plot
    
    plot([0;q(:,1,i)],[0;q(:,2,i)],'r-o','MarkerSize',10,'MarkerFaceColor','r')

    %plot([0 q(1,:,i)],[0 q(2,:,i)],'r-');
    % plot([qx(1,i) qx(2,i)],[qy(1,i) qy(2,i)],'g-');
    % plot(qx(2,i),qy(2,i),'go','MarkerSize',5,'MarkerFaceColor','g')
    % plot([qx(2,i) qx(3,i)],[qy(2,i) qy(3,i)],'b-');
    % plot(qx(3,i),qy(3,i),'go','MarkerSize',5,'MarkerFaceColor','b')
    % text(0,5,"Timer: "+num2str(t(i),2))

    xlim([-5 5])
    ylim([-5 5])
    drawnow
    hold off
end
%%

figure
hold on
plot(t,Energy)
plot(t,Kpl)
plot(t,Tpl)
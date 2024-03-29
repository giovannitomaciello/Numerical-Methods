clear 
clc
close all

int = INT;
cint = CINT;

m = [1 .2]';% .1 .05 0.025 0.0125]';
NP = length(m);
k = [100 100];% 100 100 100 100];
l = [1 1];% 1 1 0.5 0.5];
ND = 3;
g = [0 -9.81 0];


K = @(p)  sum(vecnorm(p').^2/2./m');
T = @(q) -m'*q*g' + k/2*((vecnorm(diff([zeros(1,ND);q])')-l).^2)';

dKdp = @(p) p./m;
v = ones(1,NP);
matrice = diag(v) - diag(v(2:end),1);
dTdqi = @(q,i) - m*g(i) - matrice*(k.*(diff([0;q(:,i)]))'.*(l./vecnorm(diff([zeros(1,ND);q])')-1))';
dTdq = @(q) [dTdqi(q,1) dTdqi(q,2) dTdqi(q,3)];
G = @(q,lambda) constraints(q,lambda);

%init
q0 = [1 0 0; 2 0 0];%; 3 0 0; 4 0 0; 4.5 0 0; 5 0 0];
p0 = [0 0 0; 0 0 0];%; 0 0 0; 0 0 0; 0 0 0; 0 0 0];
t = 0:0.001:1;
C = sparse(NP);
C(1,2) = 1; C(2,1) = 1;

%!!! CHOOSE A METHOD !!!%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%[q p] = int.velVerlet(q0,p0,dTdq,dKdp,t);
%[q p] = int.euleroindietro(q0,p0,dTdq,dKdp,t);
%[q p] = int.euleroavanti(q0,p0,dTdq,dKdp,t);
%[q p] = int.crankNick(q0,p0,dTdq,dKdp,t);
%[q p] = int.posVerlet(q0,p0,dTdq,dKdp,m,t);
[q, p] = int.symplecticEuler(q0,p0,dTdq,dKdp,t);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% for i = 1:length(t)
%     Tpl(i) = T(q(:,:,i));
%     Kpl(i) = K(p(:,:,i));
%     Energy(i) = T(q(:,:,i)) + K(p(:,:,i));
% end

%%

% figure
% hold on
% %plot(t,Energy - mean(Energy),"LineWidth",1.4)
% %plot(t,Kpl,"LineWidth",1.4)
% %plot(t,Tpl,"LineWidth",1.4)
% legend(["Total Energy","Kin Energy", "Pot Energy"])
% xlabel("Time [s]","FontSize",11,"FontWeight","bold")
% ylabel("Energy [J]","FontSize",11,"FontWeight","bold")
% %ylim([-22 31])

%%
figure
for i = 1:1:length(t)
    % plot
    plot([0;q(:,1,i)],[0;q(:,2,i)],'r-o','MarkerSize',10,'MarkerFaceColor','r')
    xlim([-5.5 5.5])
    ylim([-5.5 5.5])
    dist = vecnorm(q(1,:,i) - q(2,:,i)) - 1;
    % annotate time
    text(0.85,0.95,sprintf("t = %.2f s",t(i)),'Units','normalized')
    drawnow
    hold off
end


function [G,r2] = constraints(q,lambda)
    n = size(q,1);
    G = zeros(n,3);
    dx = q(:,1) - q(:,1)'; dxC = 2*dx.*lambda;
    dy = q(:,2) - q(:,2)'; dyC = 2*dy.*lambda;
    dz = q(:,3) - q(:,3)'; dzC = 2*dz.*lambda;
    r2 = dx.^2 + dy.^2 + dz.^2;

    % sum
    G(:,1) = sum(dxC,2);
    G(:,2) = sum(dyC,2);
    G(:,3) = sum(dzC,2);
end

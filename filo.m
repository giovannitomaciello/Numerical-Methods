clear 
clc
close all

int = INT;
cint = CINT;

dx = 0.125;
L = 5;
m = 1*ones(L/dx,1)/(L/dx);
NP = length(m);
k = [10000, ones(1,L/dx-1)];
l = dx*ones(1,L/dx);
ND = 3;
g = [0 -9.81 0];

K = @(p)  sum(vecnorm(p').^2/2./m');
T = @(q) -m'*q*g' + k/2*((vecnorm(diff([zeros(1,ND);q])')-l).^2)';

dKdp = @(p) p./m;
v = ones(1,NP);
matrice = diag(v) - diag(v(2:end),1);
dTdqi = @(q, i, constraintsEqZero) - m*g(i) - constraintsEqZero.*matrice*(k.*(diff([0;q(:,i)]))'...
    .*(l./vecnorm(diff([zeros(1,ND);q])')-1))';
dTdq = @(q, constraintsEqZero) [dTdqi(q,1, constraintsEqZero) dTdqi(q,2, constraintsEqZero)...
    dTdqi(q,3, constraintsEqZero)];
G = @(q,lambda) constraints(q,lambda);

%init
q0 = [linspace(dx,L,L/dx)', zeros(NP,2)];
p0 = zeros(NP,3);
t = 0:0.001:5;

%% C
C = sparse(size(q0,1),size(q0,1));
for i = 2:size(q0,1)
    C(i-1,i) = dx;
    C(i,i-1) = dx;
end

%!!! CHOOSE A METHOD !!!%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
q = cint.shake(q0,p0,dTdq,dKdp,G,C,m,t);
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
for i = 1:1:length(t)
    dist(i) = sum(vecnorm(diff(q(:,:,i))') - dx)/dx;
end

figure
plot(t,dist,"LineWidth",2)
title("relative cumulative error constraints")
xlabel("time [s]")
xlabel("time [s]")

%%
figure
for i = 1:1:length(t)
    % plot
    plot([0;q(:,1,i)],[0;q(:,2,i)])
    xlim([-5.5 5.5])
    ylim([-5.5 5.5])
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

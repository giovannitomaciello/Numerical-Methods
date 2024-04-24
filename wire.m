clear 
clc
close all

dx = .2; L2 = dx^2;
L = 2;
m = 1*ones(L/dx,1)/(L/dx);
NP = length(m);
k = ones(1,L/dx);
l = dx*ones(1,L/dx);
ND = 3;
g = [0 -9.81 0];

K = @(p)  sum(vecnorm(p').^2/2./m');
T = @(q) -m'*q*g';
S = @(q) Sfunc(q,L2);

dKdp = @(p) p./m;
v = ones(1,NP-1);
dTdq = @(q) - m*g;
G = @(q) constraints(q);
dGdt = @(G,dKdp) derivative_constraints(G,dKdp);

%init
q0 = [linspace(dx,L,L/dx)' zeros(NP,2)];
p0 = zeros(NP,3);
t = 0:0.001:5;


%!!! CHOOSE A METHOD !!!%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[q p] = cint.shake(q0,p0,dTdq,dKdp,G,S,t,m);
%[q p] = cint.rattle(q0,p0,dTdq,dKdp,G,S,t);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
for i = 1:length(t)
    Tpl(i) = T(q(:,:,i));
    Kpl(i) = K(p(:,:,i));
    Energy(i) = Tpl(i) + Kpl(i);
end

%%
figure(1)
minq = min(q(:));
maxq = max(q(:));

for i = 1:5:length(t)
    % plot
    plot([0; q(:,1,i)],[0; q(:,2,i)])
    xlim([minq maxq])
    ylim([minq maxq])
    % annotate time
    text(0.85,0.95,sprintf("t = %.2f s",t(i)),'Units','normalized')
    drawnow
    hold off
end

%%
figure(2)
hold on
plot(t,Energy,"LineWidth",2)
plot(t,Kpl,"LineWidth",2)
plot(t,Tpl,"LineWidth",2)
legend(["Total Energy","Kin Energy", "Pot Energy"])
xlabel("Time [s]","FontSize",11,"FontWeight","bold")
ylabel("Energy [J]","FontSize",11,"FontWeight","bold")
ylim([-22 31])

%%
for i = 1:1:length(t)
    dist(i) = abs(sum(vecnorm(diff([0 0 0; q(:,:,i)])') - dx)/dx)*100;
end

figure(3)
plot(t,dist,"LineWidth",2)
title("sum relative error constraints [%]")
xlabel("time [s]")
xlabel("time [s]")


function G = constraints(q)
    NP = size(q,1);
    NC = NP;
    ND = 3;
    G = zeros(NP,NP,3);

    dqDiag = diff([0 0 0; q]);
    dqDiagM1 = -diff(q);

    dx = diag(dqDiag(:,1),0) + diag(dqDiagM1(:,1),1);
    dy = diag(dqDiag(:,2),0) + diag(dqDiagM1(:,2),1);
    dz = diag(dqDiag(:,3),0) + diag(dqDiagM1(:,3),1);

    G(:,:,1) = dx;
    G(:,:,2) = dy;
    G(:,:,3) = dz;
end

function out = Sfunc(q,L2)
    out(1:numel(q)/3,1) = sum(diff([0 0 0; q]).^2,2)-L2;
end
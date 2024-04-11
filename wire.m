clear 
clc
close all

dx = 0.9; L2 = dx^2;
L = 0.9;
m = 1*ones(L/dx+1,1)/(L/dx+1);
NP = length(m);
k = ones(1,L/dx);
l = dx*ones(1,L/dx);
ND = 3;
g = [0 -9.81 0];

K = @(p)  sum(vecnorm(p').^2/2./m');
T = @(q) -m'*q*g' + k/2*((vecnorm(diff(q)')-l).^2)';
S = @(q) Sfunc(q,L2);

% to stabilize sys
 constraintsEqZero = ones(NP,NP); constraintsEqZero(2:NP+1:NP*NP) = 0;
 constraintsEqZero(NP+1:NP+1:NP*NP) = 0;

dKdp = @(p) p./[m;m;m];
v = ones(1,NP-1);
matrice = diag(v) - diag(v(2:end),1);
dTdqi = @(q, i) - m(2:end)*g(i);% - matrice*(k.*(diff(q((i-1)*length(q)/3+1:i/3*length(q))))'...
    %.*(l./vecnorm(diff([reshape(q,[],3)])')-1))';
dTdq = @(q) [0; dTdqi(q,1); -m(1)*g(2); dTdqi(q,2); 0; dTdqi(q,3)];
G = @(q) constraints(q);
dGdt = @(G,dKdp) derivative_constraints(G,dKdp);

%init
q0 = [linspace(0,L,L/dx+1)'; zeros(NP*2,1)];
p0 = zeros(NP*3,1);
t = 0:0.0001:1;


%!!! CHOOSE A METHOD !!!%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%[q p] = cint.shake(q0,p0,dTdq,dKdp,G,S,t,m);
[q p] = cint.rattle(q0,p0,dTdq,dKdp,G,S,t,dGdt);
q = reshape(q,NP,3,[]);
p = reshape(p,NP,3,[]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
for i = 1:length(t)
    Tpl(i) = T(q(:,:,i));
    Kpl(i) = K(p(:,:,i));
    Energy(i) = Tpl(i) + Kpl(i);
end

%%

figure
hold on
plot(t,Energy,"LineWidth",1.4)
plot(t,Kpl,"LineWidth",1.4)
plot(t,Tpl,"LineWidth",1.4)
legend(["Total Energy","Kin Energy", "Pot Energy"])
xlabel("Time [s]","FontSize",11,"FontWeight","bold")
ylabel("Energy [J]","FontSize",11,"FontWeight","bold")
ylim([-22 31])

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
for i = 1:10:length(t)
    % plot
    plot([0;q(:,1,i)],[0;q(:,2,i)])
    xlim([-5.5 5.5])
    ylim([-5.5 5.5])
    % annotate time
    text(0.85,0.95,sprintf("t = %.2f s",t(i)),'Units','normalized')
    drawnow
    hold off
end

function G = constraints(q)
    q = reshape(q,[],3);
    % n = size(q,1);
    % G = zeros(n,3,n);
    % dx = q(:,1) - q(:,1)'; dxC = 2*dx;
    % dy = q(:,2) - q(:,2)'; dyC = 2*dy;
    % dz = q(:,3) - q(:,3)'; dzC = 2*dz;
    % 
    % % sum
    % G(:,1,:) = dxC;
    % G(:,2,:) = dyC;
    % G(:,3,:) = dzC;
    % 
    % G = reshape(G,n*3,n);

    A = [2*q(1,:); 0*q(1,:)];
    B = [2*(q(1,:) - q(2,:)); -2*(q(1,:) - q(2,:))];
    G = [A(:), B(:)];
end

function dGdt = derivative_constraints(G,dKdp)
    dKdp = reshape(dKdp,[],3);
    n = size(dKdp,1);
    %dGdt = zeros(n,3,n);
    %dx = dKdp(:,1) - dKdp(:,1)'; 
    %dy = dKdp(:,2) - dKdp(:,2)'; 
    %dz = dKdp(:,3) - dKdp(:,3)'; 

    %dGdt(:,1,:) = dx;
    %dGdt(:,2,:) = dy;
    %dGdt(:,3,:) = dz;

    %dGdt = reshape(dGdt,n*3,n);
    %dGdt = dGdt.*G;
    dGdt = G'*reshape(dKdp',[],1);
end

function out = Sfunc(q,L2)
    q = reshape(q,[],3);
    out(1) = sum(abs(q(1,:)));
    out(2:numel(q)/3,1) = sum(diff(q).^2,2)-L2;
end
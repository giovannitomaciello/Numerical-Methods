clc; clear all; close all
int = INT;
cint = CINT;

%% constants
kb = 1.380658e-23;
epsij = 119.8*kb/(1.66e-27); 
sigmaij = 0.341; % nm

%% icosahedron 12 vertici esterni, 13 considerando anche l'atomo centrale
phi = (1+sqrt(5))/2;            % sezione aurea
dist = sigmaij*2^(1/6)/2;       % distanza iniziale tra ogni atomo e quelli adiacenti 
qx0 = [0 0 0 0 1 1 -1 -1 phi phi -phi -phi 0]*dist; % nm
qy0 = [1 1 -1 -1 phi -phi phi -phi 0 0 0 0 0]*dist; % nm
qz0 = [phi -phi phi -phi 0 0 0 0 1 -1 1 -1 0]*dist; % nm

q0prima = [qx0', qy0', qz0']; % nm

% Rotazione per vedere da varie angolazioni la figura 3D
theta = pi/7; % = 20°
Ry = [cos(theta)  0  sin(theta);
          0       1       0    ;
     -sin(theta)  0  cos(theta)];
q0sim = q0prima*Ry;    % q0sim è perfettamente simmetrico

% 
Na = length(qx0);
r = 1;% + 0.01.*randn(Na,3);
q0 = q0sim.*r;

m = 66.34/1.66*ones(1,Na)'; % Da
px0 = zeros(1,Na); % nm/ns
py0 = zeros(1,Na); % nm/ns
pz0 = zeros(1,Na); % nm/ns
p0 = [px0', py0', pz0']; % nm/ns

%% Forces
dKdp = @(p) p./m;
F = @(q) LennardJonesForceConstrained(q, sigmaij, epsij);
G = @(q,lambda) constraints(q,lambda);

%% constraints
C = sparse(Na,Na);
C(1,2) = dist*phi*2;
C(2,1) = dist*phi*2;
C(3,4) = dist*phi*2;
C(4,3) = dist*phi*2;
C(4,5) = dist*phi*2;
C(5,4) = dist*phi*2;
ind_constraints = find(triu(C));

%% init
t = 0:1e-4:0.5;

%!!! CHOOSE A METHOD !!!%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%[q p] = int.velVerlet(q0,p0,F,dKdp,t);
%[q p] = int.euleroavanti(q0,p0,F,dKdp,t);
%[q p] = int.euleroindietro(q0,p0,F,dKdp,t);
%[q p] = int.crankNick(q0,p0,F,dKdp,t);
[q p] = rattle2(q0,p0,F,dKdp,G,C,m,t);
%q = cint.shake(q0,p0,F,dKdp,G,C,m,t);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% for i = 1:length(t)
%     T(i) = sum(m.*( sum( (p(:,:,i)./m).^2 ,2) )) /Na/2/kb;
%     E(i) = Energy(q(:,:,i), p(:,:,i), m, sigmaij, epsij);
% end
% 
% %%
% figure(1)
% plot(t/1e-9,T-mean(T),"LineWidth",1.4)
% xlabel("Time [ns]","FontSize",11,"FontWeight","bold")
% ylabel("T [K]","FontSize",11,"FontWeight","bold")
% 
% %%
% figure(2)
% plot(t/1e-9,(E-mean(E))/kb,"LineWidth",1.4)
% xlabel("Time [ns]","FontSize",11,"FontWeight","bold")
% ylabel("(E - E_0)/k_b [K]","FontSize",11,"FontWeight","bold")

%% 3D plot

%figure(4)
for i = 1:50:length(t)
    % plot
    %grid on
    %hold on
    DT = delaunayTriangulation(q(:,:,i));
    %tetramesh(DT,'faceAlpha',0,'lineStyle','--')
    %scatter3(DT.Points(:,1),DT.Points(:,2),DT.Points(:,3),80,'filled')
    %scatter3(DT.Points(1:2,1),DT.Points(1:2,2),DT.Points(1:2,3),80,'filled','r')
    dx = q(:,1,i) - q(:,1,i)';
    dy = q(:,2,i) - q(:,2,i)';
    dz = q(:,3,i) - q(:,3,i)';
    r=sqrt(dx.^2 + dy.^2 + dz.^2);
    error(i) = sum(r(ind_constraints) - C(ind_constraints));
    % annotate time
    %text(0.85,0.95,sprintf("t = %.4f ns",t(i)/1e-9),'Units','normalized')
    %xlim([min(q(:,1,:),[],'all') max(q(:,1,:),[],'all')])
    %ylim([min(q(:,2,:),[],'all') max(q(:,2,:),[],'all')])
    %view(3)
    %drawnow
    %clf
end

plot(t,error)
%% Functions

function [G,r2] = constraints(q,lambda)
    n = length(q);
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

function F = LennardJonesForceConstrained(q, sigmaij, epsij)
n = length(q);
F = zeros(n,3);

% calculate distance
dx = q(:,1) - q(:,1)';
dy = q(:,2) - q(:,2)';
dz = q(:,3) - q(:,3)';
r2 = dx.^2 + dy.^2 + dz.^2 + eye(n);

% calculate sigma2 - 6
sigma2 = (sigmaij^2)./r2;
sigma6 = sigma2.*sigma2.*sigma2;

% calculate force
Fmat = 48*epsij*sigma6.*(sigma6 - 0.5)./r2;
Fx = -Fmat.*dx;
Fy = -Fmat.*dy;
Fz = -Fmat.*dz;

% sum forces
F(:,1) = sum(Fx,2);
F(:,2) = sum(Fy,2);
F(:,3) = sum(Fz,2);
end

function E = Energy(q, p, m, sigmaij, epsij)
    n = length(q);

    % calculate kinetic energy
    K = sum(m.*( sum( (p./m).^2 ,2) )) /2;

    % calculate distance
    dx = q(:,1) - q(:,1)';
    dy = q(:,2) - q(:,2)';
    dz = q(:,3) - q(:,3)';
    r2 = dx.^2 + dy.^2 + dz.^2  + eye(n);

    % calculate sigma2 - 6
    sigma2 = (sigmaij^2)./r2;
    sigma6 = sigma2.*sigma2.*sigma2;

    %set diagonal to 0
    sigma6(1:n+1:end) = 0;

    % calculate potential energy
    U = 4*epsij*sigma6.*(sigma6 - 1);

    % sum potential energy
    E = sum(sum(U))/2 + K;
end
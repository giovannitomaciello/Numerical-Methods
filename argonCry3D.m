clc; clear all; close all

%% constants
kb = 1.380658e-23;
epsij = 119.8*kb/(1.66e-27); % J
sigmaij = 0.341; % m

%% icosahedron 12 vertici esterni, 13 considerando anche l'atomo centrale
phi = (1+sqrt(5))/2;            % sezione aurea
dist = sigmaij*2^(1/6)/2;       % distanza iniziale tra ogni atomo e quelli adiacenti 
qx0 = [0 0 0 0 1 1 -1 -1 phi phi -phi -phi 0]*dist; % nm
qy0 = [1 1 -1 -1 phi -phi phi -phi 0 0 0 0 0]*dist; % nm
qz0 = [phi -phi phi -phi 0 0 0 0 1 -1 1 -1 0]*dist; % nm
q0prima = [qx0', qy0', qz0']; % nm

% Rotazione per vedere da varie angolazioni la figura 3D
theta = pi/1; % 
Ry = [cos(theta)  0  sin(theta);
          0       1       0    ;
     -sin(theta)  0  cos(theta)];
q0sim = q0prima*Ry;    % q0sim è perfettamente simmetrico

% 
Na = length(qx0);
r = 1 + 0.01.*randn(Na,3);
q0 = q0sim.*r;

m = 66.34/1.66*ones(1,Na)'; % Da
% px0 = zeros(1,Na); % nm/ns
% py0 = zeros(1,Na); % nm/ns
% pz0 = zeros(1,Na); % nm/ns
% p0 = [px0', py0', pz0']; % nm/ns
pnn = randi([-5 5], Na-1,3);       
p_centro = -sum(pnn,1);    % così da avere p0 totale  nulla
p0 = [p_centro; pnn].*m*10; % nm/ns

%% Forces
dKdp = @(p) p./m;
F = @(q) LennardJonesForce(q, sigmaij, epsij);

%% init
t = 0:20e-6:0.5;

%!!! CHOOSE A METHOD !!!%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[q p] = int.velVerlet(q0,p0,F,dKdp,t);
%[q p] = int.euleroavanti(q0,p0,F,dKdp,t);
%[q p] = int.euleroindietro(q0,p0,F,dKdp,t);
%[q p] = int.crankNick(q0,p0,F,dKdp,t);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:length(t)
    T(i) = sum(m.*( sum( (p(:,:,i)./m).^2 ,2) )) /Na/2/kb*(1.66e-27);
    E(i) = Energy(q(:,:,i), p(:,:,i), m, sigmaij, epsij)*(1.66e-27);
end

%%
figure(1)
plot(t,T-mean(T),"LineWidth",1.4)
xlabel("Time [ns]","FontSize",11,"FontWeight","bold")
ylabel("T [K]","FontSize",11,"FontWeight","bold")

%%
figure(2)
plot(t,(E-mean(E))/kb,"LineWidth",1.4)
xlabel("Time [ns]","FontSize",11,"FontWeight","bold")
ylabel("(E - E_0)/k_b [K]","FontSize",11,"FontWeight","bold")

%% 3D plot

qc = q(end,:,:); % cosi l'atomo centrale sta fermo
q = q - qc;

figure(4)
for i = 1:length(t)
    % plot
    grid on
    hold on
    DT = delaunayTriangulation(q(:,:,i));
    tetramesh(DT,'faceAlpha',0,'lineStyle','--')
    scatter3(DT.Points(:,1),DT.Points(:,2),DT.Points(:,3),80,'filled')
    % annotate time
    text(0.85,0.95,sprintf("t = %.4f ns",t(i)),'Units','normalized')
    xlim([min(q(:,1,:),[],'all') max(q(:,1,:),[],'all')])
    ylim([min(q(:,2,:),[],'all') max(q(:,2,:),[],'all')])
    view(3)
    drawnow
    clf
end

%% Functions
function F = LennardJonesForce(q, sigmaij, epsij)
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
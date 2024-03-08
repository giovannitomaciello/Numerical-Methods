clc; clear all; close all

int = INT;

%% atoms pos, vel, m (m m/s kg)
qx0 = [0 0.02 0.34 0.36 -0.02 -0.35 -0.31]*1e-9; % nm
qy0 = [0 0.39 0.17 -0.21 -0.40 -0.16 0.21]*1e-9; % nm
qz0 = zeros(1,length(qx0)); % nm
q0 = [qx0', qy0', qz0']; % nm

m = 66.34e-27*ones(1,length(qy0)); % Kg

px0 = [-30 50 -70 90 80 -40 -80].*m; % nm/ns
py0 = [-20 -90 -60 40 90 100 -60].*m; % nm/ns
pz0 = zeros(1,length(qx0)); % nm/ns
p = [px0', py0', pz0']; % nm/ns

%% constants
kb = 1.380658e-23;
epsij = 119.8*kb; % J
sigmaij = 0.341*1e-9; % m

%% init E
Energy0 = -1260.2*kb; % J
T0 = 22.72; % K

%% Forces
dKdp = @(p) p./m';
F = @(q) LennardJonesForce(q, sigmaij, epsij);

%% init
t = 0:10e-15:0.2e-9;

[q, p] = int.velVerlet(q0,p,F,dKdp,t);
%[q, p] = int.crankNick(q0,p,F,dKdp,t);
%[q, p] = int.euleroindietro(q0,p,F,dKdp,t);
%[q, p] = int.euleroavanti(q0,p,F,dKdp,t);

for i = 1:length(t)
    T(i) = sum(m'.*( sum( (p(:,:,i)./m').^2 ,2) )) /7/2/kb;
    E(i) = Energy(q(:,:,i), p(:,:,i), m, sigmaij, epsij);
end

%%
figure
plot(t,T-T0)

%%
figure
plot(t,E-Energy0)

%%
for i = 1:length(t)
    % plot
    hold on
    scatter(q(:,1,i),q(:,2,i))

    drawnow
    pause(.1)
    clf    
end

function F = LennardJonesForce(q, sigmaij, epsij)
    n = length(q);
    F = zeros(n,3);

    % calculate distance
    dx = q(:,1) - q(:,1)';
    dy = q(:,2) - q(:,2)';
    dz = q(:,3) - q(:,3)';
    r = sqrt(dx.^2 + dy.^2 + dz.^2) + eye(n);

    % calculate sigma6 - 12
    sigma6 = (sigmaij./r).^6;
    sigma12 = sigma6.^2;

    % calculate force
    Fmat = 4*epsij*(12*sigma12 - 6*sigma6)./r;
    Fx = Fmat.*dx./r;
    Fy = Fmat.*dy./r;
    Fz = Fmat.*dz./r;

    % set diagonal to 0
    Fx(1:n+1:end) = 0;
    Fy(1:n+1:end) = 0;
    Fz(1:n+1:end) = 0;

    % sum forces
    F(:,1) = sum(Fx,2);
    F(:,2) = sum(Fy,2);
    F(:,3) = sum(Fz,2);
end

function E = Energy(q, p, m, sigmaij, epsij)
    n = length(q);
    E = 0;

    % calculate kinetic energy
    K = sum(m'.*( sum( (p./m').^2 ,2) )) /2;

    % calculate distance
    dx = q(:,1) - q(:,1)';
    dy = q(:,2) - q(:,2)';
    dz = q(:,3) - q(:,3)';
    r = sqrt(dx.^2 + dy.^2 + dz.^2)  + eye(n);

    % calculate sigma6 - 12
    sigma6 = (sigmaij./r).^6;
    sigma12 = sigma6.^2;

    %set diagonal to 0
    sigma6(1:n+1:end) = 0;

    % calculate potential energy
    U = 4*epsij*(sigma12 - sigma6);

    % sum potential energy
    E = sum(sum(U))/2 + K;
end



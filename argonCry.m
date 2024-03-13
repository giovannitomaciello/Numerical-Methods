clc; clear all; close all

int = INT;

%% atoms pos, vel, m (m m/s kg)
qx0 = [0 0.02 0.34 0.36 -0.02 -0.35 -0.31]*1e-9; % nm
qy0 = [0 0.39 0.17 -0.21 -0.40 -0.16 0.21]*1e-9; % nm
qz0 = zeros(1,length(qx0)); % nm
q0 = [qx0', qy0', qz0']; % nm

m = 66.34e-27*ones(1,length(qy0))'; % Kg

px0 = [-30 50 -70 90 80 -40 -80].*m'; % nm/ns
py0 = [-20 -90 -60 40 90 100 -60].*m'; % nm/ns
pz0 = zeros(1,length(qx0)); % nm/ns
p0 = [px0', py0', pz0']; % nm/ns

%% constants
kb = 1.380658e-23;
epsij = 119.8*kb; % J
sigmaij = 0.341*1e-9; % m

%% init E
Energy0 = -1260.2*kb; % J
T0 = 22.72; % K

%% Forces
dKdp = @(p) p./m;
F = @(q) LennardJonesForce(q, sigmaij, epsij);

%% init
t = 0:10e-15:0.2e-9;

%!!! CHOOSE A METHOD !!!%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[q p] = int.velVerlet(q0,p0,F,dKdp,t);
%[q p] = int.euleroavanti(q0,p0,F,dKdp,t);
%[q p] = int.crankNick(q0,p0,F,dKdp,t);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:length(t)
    T(i) = sum(m.*( sum( (p(:,:,i)./m).^2 ,2) )) /7/2/kb;
    E(i) = Energy(q(:,:,i), p(:,:,i), m, sigmaij, epsij);
end

%%
figure
plot(t/1e-9,T-T0,"LineWidth",1.4)
ylim([-30 30])
xlabel("Time [ns]","FontSize",11,"FontWeight","bold")
ylabel("T [K]","FontSize",11,"FontWeight","bold")

%%
figure
plot(t/1e-9,(E-Energy0)/kb,"LineWidth",1.4)
ylim([-30 30])
xlabel("Time [ns]","FontSize",11,"FontWeight","bold")
ylabel("(E - E_0)/k_b [K]","FontSize",11,"FontWeight","bold")

%%
figure
axis square

xlabel("x [m]","FontSize",11,"FontWeight","bold")
ylabel("y [m]","FontSize",11,"FontWeight","bold")

for i = 1:length(t)
    % plot
    scatter(q(:,1,i),q(:,2,i))
    % annotate time
    text(0.85,0.95,sprintf("t = %.4f ns",t(i)/1e-9),'Units','normalized')
    xlim([min(q(:,1,:),[],'all') max(q(:,1,:),[],'all')])
    ylim([min(q(:,2,:),[],'all') max(q(:,2,:),[],'all')])
    drawnow
    pause(.01) 
end

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



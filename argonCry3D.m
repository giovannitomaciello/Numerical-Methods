clc; clear all; close all
int = INT;

%% constants
kb = 1.380658e-23;
epsij = 119.8*kb; % J
sigmaij = 0.341*1e-9; % m

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
r = 1 + 0.01.*randn(Na,3);
q0 = q0sim.*r;

m = 66.34e-27*ones(1,Na)'; % Kg
px0 = zeros(1,Na); % nm/ns
py0 = zeros(1,Na); % nm/ns
pz0 = zeros(1,Na); % nm/ns
p0 = [px0', py0', pz0']; % nm/ns

%% Forces
dKdp = @(p) p./m;
F = @(q) LennardJonesForce(q, sigmaij, epsij);

%% init
t = 0:20e-15:0.5e-9;

%!!! CHOOSE A METHOD !!!%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%[q p] = int.velVerlet(q0,p0,F,dKdp,t);
%[q p] = int.euleroavanti(q0,p0,F,dKdp,t);
%[q p] = int.euleroindietro(q0,p0,F,dKdp,t);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:length(t)
    T(i) = sum(m.*( sum( (p(:,:,i)./m).^2 ,2) )) /Na/2/kb;
    E(i) = Energy(q(:,:,i), p(:,:,i), m, sigmaij, epsij);
end

%%
figure(1)
plot(t/1e-9,T-mean(T),"LineWidth",1.4)
xlabel("Time [ns]","FontSize",11,"FontWeight","bold")
ylabel("T [K]","FontSize",11,"FontWeight","bold")

%%
figure(2)
plot(t/1e-9,(E-mean(E))/kb,"LineWidth",1.4)
xlabel("Time [ns]","FontSize",11,"FontWeight","bold")
ylabel("(E - E_0)/k_b [K]","FontSize",11,"FontWeight","bold")

%% 3D plot
% Segmenti che uniscono gli atomi adiacenti
edges = [1, 3; 1, 5; 1, 7; 1, 9; 1, 11; 2, 4; 2, 5; 2, 7; 2, 10; 2, 12;
         3, 6; 3, 8; 3, 9; 3, 11; 4, 6; 4, 8; 4, 10; 4, 12; 5, 7; 5, 9;
         5, 10; 6, 8; 6, 9; 6, 10; 7, 11; 7, 12; 8, 11; 8, 12; 9,10; 11,12];

figure(4)
for i = 1:length(t)
    % plot
    grid on
    hold on
    for j = 1:30
        plot3(q(edges(j, :), 1, i), q(edges(j, :), 2, i), q(edges(j, :), 3, i), 'k-.', 'LineWidth', 1);
    end
    scatter3(q(:, 1, i), q(:, 2, i), q(:, 3, i), 100, 'red', 'filled');
    % annotate time
    text(0.85,0.95,sprintf("t = %.4f ns",t(i)/1e-9),'Units','normalized')
    xlim([min(q(:,1,:),[],'all') max(q(:,1,:),[],'all')])
    ylim([min(q(:,2,:),[],'all') max(q(:,2,:),[],'all')])
    view(3)
    drawnow
    pause(.01) 
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

    % set diagonal to 0
    Fx = Fx - diag(diag(Fx));
    Fy = Fy - diag(diag(Fy));
    Fz = Fz - diag(diag(Fz));

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
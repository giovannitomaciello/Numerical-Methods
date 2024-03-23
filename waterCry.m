clc; clear all; close all;

cint = CINT;
int = INT;

%% water crystal parameters
HOdist = 0.96; % oxygen atom
HOangle = 104.5; % hydrogen atoms
HHdist = 0.96*sind(HOangle/2)*2;
sigmaij = 2.725;
dist = sigmaij*2^(1/6)/2;
kb = 1.380658e-23;
epsij = 355.8060*kb/(1.66e-27); % J

%% 3D water molecules
positionNoRot = [0 HOdist*cosd(HOangle/2); -HOdist*sind(HOangle/2) 0; HOdist*sind(HOangle/2) 0];
shift = [0 HOdist*cosd(HOangle/2)];
positionNoRot = positionNoRot - shift;

%% rotate by 60 deg and add to the original position
position1layer = [];
for i = 1:6
    positionRot = positionNoRot*[cosd(120*i + 120) -sind(120*i + 120); sind(120*i + 120) cosd(120*i + 120)];
    positionTrasl = positionRot + [0 dist; 0 dist; 0 dist]*2;
    position1layer = [positionTrasl*[cosd(60*i) -sind(60*i); sind(60*i) cosd(60*i)]; position1layer];
end
position1layer = [position1layer , -ones(length(position1layer),1)*dist];

position2layer = [];
for i = 1:6
    positionRot = positionNoRot*[cosd(120*i - 120) -sind(120*i - 120); sind(120*i - 120) cosd(120*i - 120)];
    positionTrasl = positionRot + [0 dist; 0 dist; 0 dist]*2;
    position2layer = [positionTrasl*[cosd(60*i) -sind(60*i); sind(60*i) cosd(60*i)]; position2layer];
end
position2layer = [position2layer , zeros(length(position2layer),1)*dist];

position3layer = [];
for i = 1:6
    positionRot = positionNoRot*[cosd(120*i) -sind(120*i); sind(120*i) cosd(120*i)];
    positionTrasl = positionRot + [0 dist; 0 dist; 0 dist]*2;
    position3layer = [positionTrasl*[cosd(60*i) -sind(60*i); sind(60*i) cosd(60*i)]; position3layer];
end
position3layer = [position3layer , ones(length(position3layer),1)*dist];

q0 = [position1layer];% position2layer; position3layer];
p0 = zeros(size(q0));

%% constraint matrix
C = sparse(size(q0,1),size(q0,1));
connectivity = C;
for i = 3:3:size(q0,1)
    C(i-2,i-1) = HOdist;
    C(i-2,i) = HOdist;
    C(i-1,i) = HHdist;

    connectivity(i-2,i-1) = HOdist;
    connectivity(i-2,i) = HOdist;
end

C = C + C';

gplot3(connectivity,q0,'lineWidth',2,'Color','k')
hold on
scatter3(q0(1:3:end,1),q0(1:3:end,2),q0(1:3:end,3),80,"red","filled")
scatter3(q0(2:3:end,1),q0(2:3:end,2),q0(2:3:end,3),60,"blue","filled")
scatter3(q0(3:3:end,1),q0(3:3:end,2),q0(3:3:end,3),60,"blue","filled")

%%
Na = size(q0,1);
m = repmat([16 1 1],[1,Na/3])'; % Da
px0 = zeros(1,Na); % nm/ns
py0 = zeros(1,Na); % nm/ns
pz0 = zeros(1,Na); % nm/ns
p0 = [px0', py0', pz0']; % nm/ns

%% Forces
dKdp = @(p) p./m;
F = @(q, constraintsEqZero) LennardJonesForceConstrained(q, sigmaij, epsij, constraintsEqZero);
G = @(q,lambda) constraints(q,lambda);

%% init
t = 0:1e-8:0.0005;

%!!! CHOOSE A METHOD !!!%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
q = cint.shake(q0,p0,F,dKdp,G,C,m,t);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% plot
figure
for i = 1:size(q,1)
    gplot3(connectivity,q(:,:,i),'lineWidth',2,'Color','k')
    hold on
    scatter3(q(1:3:end,1,i),q(1:3:end,2,i),q(1:3:end,3,i),80,"red","filled")
    scatter3(q(2:3:end,1,i),q(2:3:end,2,i),q(2:3:end,3,i),60,"blue","filled")
    scatter3(q(3:3:end,1,i),q(3:3:end,2,i),q(3:3:end,3,i),60,"blue","filled")
    axis equal
    drawnow
    hold off
end
%%
for i = 1:length(t)
    T(i) = sum(m.*( sum( (p(:,:,i)./m).^2 ,2) )) /Na/2/kb*(1.66e-27);
end

%%
figure(1)
plot(t,T - mean(T),"LineWidth",1.4)
xlabel("Time [ns]","FontSize",11,"FontWeight","bold")
ylabel("T [K]","FontSize",11,"FontWeight","bold")

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

function F = LennardJonesForceConstrained(q, sigmaij, epsij, constraintsEqZero)
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
Fx = -Fmat.*dx.*constraintsEqZero;
Fy = -Fmat.*dy.*constraintsEqZero;
Fz = -Fmat.*dz.*constraintsEqZero;

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
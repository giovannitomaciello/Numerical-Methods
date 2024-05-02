clc; clear all; close all;


sigmaij = 2.725;
dist = 3.5;
kb = 1.380658e-23;
epsij = 355.8060*kb/(1.66e-27); % Da A^2/10^-10 s

%% water crystal parameters
HOdist = 0.96; % A
HOangle = 104.5; 
HHdist = 0.96*sind(HOangle/2)*2; %A
kb = 1.380658e-23; %J/K

sigma_OO = 2.75; %nm
epsi_OO = 80*kb/(1.66e-27); %

sigma_HH = sigma_OO; %nm
epsi_HH = epsi_OO ; %


sigma_HO = (sigma_OO + sigma_HH )/2 ; %A
epsi_HO = sqrt(epsi_HH*epsi_OO); %

% 
%dist = 2.725*2^(1/6)/2; %A

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

q0 = [position1layer; position2layer; position3layer];
p0 = zeros(size(q0));


%% constraint matrix
C = zeros(size(q0,1),size(q0,1));
connectivity = C;
for i = 3:3:size(q0,1)
    C(i-2,i-1) = HOdist;
    C(i-2,i) = HOdist;
    C(i-1,i) = HHdist;

    connectivity(i-2,i-1) = HOdist;
    connectivity(i-2,i) = HOdist;
end

C = C + C';
ind_constraints = find(triu(C));

% gplot3(connectivity,q0,'lineWidth',2,'Color','k')
% hold on
% scatter3(q0(1:3:end,1),q0(1:3:end,2),q0(1:3:end,3),80,"red","filled")
% scatter3(q0(2:3:end,1),q0(2:3:end,2),q0(2:3:end,3),60,"blue","filled")
% scatter3(q0(3:3:end,1),q0(3:3:end,2),q0(3:3:end,3),60,"blue","filled")

%%
Na = size(q0,1);
m = repmat([16 1 1],[1,Na/3])'; % Da
m_ind = m+m';

LJ_sigma = zeros(Na,Na);
LJ_sigma(m_ind==16*2) = sigma_OO;
LJ_sigma(m_ind==1*2) = sigma_HH;
LJ_sigma(m_ind==1+16) = sigma_HO;
stessa_molecola = C>0;
LJ_sigma(stessa_molecola) = 0;

LJ_epsi = zeros(Na,Na);
LJ_epsi(m_ind==16*2) = epsi_OO;
LJ_epsi(m_ind==1*2) = epsi_HH;
LJ_epsi(m_ind==1+16) = epsi_HO;
LJ_epsi(stessa_molecola) = 0;

r = LennardJonesForce(q0, LJ_sigma, LJ_epsi);

px0 = randn(1,Na); % A/ps
py0 = randn(1,Na); % A/ps
pz0 = randn(1,Na); % A/ps
p0 = [px0', py0', pz0']; % A/ps

%% Forces
dKdp = @(p) p./m;
F = @(q) LennardJonesForce(q, LJ_sigma, LJ_epsi); %Da A/ps^2
G = @(q) constraints(q,stessa_molecola);
S = @(q) Sfunc(q,C,ind_constraints);


r = G(q0);

%% init
t = 0:0.5e-5:0.001; %

%!!! CHOOSE A METHOD !!!%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[q p] = cint.shake(q0,p0,F,dKdp,G,S,t,m);
%[q p] = cint.rattle(q0,p0,F,dKdp,G,S,t);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% plot
figure
for i = 1:size(q,3)
    gplot3(connectivity,q(:,:,i),'lineWidth',2,'Color','k')
    hold on
    scatter3(q(1:3:end,1,i),q(1:3:end,2,i),q(1:3:end,3,i),80,"red","filled")
    scatter3(q(2:3:end,1,i),q(2:3:end,2,i),q(2:3:end,3,i),60,"blue","filled")
    scatter3(q(3:3:end,1,i),q(3:3:end,2,i),q(3:3:end,3,i),60,"blue","filled")
    axis equal
    drawnow

    hold off
end


%% Energy and error

E = zeros(1,numel(t));
U = zeros(1,numel(t));
K = zeros(1,numel(t));
error = zeros(1,numel(t));

for i = 1:length(t)
    [E(i), U(i), K(i)] = Energy(q(:,:,i), p(:,:,i), m, LJ_sigma, LJ_epsi);
    dx = q(:,1,i) - q(:,1,i)';
    dy = q(:,2,i) - q(:,2,i)';
    dz = q(:,3,i) - q(:,3,i)';
    r = (dx.^2 + dy.^2 + dz.^2).^(1/2);
    error_vet = abs(r(ind_constraints) - C(ind_constraints))./C(ind_constraints);
    error(i) = sum(error_vet,"all")*100;
end

figure(2)
plot(t,E,"LineWidth",2)
hold on
plot(t,K,"LineWidth",2)
plot(t,U,"LineWidth",2)
legend(["Total Energy","Kin Energy", "Pot Energy"])
xlabel("Time [s]","FontSize",11,"FontWeight","bold")
ylabel("Energy [J]","FontSize",11,"FontWeight","bold")


figure(3)
plot(t,error,"LineWidth",2)
xlabel("Time [s]","FontSize",11,"FontWeight","bold")
ylabel("Sum relative error constraints [%]","FontSize",11,"FontWeight","bold")


%%
% for i = 1:length(t)
%     T(i) = sum(m.*( sum( (p(:,:,i)./m).^2 ,2) )) /Na/2/kb*(1.66e-27);
% end

%%
% figure(1)
% plot(t,T - mean(T),"LineWidth",1.4)
% xlabel("Time [ns]","FontSize",11,"FontWeight","bold")
% ylabel("T [K]","FontSize",11,"FontWeight","bold")

%% Functions

function G = constraints(q,stessa_molecola)
    NP = size(q,1);
    NC = NP;
    ND = 3;
    G = zeros(NP,NP,3);
    dx = q(:,1) - q(:,1)'; dxC = 2*dx;
    dy = q(:,2) - q(:,2)'; dyC = 2*dy;
    dz = q(:,3) - q(:,3)'; dzC = 2*dz;

    Gx = zeros(NP,NP);
    Gy = Gx;
    Gz = Gx;

    Gx = stessa_molecola.*dxC;
    Gy = stessa_molecola.*dyC;
    Gz = stessa_molecola.*dzC;

    G(:,:,1) = Gx;
    G(:,:,2) = Gy;
    G(:,:,3) = Gz;

    
end

function out = Sfunc(q,C,ind_constraints)

 dx = q(:,1) - q(:,1)';
 dy = q(:,2) - q(:,2)';
 dz = q(:,3) - q(:,3)';
 r = dx.^2 + dy.^2 + dz.^2;

 out(1:numel(ind_constraints))= r(ind_constraints) - C(ind_constraints).^2;

end

function F = LennardJonesForce(q, LJ_sigma, LJ_epsi)
n = length(q);
F = zeros(n,3);

% calculate distance
dx = q(:,1) - q(:,1)';
dy = q(:,2) - q(:,2)';
dz = q(:,3) - q(:,3)';
r2 = dx.^2 + dy.^2 + dz.^2 + eye(n);

% calculate sigma2 - 6
sigma2 = (LJ_sigma.^2)./r2;
sigma6 = sigma2.*sigma2.*sigma2;

% calculate force
Fmat = 48*LJ_epsi.*sigma6.*(sigma6 - 0.5)./r2;
Fx = -Fmat.*dx;
Fy = -Fmat.*dy;
Fz = -Fmat.*dz;

% sum forces
F(:,1) = sum(Fx,2);
F(:,2) = sum(Fy,2);
F(:,3) = sum(Fz,2);
end

function [E,U,K] = Energy(q, p, m, LJ_sigma, LJ_epsi)
n = length(q);

% calculate kinetic energy
K = sum(m.*( sum( (p./m).^2 ,2) )) /2;

% calculate distance
dx = q(:,1) - q(:,1)';
dy = q(:,2) - q(:,2)';
dz = q(:,3) - q(:,3)';
r2 = dx.^2 + dy.^2 + dz.^2  + eye(n);

% calculate sigma2 - 6
sigma2 = (LJ_sigma.^2)./r2;
sigma6 = sigma2.*sigma2.*sigma2;

%set diagonal to 0
sigma6(1:n+1:end) = 0;

% calculate potential energy
U_mat = 4*LJ_epsi.*sigma6.*(sigma6 - 1);

U = sum(sum(U_mat))/2;

% sum potential energy
E = U + K;
end
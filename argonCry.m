clc; clear all; close all

int = INT;

%% atoms pos, vel, m (m m/s kg)
qx0 = [0 0.02 0.34 0.36 -0.02 -0.35 -0.31]*1e-9; 
qy0 = [0 0.39 0.17 -0.21 -0.40 -0.16 0.21]*1e-9;
qz0 = zeros(1,length(qx0));
q0 = [qx0', qy0', qz0'];

m = 66.34e-27*ones(1,length(qy0));

px0 = [-30 50 -70 90 80 -40 -80].*m;
py0 = [-20 -90 -60 40 90 100 -60].*m;
pz0 = zeros(1,length(qx0));
p = [px0', py0', pz0'];

%% constants
kb = 1.380658e-23;
epsij = 119.8*kb; % J
sigmaij = 0.341*1e-9; % m

%% init E
Energy0 = -1260.2*kb;

%% Forces
dKdp = @(p) p./m';
F = @(q) LennardJonesForce(q, sigmaij, epsij);

%% init
t = 0:10e-14:1e-11;

%[q, p] = int.velVerlet(q0,p,F,dKdp,t);
[q, p] = int.crankNick(q0,p,F,dKdp,t);
%[q, p] = int.euleroindietro(q0,p,F,dKdp,t);

for i = 1:length(t)
    T(i) = sum(m'.*( (sum(p(:,:,i)./m',2)).^2 )) /7/2/kb;
    E(i) = Energy(q(:,:,i), p(:,:,i), m, sigmaij, epsij);
end
%%
figure
plot(t,T-T(1))
figure
plot(t,E-E(1))

%%
for i = 1:length(t)
    % plot
    hold on
    scatter(q(:,1,i),q(:,2,i))

    drawnow
    pause(.1)
    clf    
end

function F = LennardJonesForce(positions, sigma, epsilon)

  % Calculate pairwise distances
  dx = bsxfun(@minus, positions(:, 1), permute(positions(:, 1), [2 1]));
  dy = bsxfun(@minus, positions(:, 2), permute(positions(:, 2), [2 1]));
  dz = bsxfun(@minus, positions(:, 3), permute(positions(:, 3), [2 1]));

  dx(isnan(dx)) = 0;
  dy(isnan(dy)) = 0;
  dy(isnan(dz)) = 0;

  r = sqrt(dx.^2 + dy.^2 + dz.^2);
  for i = 1:size(r,1)
      r(i,i) = 1e100;
  end

  sigmar6 = (sigma./r).^6;
  Fonr = 48*epsilon*(sigmar6).*(sigmar6 - 0.5);
  Fmat = Fonr.*r;

  % Assign forces in x and y directions
  F(:, 1) = sum(Fmat .* dx ./ r);
  F(:, 2) = sum(Fmat .* dy ./ r);
  F(:, 3) = sum(Fmat .* dz ./ r);
end

function out = Energy(positions, moments, m, sigma, epsilon)

  % Calculate pairwise distances
  dx = bsxfun(@minus, positions(:, 1), permute(positions(:, 1), [2 1]));
  dy = bsxfun(@minus, positions(:, 2), permute(positions(:, 2), [2 1]));
  dz = bsxfun(@minus, positions(:, 3), permute(positions(:, 3), [2 1]));
  dx(isnan(dx)) = 0;
  dy(isnan(dy)) = 0;
  dz(isnan(dy)) = 0;

  r = sqrt(dx.^2 + dy.^2 + dz.^2);
  for i = 1:size(r,1)
      r(i,i) = 1e100;
  end

  sigmar6 = (sigma./r).^6;
  U = 4*epsilon*(sigmar6.^2 -sigmar6);
  Epot = sum(triu(U).*r,"all");

  Ekin = 1/2*sum(( sum(moments.^2, 2) )./m');

  out = Ekin + Epot;

end



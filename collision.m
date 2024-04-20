clc; clear all; close all

%% distribute particles in a 2D domain as two colliding circles

% number of particles
N = 1000;

% radius of the circles
r = 20;

% center of the circles
x = [0 0];
y = [-50 50];

% generate random points in the circles
theta1 = 2*pi*rand(N/2,1);
theta2 = 2*pi*rand(N/2,1);
r1 = r*sqrt(rand(N/2,1));
r2 = r*sqrt(rand(N/2,1));
x1 = x(1) + r1.*cos(theta1);
y1 = y(1) + r1.*sin(theta1);
x2 = x(2) + r2.*cos(theta2);
y2 = y(2) + r2.*sin(theta2);

%% parameters of the simulation
L1 = 250; % length of the domain
L2 = 200; % width of the domain
epsilon = 5;
sigma = 1;
rCut = 2.5*sigma;
m = 1;
dt = 0.001;

% number of time steps
tFinal = 10;
nTime = round(tFinal/dt);

% create the meshgrid for the domain
grd.ncy = L2/rCut;
grd.ncx = L1/rCut;
grd.x = linspace (0, L1, grd.ncx+1);
grd.y = linspace (0, L2, grd.ncy+1);

% initialize the particles
ptcls.x = [[x1(:);x2(:)], [y1(:);y2(:)]]';
ptcls.v = 0*ptcls.x + 0.01*randn(size(ptcls.x));
% speed of the first circle
ptcls.v(1:N/2,1) = -1;

grd_to_ptcl = sint.init_ptcl_mesh(grd, ptcls);

%% FUNCTIONS

function F = forceCells(grd, ptcls, grd_to_ptcl, epsilon, sigma, lennardJonesForce)
    F = zeros(size(ptcls.x,2));

    % size of the grid
    nL1 = size(grd.x,1);
    nL2 = size(grd.x,2);
    Ncells = nL1*nL2;

    index = cellfun(@numel, grd_to_ptcl, 'UniformOutput', true);
    indexNonZero = find(index);

    % loop over all cells
    for i = indexNonZero(:)'
        % nw vertex
        if i == 1
            adCells = [i, i+1, i+nL1, i+nL1+1, i+nL1-1];
        % ne vertex
        elseif i == nL1
            adCells = [i, i-1, i+nL1, i+nL1-1, i+nL1-2];
        % sw vertex
        elseif i == Ncells-nL1+1
            adCells = [i, i+1, i-nL1, i-nL1+1, i-nL1-1];
        % se vertex
        elseif i == Ncells
            adCells = [i, i-1, i-nL1, i-nL1-1, i-nL1+1];
        % n edge
        elseif i > 1 && i < nL1
            adCells = [i, i+1, i-1, i+nL1, i+nL1+1, i+nL1-1];
        % s edge
        elseif i > Ncells-nL1 && i < Ncells
            adCells = [i, i+1, i-1, i-nL1, i-nL1+1, i-nL1-1];
        % w edge
        elseif mod(i-1,nL1) == 0
            adCells = [i, i+1, i+nL1, i-nL1, i+nL1+1, i-nL1+1];
        % e edge
        elseif mod(i,nL1) == 0
            adCells = [i, i-1, i+nL1, i-nL1, i+nL1-1, i-nL1-1];
        % interior cells
        else
            adCells = [i, i+1, i-1, i+nL1, i-nL1, i+nL1+1, i+nL1-1, i-nL1+1, i-nL1-1];
        end
    end

    particlesToCalc = grd_to_ptcl(adCells);
    particlesToCalc = particlesToCalc{:};

    % calculate forces between particles
    Flocal = F;
    Flocal(particlesToCalc) = lennardJonesForce(ptcls.x(:,[particlesToCalc particlesToCalc+size(ptcls.x,2)/2]), sigma, epsilon);
    F = F + Flocal;
end

function F = lennardJonesForce(q, sigmaij, epsij)
    n = length(q);

    % calculate distance
    dx = q(1:n/2) - q(n/2+1:end);
    dy = q(n/2+1:end) - q(n/2+1:end);
    r2 = dx.^2 + dy.^2 + eye(n);

    % calculate sigma2 - 6
    sigma2 = (sigmaij^2)./r2;
    sigma6 = sigma2.*sigma2.*sigma2;

    % calculate force
    Fmat = 48*epsij*sigma6.*(sigma6 - 0.5)./r2;
    Fx = -Fmat.*dx;
    Fy = -Fmat.*dy;

    % sum forces as a column vector
    F = [sum(Fx,2); sum(Fy,2)];
end

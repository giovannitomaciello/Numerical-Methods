clear all
close all
clc

M       = 128;              % --- Number of Fourier harmonics along x (should be a multiple of 2)  
N       = 128;              % --- Number of Fourier harmonics along y (should be a multiple of 2)  
P       = 128;              % --- Number of Fourier harmonics along z (should be a multiple of 2)  
Lx      = 2;               % --- Domain size along x
Ly      = 2;               % --- Domain size along y
Lz      = 2;               % --- Domain size along z
sigma   = 0.1;             % --- Characteristic width of f (make << 1)

% --- Wavenumbers
kx = (2 * pi / Lx) * [0 : (M / 2 - 1) (- M / 2) : (-1)]; % --- Wavenumbers along x
ky = (2 * pi / Ly) * [0 : (N / 2 - 1) (- N / 2) : (-1)]; % --- Wavenumbers along y
kz = (2 * pi / Lz) * [0 : (P / 2 - 1) (- P / 2) : (-1)]; % --- Wavenumbers along z
[Kx, Ky, Kz] = ndgrid(kx, ky, kz); 

% --- Right-hand side of differential equation
hx              = Lx / M;                   % --- Grid spacing along x
hy              = Ly / N;                   % --- Grid spacing along y
hz              = Lz / P;                   % --- Grid spacing along z
x               = (0 : (M - 1)) * hx;
y               = (0 : (N - 1)) * hy;
z               = (0 : (P - 1)) * hz;
[X, Y, Z]       = ndgrid(x, y, z);
rSquared        = (X - 0.5 * Lx).^2 + (Y - 0.5 * Ly).^2 + (Z - 0.5 * Lz).^2;
sigmaSquared    = sigma^2;
f               = exp(-rSquared / (2 * sigmaSquared)) .* (rSquared - 3 * sigmaSquared) / (sigmaSquared^2);
fHat            = fftn(f);

% --- Denominator of the unknown spectrum
den             = -(Kx.^2 + Ky.^2 + Kz.^2); 
den(1, 1, 1)    = 1;            % --- Avoid division by zero at wavenumber (0, 0, 0)

% --- Unknown determination
uHat            = ifftn(fHat ./ den);
u               = real(uHat);
u               = u - u(1,1,1);   % --- Force arbitrary constant to be zero by forcing u(1, 1, 1) = 0

% --- Plots
uRef    = exp(-rSquared / (2 * sigmaSquared));
err     = 100 * sqrt(sum(abs(u(:) - uRef(:)).^2) / sum(abs(uRef(:)).^2));
errMax  = norm(u(:)-uRef(:),inf);
fprintf('Percentage root mean square error = %f\n', err);
fprintf('Maximum error = %f\n', errMax);

% --- Plotting the solution at a fixed z-plane (e.g., z = Lz/2)
zPlane = round(P / 2);
[X2D, Y2D] = meshgrid(x, y);
surf(X2D, Y2D, squeeze(u(:,:,zPlane)))
xlabel('x')
ylabel('y')
zlabel('u')
title('Solution of 3D Poisson equation by spectral method at z = Lz/2')

% --- Plotting the solution at a fixed z-plane (e.g., z = Lz/2)
figure
xPlane = round(P / 2);
[Y2D, Z2D] = meshgrid(z, y);
surf(Y2D, Z2D, squeeze(u(:,:,xPlane)))
xlabel('x')
ylabel('y')
zlabel('u')
title('Solution of 3D Poisson equation by spectral method at x = Lz/2')

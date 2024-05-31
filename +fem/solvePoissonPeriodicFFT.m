function u = solvePoissonPeriodicFFT(A, RHS)
    % Solve the Poisson equation using the finite element method
    % A: the matrix of coefficients
    % RHS: the right-hand side vector

    M       = A.Nx;              % --- Number of Fourier harmonics along x (should be a multiple of 2)  
    N       = A.Ny;              % --- Number of Fourier harmonics along y (should be a multiple of 2)  
    P       = A.Nz;              % --- Number of Fourier harmonics along z (should be a multiple of 2)  
    Lx      = A.Lx;               % --- Domain size along x
    Ly      = A.Ly;               % --- Domain size along y
    Lz      = A.Lz;               % --- Domain size along z

    % --- Wavenumbers
    kx = (2 * pi / Lx) * [0 : (M / 2 - 1) (- M / 2) : (-1)]; % --- Wavenumbers along x
    ky = (2 * pi / Ly) * [0 : (N / 2 - 1) (- N / 2) : (-1)]; % --- Wavenumbers along y
    kz = (2 * pi / Lz) * [0 : (P / 2 - 1) (- P / 2) : (-1)]; % --- Wavenumbers along z
    [Kx, Ky, Kz] = ndgrid(kx, ky, kz); 

    fHat            = fftn(reshape(RHS,A.Nx,A.Ny,A.Nz));

    % --- Denominator of the unknown spectrum
    den             = -(Kx.^2 + Ky.^2 + Kz.^2); 
    den(1, 1, 1)    = 1;            % --- Avoid division by zero at wavenumber (0, 0, 0)

    % --- Unknown determination
    uHat            = ifftn(fHat ./ den);
    u               = real(uHat);
    u               = u - u(1,1,1);   % --- Force arbitrary constant to be zero by forcing u(1, 1, 1) = 0
end
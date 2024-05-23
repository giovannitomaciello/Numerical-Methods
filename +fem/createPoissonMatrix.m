function out = createPoissonMatrix(Nx, Ny, Nz, hx, hy, hz)
    % Create a Poisson matrix for a 3D grid with Nx, Ny, Nz points in each
    % direction.
    % define lenghts
    M = Nx;
    L = Nx*Ny;
    K = Nx*Ny*Nz;
    e = ones(K,1);

    % Define the coefficients
    Cw = -1/hx^2*e; % Coefficient for the west point
    Ce = -1/hx^2*e; % Coefficient for the east point
    Cs = -1/hy^2*e; % Coefficient for the south point
    Cn = -1/hy^2*e; % Coefficient for the north point
    Cb = -1/hz^2*e; % Coefficient for the back point
    Cf = -1/hz^2*e; % Coefficient for the front point
    Cp = 2/hx^2*e + 2/hy^2*e + 2/hz^2*e; % Coefficient for the central point

    % Create the matrix
    out.A = spdiags([Cf, Cw, Cn, Cp, Cs, Ce, Cb], [-L, -M, -1:1, M, L], K, K);
    out.Nx = Nx; out.Ny = Ny; out.Nz = Nz;
    out.hx = hx; out.hy = hy; out.hz = hz;
end
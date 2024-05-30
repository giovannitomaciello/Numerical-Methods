function u = solvePoissonPeriodic(A, RHS, u)
    % Solve the Poisson equation using the finite element method
    % A: the matrix of coefficients
    % RHS: the right-hand side vector

    % Solve the system of equations
    uu = A.A \ RHS;
    uu = reshape(uu,A.Nx-1,A.Ny-1,A.Nz-1);
    u(1:A.Nx-1,1:A.Ny-1,1:A.Nz-1) = uu;
    u(A.Nx,:,:) = u(1,:,:);
    u(:,A.Ny,:) = u(:,1,:);
    u(:,:,A.Nz) = u(:,:,1);
end
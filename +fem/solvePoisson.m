function u = solvePoisson(A, RHS, inodes, bnodes, u)
    % Solve the Poisson equation using the finite element method
    % A: the matrix of coefficients
    % RHS: the right-hand side vector
    % inodes: the interior nodes
    % bnodes: the boundary nodes
    % bvals: the boundary values
    
    % Solve the system of equations
    u(inodes) = A.A(inodes, inodes) \ (RHS(inodes) - A.A(inodes, bnodes) * u(bnodes));
end
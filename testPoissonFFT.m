n = 200; % Size of the grid for each dimension

tic
% Compute eigenvalues for the 1D Laplacian
lambda = zeros(n, 1);
for i = 1:n
    lambda(i) = 2 * (1 - cos(pi * i / (n + 1)));
end

% Construct the 3D Laplacian eigenvalues
lambda3D = bsxfun(@plus, bsxfun(@plus, lambda, permute(lambda, [2, 3, 1])), permute(lambda, [3, 1, 2]));

% Generate a non-random right-hand side (Gaussian distribution centered in the domain)
% 9x9x9 domain
[X, Y, Z] = ndgrid(linspace(1, 10, n), linspace(1, 10, n), linspace(1, 10, n));
sigma = 1 / 5; % Standard deviation of the Gaussian
center = 5.5; % Center of the Gaussian
B = exp(-((X - center).^2 + (Y - center).^2 + (Z - center).^2) / (2 * sigma^2));

% Apply FFT to the right-hand side with appropriate padding and scaling
c = -sqrt(2 / (n + 1));
B_ext = padarray(B, [1, 1, 1], 0, 'both');
B_prime = c^3 * imag(fft(B_ext, [], 1));
B_prime = imag(fft(B_prime(2:n+1, :, :), [], 2));
B_prime = imag(fft(B_prime(:, 2:n+1, :), [], 3));
B_prime = B_prime(:, :, 2:n+1);

% Solve the diagonal system
U_prime = B_prime ./ lambda3D;

% Apply inverse FFT to get back to spatial domain
U_prime_ext = padarray(U_prime, [1, 1, 1], 0, 'both');
U = c^3 * imag(fft(U_prime_ext, [], 1));
U = imag(fft(U(2:n+1, :, :), [], 2));
U = imag(fft(U(:, 2:n+1, :), [], 3));
U = U(:, :, 2:n+1);
toc
% Plot a slice of the solution
slice = U(:, :, round(n / 2));
imagesc(slice);
title('Slice of the solution U at the middle plane z');
colorbar;

% Plot a slice of the solution
figure
slice = squeeze(U(:, round(n / 2), :));
imagesc(slice);
title('Slice of the solution U at the middle plane y');
colorbar;

% Plot a slice of RHS
slice = B(:, :, round(n / 2));
imagesc(slice);
title('Slice of the RHS at the middle plane z');
colorbar;

% Plot a slice of RHS
figure
slice = squeeze(B(:, round(n / 2), :));
imagesc(slice);
title('Slice of the RHS at the middle plane y');
colorbar;
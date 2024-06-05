function out = Interp3ptcls(X,Y,Z,V,x,y,z,dim)

if dim == 'x'
    [~, idx] = min(abs(X(1,:,1) - x'), [], 2); % Find indices of nearest X grid points
    [Y_grid,Z_grid] = meshgrid(linspace(min(Y(:)),max(Y(:)),size(X,2)), linspace(min(Z(:)),max(Z(:)),size(X,3)));
    out = zeros(size(x));

    for i = 1:length(x)
        % Extract a 2D slice at the nearest X.
        V_slice = squeeze(V(idx(i),:,:));
        
        % Linear interpolation on the slice along Y and Z
        out(i) = interp2(Y_grid, Z_grid, V_slice, y(i), z(i), 'linear');
    end
elseif dim == 'y'
    [~, idx] = min(abs(Y(:,1,1) - y), [], 1); % Find indices of nearest Y grid points
    [X_grid,Z_grid] = meshgrid(linspace(min(X(:)),max(X(:)),size(X,1)), linspace(min(Z(:)),max(Z(:)),size(X,3)));
    out = zeros(size(y));

    for i = 1:length(y)
        % Extract a 2D slice at the nearest Y.
        V_slice = squeeze(V(:,idx(i),:));
        
        % Linear interpolation on the slice along X and Z
        out(i) = interp2(X_grid, Z_grid, V_slice, x(i), z(i), 'linear');
    end
elseif dim == 'z'
    [~, idx] = min(abs(Z(1,1,:) - z'), [], 2); % Find indices of nearest Z grid points
    [X_grid,Y_grid] = meshgrid(linspace(min(X(:)),max(X(:)),size(X,1)), linspace(min(Y(:)),max(Y(:)),size(X,2)));
    out = zeros(size(z));

    for i = 1:length(z)
        % Extract a 2D slice at the nearest Z.
        V_slice = squeeze(V(:, :, idx(i)));
        
        % Linear interpolation on the slice along X and Y
        out(i) = interp2(X_grid, Y_grid, V_slice, x(i), y(i), 'linear');
    end
end
end


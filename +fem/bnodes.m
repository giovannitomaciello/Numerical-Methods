function index = bnodes(face, X, Y, Z)
    % faces w, e, n, s, b, f
    switch face
        case 'w'
            indexw = find(X == min(X,[],"all"));
            index = indexw(:);
        case 'e'
            indexe = find(X == max(X,[],"all"));
            index = indexe(:);
        case 'n'
            indexn = find(Y == max(Y,[],"all"));
            index = indexn(:);
        case 's'
            indexs = find(Y == min(Y,[],"all"));
            index = indexs(:);
        case 'b'
            indexb = find(Z == min(Z,[],"all"));
            index = indexb(:);
        case 'f'
            indexf = find(Z == max(Z,[],"all"));
            index = indexf(:);
    end
end
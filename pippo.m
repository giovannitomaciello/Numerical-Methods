N = 100;
M = (N+1)^3;
x = linspace(0,1,N+1);
y = x; z=x;

[X , Y, Z] = meshgrid(x,y,z);

u0 = zeros((N+1)^3,1);
u0(Z==0) = 1;
u0_slice = reshape(u0,[N+1,N+1,N+1]);
u0_slice = squeeze(u0_slice(:,5,:));
mesh(squeeze(Y(:,5,:)),squeeze(Z(:,5,:)),u0_slice)

A =spdiags([-1 -1 8 -1 -1],[-N-1,-1:1,N+1],M,M);
spy(A)

u = ones(M,1); b = A*u;
r = chol(A);
tic
y = r'\b; 
u1 = r\y; 
toc

tic
u2 = A\b;
toc


clc; clear all; close all;
set(0,'DefaultFigureColor',[1 1 1]);
%gunzip poissonCollisionData\c.mat.gz
%gunzip poissonCollisionData\q.mat.gz
%gunzip poissonCollisionData\p.mat.gz
load("./poissonCollisionData/coordLarge/p.mat")
load("./poissonCollisionData/coordLarge/c.mat")
load("./poissonCollisionData/coordLarge/q.mat")
%%
img = imread("./poissonCollisionData/theEnd.png");
% convert to grey scale
img = rgb2gray(img);
% convert to double
img = double(img);
% flip the image
imgTE = flipud(img);

img = imread("./poissonCollisionData/ThisIs.png");
% convert to grey scale
img = rgb2gray(img);
% convert to double
img = double(img);
% flip the image
imgTI = flipud(img);

% interpolate latest time of the simulation from the image 1
qLatest = q(:,:,end);
x = linspace(0,64,size(imgTE,2));
y = linspace(0,64,size(imgTE,1));
[X,Y] = meshgrid(y,x);
qLatestInterp1 = interp2(X,Y,imgTE,q(1,:,end),q(2,:,end));

% interpolate latest time of the simulation from the image 2
x = linspace(0,64,size(imgTI,2));
y = linspace(0,64,size(imgTI,1));
[X,Y] = meshgrid(y,x);
qLatestInterp2 = interp2(X,Y,imgTI,q(1,:,end),q(3,:,end));

qLatestInterp = sqrt(abs(qLatestInterp2 + qLatestInterp1));
qLatestInterp(qLatestInterp>mean(qLatestInterp)) = mean(qLatestInterp);

%% plot the simulation and save video for the end
figure
v = VideoWriter('poissonCollisionTheEnd.avi');
open(v)
for i = 1:size(q,3)
    scatter3(q(1,:,i),q(2,:,i),q(3,:,i),20,qLatestInterp,'filled','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2)
    title("Initial Random Value?")
    xlabel("x")
    ylabel("y")
    zlabel("z")
    view(0,90)
    axis equal
    colormap("prism")
    % cut the image to the size of the simulation
    xlim([0 64])
    ylim([0 64])
    zlim([0 64])
    drawnow
    frame = getframe(gcf);
    writeVideo(v,frame);
end
close(v)

%% plot the simulation and save video for the end
figure
v = VideoWriter('poissonCollisionThisIs.avi');
open(v)
for i = 1:size(q,3)
    scatter3(q(1,:,i),q(2,:,i),q(3,:,i),20,qLatestInterp,'filled','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2)
    title("Initial Random Value?")
    xlabel("x")
    ylabel("y")
    zlabel("z")
    view(0,0)
    axis equal
    colormap("prism")
    % cut the image to the size of the simulation
    xlim([0 64])
    ylim([0 64])
    zlim([0 64])
    drawnow
    frame = getframe(gcf);
    writeVideo(v,frame);
end
close(v)
close(v)

%% plot the simulation and save video for the charge for x and y
figure
v = VideoWriter('poissonCollisionChargeXY.avi');
open(v)
for i = 1:1:size(q,3)
    scatter(q(1,:,i),q(2,:,i),1.5,c,'filled')
    title("Charge Distribution XY Plane")
    xlabel("x")
    ylabel("y")
    axis equal
    colormap(winter(2))
    colorbar
    % cut the image to the size of the simulation
    xlim([0 64])
    ylim([0 64])
    drawnow
    frame = getframe(gcf);
    writeVideo(v,frame);
end
close(v)

%% plot the simulation and save video for the charge for x and z
figure
v = VideoWriter('poissonCollisionChargeXZ.avi');
open(v)
for i = 1:1:size(q,3)
    scatter(q(1,:,i),q(3,:,i),1.5,c,'filled')
    title("Charge Distribution XZ Plane")
    xlabel("x")
    ylabel("z")
    axis equal
    colormap(winter(2))
    colorbar
    % cut the image to the size of the simulation
    xlim([0 64])
    ylim([0 64])
    drawnow
    frame = getframe(gcf);
    writeVideo(v,frame);
end
close(v)

%% plot the simulation and save video for the charge for y and z
figure
v = VideoWriter('poissonCollisionChargeYZ.avi');
open(v)
for i = 1:1:size(q,3)
    scatter(q(2,:,i),q(3,:,i),1.5,c,'filled')
    title("Charge Distribution YZ Plane")
    xlabel("y")
    ylabel("z")
    axis equal
    colormap(winter(2))
    colorbar
    % cut the image to the size of the simulation
    xlim([0 64])
    ylim([0 64])
    drawnow
    frame = getframe(gcf);
    writeVideo(v,frame);
end
close(v)

%% plot the simulation and save video for the momentum for x and y
figure
v = VideoWriter('poissonCollisionMomentumXY.avi');
open(v)
for i = 1:1:size(q,3)
    scatter(q(1,:,i),q(2,:,i),1.5,vecnorm(p(:,:,i)),'filled')
    title("Momentum Distribution XY Plane")
    xlabel("x")
    ylabel("y")
    axis equal
    colormap("parula")
    colorbar
    % cut the image to the size of the simulation
    xlim([0 64])
    ylim([0 64])
    drawnow
    frame = getframe(gcf);
    writeVideo(v,frame);
end
close(v)

%% plot the simulation and save video for the momentum for x and z
figure
v = VideoWriter('poissonCollisionMomentumXZ.avi');
open(v)
for i = 1:1:size(q,3)
    scatter(q(1,:,i),q(3,:,i),1.5,vecnorm(p(:,:,i)),'filled')
    title("Momentum Distribution XZ Plane")
    xlabel("x")
    ylabel("z")
    axis equal
    colormap("parula")
    colorbar
    % cut the image to the size of the simulation
    xlim([0 64])
    ylim([0 64])
    drawnow
    frame = getframe(gcf);
    writeVideo(v,frame);
end
close(v)

%% plot the simulation and save video for the momentum for y and z
figure
v = VideoWriter('poissonCollisionMomentumYZ.avi');
open(v)
for i = 1:1:size(q,3)
    scatter(q(2,:,i),q(3,:,i),1.5,vecnorm(p(:,:,i)),'filled')
    title("Momentum Distribution YZ Plane")
    xlabel("y")
    ylabel("z")
    axis equal
    colormap("parula")
    colorbar
    % cut the image to the size of the simulation
    xlim([0 64])
    ylim([0 64])
    drawnow
    frame = getframe(gcf);
    writeVideo(v,frame);
end
close(v)

%% plot the mesh of the simulation for phi
x = linspace(0,64,64);
y = linspace(0,64,64);
z = linspace(0,64,64);
[X,Y,Z] = meshgrid(x,y,z);
figure
v = VideoWriter('poissonCollisionMeshXYphi.avi');
open(v)
A.Nx = 64; A.Ny = 64; A.Nz = 64;
A.Lx = 64; A.Ly = 64; A.Lz = 64;
rCut = 2; epsilon = .1;
for i = 1:size(q,3)
    % calculate charge density
    % higly optimized with eigen3, hashing and IntelTBB
    rho_lr = sint.ptclsToMeshInterp(X, Y, Z, c, q(:,:,i), rCut);
    RHS = reshape(rho_lr,[],1)/epsilon;
    phi = fem.solvePoissonPeriodicFFT(A, RHS);
    phi = reshape(phi, size(X));
    % plot surface
    imagesc(phi(:,:,32))
    colormap("parula")
    colorbar
    title("phi XY Plane")
    xlabel("x")
    ylabel("y")
    axis equal
    shading interp
    %cut the image to the size of the simulation
    drawnow
    frame = getframe(gcf);
    writeVideo(v,frame);
end
close(v)

%% plot the mesh of the simulation for phi
x = linspace(0,64,64);
y = linspace(0,64,64);
z = linspace(0,64,64);
[X,Y,Z] = meshgrid(x,y,z);
figure
v = VideoWriter('poissonCollisionMeshXZphi.avi');
open(v)
A.Nx = 64; A.Ny = 64; A.Nz = 64;
A.Lx = 64; A.Ly = 64; A.Lz = 64;
rCut = 2; epsilon = .1;
for i = 1:1:size(q,3)
    % calculate charge density
    % higly optimized with eigen3, hashing and IntelTBB
    rho_lr = sint.ptclsToMeshInterp(X, Y, Z, c, q(:,:,i), rCut);
    RHS = reshape(rho_lr,[],1)/epsilon;
    phi = fem.solvePoissonPeriodicFFT(A, RHS);
    phi = reshape(phi, size(X));
    % plot surface
    imagesc(squeeze(phi(:,32,:)))
    title("phi XZ Plane")
    xlabel("x")
    ylabel("z")
    axis equal
    colormap("parula")
    colorbar
    shading interp
    % cut the image to the size of the simulation
    drawnow
    frame = getframe(gcf);
    writeVideo(v,frame);
end
close(v)

%% plot the mesh of the simulation for phi
x = linspace(0,64,64);
y = linspace(0,64,64);
z = linspace(0,64,64);
[X,Y,Z] = meshgrid(x,y,z);
figure
v = VideoWriter('poissonCollisionMeshYZphi.avi');
open(v)
A.Nx = 64; A.Ny = 64; A.Nz = 64;
A.Lx = 64; A.Ly = 64; A.Lz = 64;
rCut = 2; epsilon = .1;
for i = 1:1:size(q,3)
    % calculate charge density
    % higly optimized with eigen3, hashing and IntelTBB
    rho_lr = sint.ptclsToMeshInterp(X, Y, Z, c, q(:,:,i), rCut);
    RHS = reshape(rho_lr,[],1)/epsilon;
    phi = fem.solvePoissonPeriodicFFT(A, RHS);
    phi = reshape(phi, size(X));
    % plot surface
    imagesc(squeeze(phi(32,:,:)))
    title("phi YZ Plane")
    xlabel("y")
    ylabel("z")
    axis equal
    colormap("parula")
    colorbar
    shading interp
    % cut the image to the size of the simulation
    drawnow
    frame = getframe(gcf);
    writeVideo(v,frame);
end
close(v)

%% plot the mesh of the simulation for rho
x = linspace(0,64,64);
y = linspace(0,64,64);
z = linspace(0,64,64);
[X,Y,Z] = meshgrid(x,y,z);
figure
v = VideoWriter('poissonCollisionMeshXYcharge.avi');
open(v)
A.Nx = 64; A.Ny = 64; A.Nz = 64;
A.Lx = 64; A.Ly = 64; A.Lz = 64;
rCut = 2; epsilon = .1;
for i = 1:1:size(q,3)
    % calculate charge density
    % higly optimized with eigen3, hashing and IntelTBB
    rho_lr = sint.ptclsToMeshInterp(X, Y, Z, c, q(:,:,i), rCut);
    rho_lr = reshape(rho_lr,64,64,64);
    % plot surface
    imagesc(rho_lr(:,:,32))
    title("rho XY Plane")
    xlabel("x")
    ylabel("y")
    axis equal
    colormap("jet") % Change colormap to jet
    colorbar
    shading interp
    % cut the image to the size of the simulation
    drawnow
    frame = getframe(gcf);
    writeVideo(v,frame);
end
close(v)

%% plot the mesh of the simulation for rho
x = linspace(0,64,64);
y = linspace(0,64,64);
z = linspace(0,64,64);
[X,Y,Z] = meshgrid(x,y,z);
figure
v = VideoWriter('poissonCollisionMeshXZcharge.avi');
open(v)
A.Nx = 64; A.Ny = 64; A.Nz = 64;
A.Lx = 64; A.Ly = 64; A.Lz = 64;
rCut = 2; epsilon = .1;
for i = 1:1:size(q,3)
    % calculate charge density
    % higly optimized with eigen3, hashing and IntelTBB
    rho_lr = sint.ptclsToMeshInterp(X, Y, Z, c, q(:,:,i), rCut);
    rho_lr = reshape(rho_lr,64,64,64);
    % plot surface
    imagesc(squeeze(rho_lr(:,32,:)))
    title("rho XZ Plane")
    xlabel("x")
    ylabel("z")
    axis equal
    colormap("jet") % Change colormap to jet
    colorbar
    shading interp
    % cut the image to the size of the simulation
    drawnow
    frame = getframe(gcf);
    writeVideo(v,frame);
end
close(v)

%% plot the mesh of the simulation for rho
x = linspace(0,64,64);
y = linspace(0,64,64);
z = linspace(0,64,64);
[X,Y,Z] = meshgrid(x,y,z);
figure
v = VideoWriter('poissonCollisionMeshYZcharge.avi');
open(v)
A.Nx = 64; A.Ny = 64; A.Nz = 64;
A.Lx = 64; A.Ly = 64; A.Lz = 64;
rCut = 2; epsilon = .1;
for i = 1:1:size(q,3)
    % calculate charge density
    % higly optimized with eigen3, hashing and IntelTBB
    rho_lr = sint.ptclsToMeshInterp(X, Y, Z, c, q(:,:,i), rCut);
    rho_lr = reshape(rho_lr,64,64,64);
    % plot surface
    imagesc(squeeze(rho_lr(32,:,:)))
    title("rho YZ Plane")
    xlabel("y")
    ylabel("z")
    axis equal
    colormap("jet") % Change colormap to jet
    colorbar
    shading interp
    % cut the image to the size of the simulation
    drawnow
    frame = getframe(gcf);
    writeVideo(v,frame);
end
close(v)

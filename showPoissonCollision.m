clc; clear all; close all;
set(0,'DefaultFigureColor',[1 1 1]);
gunzip poissonCollisionData\c.mat.gz
gunzip poissonCollisionData\q.mat.gz
gunzip poissonCollisionData\p.mat.gz
load("./poissonCollisionData/p.mat")
load("./poissonCollisionData/q.mat")
load("./poissonCollisionData/c.mat")
img = imread("./poissonCollisionData/theEnd.png");
% convert to grey scale
img = rgb2gray(img);
% convert to double
img = double(img);
% flip the image
imgTE = flipud(img);

% interpolate latest time of the simulation from the image
qLatest = q(:,:,end);
x = linspace(0,64,size(imgTE,2));
y = linspace(0,64,size(imgTE,1));
[X,Y] = meshgrid(y,x);
qLatestInterp = interp2(X,Y,imgTE,q(1,:,end),q(2,:,end));

%% plot the simulation and save video for the end
figure
v = VideoWriter('poissonCollisionTheEnd.avi');
open(v)
for i = 1:1:size(q,3)
    scatter(q(1,:,i),q(2,:,i),3,qLatestInterp,'filled')
    title("Initial Random Value?")
    xlabel("x")
    ylabel("y")
    axis equal
    colormap("prism")
    % cut the image to the size of the simulation
    xlim([0 64])
    ylim([0 64])
    drawnow
    frame = getframe(gcf);
    writeVideo(v,frame);
end
close(v)

%% plot the simulation and save video for the charge for x and y
figure
v = VideoWriter('poissonCollisionChargeXY.avi');
open(v)
for i = 1:1:size(q,3)
    scatter(q(1,:,i),q(2,:,i),3,c,'filled')
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
    scatter(q(1,:,i),q(3,:,i),3,c,'filled')
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
    scatter(q(2,:,i),q(3,:,i),3,c,'filled')
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
    scatter(q(1,:,i),q(2,:,i),3,vecnorm(p(:,:,i)),'filled')
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
    scatter(q(1,:,i),q(3,:,i),3,vecnorm(p(:,:,i)),'filled')
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
    scatter(q(2,:,i),q(3,:,i),3,vecnorm(p(:,:,i)),'filled')
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




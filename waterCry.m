clc; clear all; close all;

%% water crystal parameters
HOdist = 0.96; % oxygen atom
HOangle = 104.5; % hydrogen atom
HHdist = 0.96*sind(HOangle/2)*2;

%% 3D water molecules
positionNoRot = [0 HOdist*cosd(HOangle/2); -HOdist*sind(HOangle/2) 0; HOdist*sind(HOangle/2) 0];
shift = [0 HOdist*cosd(HOangle/2)];
positionNoRot = positionNoRot - shift;

%% rotate by 60 deg and add to the original position
position1layer = [];
for i = 1:6
    positionRot = positionNoRot*[cosd(120*i + 120) -sind(120*i + 120); sind(120*i + 120) cosd(120*i + 120)];
    positionTrasl = positionRot + [0 2.2; 0 2.2; 0 2.2];
    position1layer = [positionTrasl*[cosd(60*i) -sind(60*i); sind(60*i) cosd(60*i)]; position1layer];
end
position1layer = [position1layer , -ones(length(position1layer),1)];

position2layer = [];
for i = 1:6
    positionRot = positionNoRot*[cosd(120*i - 120) -sind(120*i - 120); sind(120*i - 120) cosd(120*i - 120)];
    positionTrasl = positionRot + [0 2.2; 0 2.2; 0 2.2];
    position2layer = [positionTrasl*[cosd(60*i) -sind(60*i); sind(60*i) cosd(60*i)]; position2layer];
end
position2layer = [position2layer , zeros(length(position2layer),1)];

position3layer = [];
for i = 1:6
    positionRot = positionNoRot*[cosd(120*i) -sind(120*i); sind(120*i) cosd(120*i)];
    positionTrasl = positionRot + [0 2.2; 0 2.2; 0 2.2];
    position3layer = [positionTrasl*[cosd(60*i) -sind(60*i); sind(60*i) cosd(60*i)]; position3layer];
end
position3layer = [position3layer , ones(length(position3layer),1)];

q0 = [position1layer; position2layer; position3layer];
p0 = zeros(size(q0));

%% constraint matrix
C = sparse(size(q0,1),size(q0,1));
connectivity = C;
for i = 3:3:size(q0,1)
    C(i-2,i-1) = HOdist;
    C(i-2,i) = HOdist;
    C(i-1,i) = HHdist;

    connectivity(i-2,i-1) = HOdist;
    connectivity(i-2,i) = HOdist;
end

 C = C + C';

gplot3(connectivity,q0,'lineWidth',2,'Color','k')
hold on
scatter3(q0(1:3:end,1),q0(1:3:end,2),q0(1:3:end,3),80,"red","filled")
scatter3(q0(2:3:end,1),q0(2:3:end,2),q0(2:3:end,3),60,"blue","filled")
scatter3(q0(3:3:end,1),q0(3:3:end,2),q0(3:3:end,3),60,"blue","filled")
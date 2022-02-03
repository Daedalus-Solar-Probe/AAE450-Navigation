%% Clear

clear
clc

%% Setup

launch = juliandate(datetime("now")); % [JD]
tof = 200; % [days]
arrive = launch + tof; % [JD]

%% Calc

% earth
[Re,Ve] = planetEphemeris(launch,'Sun','Earth'); % [km,km/s]

% mars
[Rm,Vm] = planetEphemeris(launch,'Sun','Mars'); % [km,km/s]

%% Plot

figure(1)
plot3(0,0,0,'ko','MarkerSize',20,'MarkerFaceColor','y')
hold on
plot3(Re(1),Re(2),Re(3),'ko','MarkerSize',10,'MarkerFaceColor','b')
plot3(Rm(1),Rm(2),Rm(3),'ko','MarkerSize',10,'MarkerFaceColor','r')
hold off
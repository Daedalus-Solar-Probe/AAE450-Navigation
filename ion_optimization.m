%% Clear

clear
clc

%% Constants

% solar standard gravitational parameter
mu = 1.327e11; % [km^3/s^2]

% AU in km
AU = 4.496e8; % [km/AU]

% %% Inputs
% 
% % semimajor axis
% a = 1*AU; % [km]
% 
% % 15 degree inclination change
% di = deg2rad(15); % [rad]
% 
% %% Inclination Change
% 
% % circular orbit
% e = 0; % [-]
% 
% % argument of periapsis
% omega = 0; % [rad]
% 
% % true anomaly
% % f = 0; % [rad]
% f = linspace(0,2*pi); % [rad]
% 
% % mean motion
% n = sqrt(mu/a^3); % [rad/s]
% 
% % case 1
% dV1 = 2*sin(di/2)*sqrt(1-e^2)*cos(omega+f)*n*a./(1+e*cos(f)); % [km/s]
% 
% % case 2
% dV2 = 2*sqrt(mu/a)*sin(di/2); % [km/s]



%%

% start at Earth
r1 = 1*AU; % [km]

% end at venus
r2 = 0.723*AU; % [km]

% delta-V
dv12 = sqrt(mu/r2) - sqrt(mu/r1); % [km/s]

%% Clear

clear
clc

%% Constants

% solar standard gravitational parameter
mu = 1.327e11; % [km^3/s^2]

% AU in km
AU = 1.496e+8; % [km/AU]

%% Inputs

r0 = [1*AU; 0; 0]; % [km]
v0 = [0; sqrt(mu/AU); 0]; % [km/s]

y0 = [r0; v0];

tspan = [0 inf];

accel = 8e-8; % [km/s^2]

r_desired = 0.723*AU;

opts = odeset("RelTol",1e-6,"AbsTol",1e-6,"Events",@(t,y)eventsfun1(t,y,r_desired));
[t1,y1] = ode45(@(t,y)odefun1(t,y,mu,accel),tspan,y0,opts);

dV1 = accel*max(t1) % [km/s]
tof1 = max(t1)/86400 % [days]



%% Plots

% trajectory
figure(1)
plot3(y1(:,1),y1(:,2),y1(:,3),'r')
grid on
xlim([-5e8 5e8])
ylim([-5e8 5e8])
zlim([-5e8 5e8])

% %% 
% 
% i_desired = deg2rad(90);
% 
% tspan = [0 350*86400];
% 
% % y0 = y1(end,:)';
% 
% opts = odeset("RelTol",1e-6,"AbsTol",1e-6,"Events",@(t,y)eventsfun2(t,y,i_desired));
% [t2,y2] = ode45(@(t,y)odefun2(t,y,mu,accel),tspan,y0,opts);
% 
% dV1 = accel*max(t1) % [km/s]
% tof1 = max(t1)/86400 % [days]
% 
% %%
% 
% figure(2)
% plot3(y2(:,1),y2(:,2),y2(:,3),'b')
% grid on
% xlim([-1e9 1e9])
% ylim([-1e9 1e9])
% zlim([-1e9 1e9])
% 
% dV2 = accel*max(t2) % [km/s]
% tof2 = max(t2)/86400 % [days]

%%

function dydt = odefun1(~,y,mu,accel)

% initialize output
dydt = zeros(size(y));

% position
r = y(1:3); % [km]

% velocity
v = y(4:6); % [km/s]

% % specific angular momentum
% h = cross(r,v); % [km^2/s]

% acceleration due to gravity
a_g = -mu*r/norm(r)^3; % [km/s^2]

% acceleration due to propulsion
a_p = -accel*v/norm(v); % [km/s]

dydt(1:3) = v; % [km/s] r_dot = v
dydt(4:6) = a_g + a_p; % [km/s^2] v_dot = a

end % function

function dydt = odefun2(~,y,mu,accel)

% initialize output
dydt = zeros(size(y));

% position
r = y(1:3); % [km]

% velocity
v = y(4:6); % [km/s]

% specific angular momentum
h = cross(r,v); % [km^2/s]

% acceleration due to gravity
a_g = -mu*r/norm(r)^3; % [km/s^2]

% acceleration due to propulsion
a_p = accel*h/norm(h); % [km/s]

dydt(1:3) = v; % [km/s] r_dot = v
dydt(4:6) = a_g + a_p; % [km/s^2] v_dot = a

end % function

function [val,isterm,dir] = eventsfun1(~,y,r_desired)

r = y(1:3);

val = norm(r) - r_desired;
isterm = 1;
dir = [];

end % function

function [val,isterm,dir] = eventsfun2(~,y,i_desired)

r = y(1:3);
v = y(4:6);

h = cross(r,v);

val = i_desired - h(3)/norm(h);
isterm = 1;
dir = [];

end % function
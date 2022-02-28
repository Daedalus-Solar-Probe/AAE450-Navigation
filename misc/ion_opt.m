%% Clear

clear
clc

%% Constants

% solar standard gravitational parameter
mu = 1.327e11; % [km^3/s^2]

% AU in km
AU = 1.496e+8; % [km/AU]

%% Inputs

% initial semi-major axis
a0 = 1*AU; % [km]

% maneuver semi-major axis
a1 = 1.5*AU; % [km]

% inclination change
i2 = deg2rad(60);

% final semi-major axis
a3 = 0.48*AU; % [km]

% propulsion acceleration
accel = 100e-9; % [km/s]

%% Simulate Spiral Outwards

% parameters
P.mu = mu; % [km^3/s^2]
P.a = a1; % [km]
P.dir = 1; % [-] outwards
P.accel = accel; % [km/s^2]

% initial conditions
t0 = 0; % [s]
r0 = [a0; 0; 0]; % [km]
v0 = [0; sqrt(mu/a0); 0]; % [km/s]

[t1,r1,v1] = propSpiral(t0,r0,v0,P);

%% Simulate Inclination Cranking

% initial conditions
t0 = t1(end); % [s]
r0 = r1(end,:)'; % [km]
v0 = v1(end,:)'; % [km/s]

P.inc = i2; % [rad]

[t2,r2,v2] = propCrank(t0,r0,v0,P);

%% Simulate Spiral Inwards

% initial conditions
t0 = t2(end); % [s]
r0 = r2(end,:)'; % [km]
v0 = v2(end,:)'; % [km/s]

P.a = a3; % [km]
P.dir = -1; % [-] inwards

[t3,r3,v3] = propSpiral(t0,r0,v0,P);

%% Post

% combine data
t = [t1; t2; t3]; % [s] time
r = [r1; r2; r3]; % [km] position
v = [v1; v2; v3]; % [km/s] velocity
h = cross(r,v); % [km^2/s] specific angular momentum

%% Plot

figure(1)
plot3(r(:,1)/AU,r(:,2)/AU,r(:,3)/AU,'b')

xlim([-2 2])
ylim([-2 2])
zlim([-2 2])

title("Ion Trajectory")
xlabel("X [AU]")
ylabel("Y [AU]")
zlabel("Z [AU]")
grid on

figure(2)
quiver3(r(:,1),r(:,2),r(:,3),h(:,1),h(:,2),h(:,3))
grid on

%%

% opts = optimset('Display','iter');
% x = fminbnd(@objfun,1.01,10,opts);

%%

x = linspace(1,5);

for n = 1:length(x)
    x(n)
    f(n) = objfun(x(n));
end

figure(3)
plot(x,f)

%% Objective Function
function f = objfun(x)

mu = 1.327e11; % [km^3/s^2]
a0 = 149600000; % [km]

P.mu = mu; % [km^3/s^2]
P.a = x*149600000; % [km]
P.dir = 1; % [-]
P.accel = 100e-9; % [km/s]

% initial conditions
t0 = 0; % [s]
r0 = [a0; 0; 0]; % [km]
v0 = [0; sqrt(mu/a0); 0]; % [km/s]

[t1,r1,v1] = propSpiral(t0,r0,v0,P);

% initial conditions
t0 = t1(end); % [s]
r0 = r1(end,:)'; % [km]
v0 = v1(end,:)'; % [km/s]

P.inc = deg2rad(65); % [rad]

[t2,r2,v2] = propCrank(t0,r0,v0,P);

% initial conditions
t0 = t2(end); % [s]
r0 = r2(end,:)'; % [km]
v0 = v2(end,:)'; % [km/s]

P.a = 71808000; % [km]
P.dir = -1; % [-] inwards

[t3,~,~] = propSpiral(t0,r0,v0,P);

if t3(end) - 1.01*t0 < 1e9
    disp("oof")
end

f = P.accel*t3(end);

end % function

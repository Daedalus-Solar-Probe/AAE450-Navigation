%% Clear

clear
clc

%% Constants

% solar standard gravitational parameter
mu = 1.327e11; % [km^3/s^2]

% AU in km
AU = 1.496e+8; % [km/AU]

%% Inputs

% initial condition (at Earth in Heliocentric orbit)
R0 = [AU; 0; 0]; % [km]
V0 = [0; sqrt(mu/AU); 0]; % [km/s]

y0 = [R0; V0];

% spacecraft acceleration (constant)
accel = 1e-07; % [km/s]

% simulation time span
tspan = [0 1e9]; % [s]

% integration error
epsilon = 1e-9;

% Parameters
P.mu = mu;
P.a = 4*AU;
P.dir = 1;
P.accel = accel;
P.inc = deg2rad(90);

P.m_dry = 100; % [kg]
P.Isp = 4100; % [s]
P.a0 = accel; % [km/s^2]

%% Initial Orbit

% integration
[t0,r0,v0] = propPeriod(0,R0,V0,mu);

% orbital elements
a0 = zeros(size(t0)); % [km] semi-major axis
e0 = zeros(size(t0)); % [-] eccentricity
omega0 = zeros(size(t0)); % [rad] argument of periapsis
OMEGA0 = zeros(size(t0)); % [rad] longitude of periapsis
i0 = zeros(size(t0)); % [rad] inclination
M0 = zeros(size(t0)); % [rad] mean anomaly

for n = 1:length(t0)

    [a0(n),e0(n),omega0(n),OMEGA0(n),i0(n),M0(n)] = rv2coes(r0(n,:)',v0(n,:)',mu);

end % for

%% Spiral

[t1,r1,v1] = propSpiral(t0(end),r0(end,:)',v0(end,:)',P);

% orbital elements
a1 = zeros(size(t0)); % [km] semi-major axis
e1 = zeros(size(t0)); % [-] eccentricity
omega1 = zeros(size(t0)); % [rad] argument of periapsis
OMEGA1 = zeros(size(t0)); % [rad] longitude of periapsis
i1 = zeros(size(t0)); % [rad] inclination
M1 = zeros(size(t0)); % [rad] mean anomaly

for n = 1:length(t1)

    [a1(n),e1(n),omega1(n),OMEGA1(n),i1(n),M1(n)] = rv2coes(r1(n,:)',v1(n,:)',mu);

end % for

%% Orbit After Spiral

% integration
[t2,r2,v2] = propPeriod(t1(end),r1(end,:)',v1(end,:)',mu);

% orbital elements
a2 = zeros(size(t2)); % [km] semi-major axis
e2 = zeros(size(t2)); % [-] eccentricity
omega2 = zeros(size(t2)); % [rad] argument of periapsis
OMEGA2 = zeros(size(t2)); % [rad] longitude of periapsis
i2 = zeros(size(t2)); % [rad] inclination
M2 = zeros(size(t2)); % [rad] mean anomaly

for n = 1:length(t2)

    [a2(n),e2(n),omega2(n),OMEGA2(n),i2(n),M2(n)] = rv2coes(r2(n,:)',v2(n,:)',mu);

end % for

%% Cranking

% integration
[t3,r3,v3] = propCrank(t2(end),r2(end,:)',v2(end,:)',P);

% orbital elements
a3 = zeros(size(t3)); % [km] semi-major axis
e3 = zeros(size(t3)); % [-] eccentricity
omega3 = zeros(size(t3)); % [rad] argument of periapsis
OMEGA3 = zeros(size(t3)); % [rad] longitude of periapsis
i3 = zeros(size(t3)); % [rad] inclination
M3 = zeros(size(t3)); % [rad] mean anomaly

for n = 1:length(t3)

    [a3(n),e3(n),omega3(n),OMEGA3(n),i3(n),M3(n)] = rv2coes(r3(n,:)',v3(n,:)',mu);

end % for

%% Orbit After Cranking

% integration
[t4,r4,v4] = propPeriod(t3(end),r3(end,:)',v3(end,:)',mu);

% orbital elements
a4 = zeros(size(t4)); % [km] semi-major axis
e4 = zeros(size(t4)); % [-] eccentricity
omega4 = zeros(size(t4)); % [rad] argument of periapsis
OMEGA4 = zeros(size(t4)); % [rad] longitude of periapsis
i4 = zeros(size(t4)); % [rad] inclination
M4 = zeros(size(t4)); % [rad] mean anomaly

for n = 1:length(t4)

    [a4(n),e4(n),omega4(n),OMEGA4(n),i4(n),M4(n)] = rv2coes(r4(n,:)',v4(n,:)',mu);

end % for

%% Spiral Again!

P.dir = -1;
P.target = 0.48*AU;

% integration
[t5,r5,v5] = propSpiral(t4(end),r4(end,:)',v4(end,:)',P);

% orbital elements
a5 = zeros(size(t5)); % [km] semi-major axis
e5 = zeros(size(t5)); % [-] eccentricity
omega5 = zeros(size(t5)); % [rad] argument of periapsis
OMEGA5 = zeros(size(t5)); % [rad] longitude of periapsis
i5 = zeros(size(t5)); % [rad] inclination
M5 = zeros(size(t5)); % [rad] mean anomaly

for n = 1:length(t5)

    [a5(n),e5(n),omega5(n),OMEGA5(n),i5(n),M5(n)] = rv2coes(r5(n,:)',v5(n,:)',mu);

end % for

%% Orbit After Spiraling Again!

% integration
[t6,r6,v6] = propPeriod(t5(end),r5(end,:)',v5(end,:)',mu);

% orbital elements
a6 = zeros(size(t6)); % [km] semi-major axis
e6 = zeros(size(t6)); % [-] eccentricity
omega6 = zeros(size(t6)); % [rad] argument of periapsis
OMEGA6 = zeros(size(t6)); % [rad] longitude of periapsis
i6 = zeros(size(t6)); % [rad] inclination
M6 = zeros(size(t6)); % [rad] mean anomaly

for n = 1:length(t6)

    [a6(n),e6(n),omega6(n),OMEGA6(n),i6(n),M6(n)] = rv2coes(r6(n,:)',v6(n,:)',mu);

end % for

%% Everything

t = [t1; t2; t3; t4; t5; t6];
r = [r1; r2; r3; r4; r5; r6];
v = [v1; v2; v3; v4; v5; v6];
h = cross(r,v);


%% Trajectory Plot

% trajectory
figure(1)
plot3(r0(:,1),r0(:,2),r0(:,3),'k', ...
    r1(:,1),r1(:,2),r1(:,3),'g-.', ...
    r2(:,1),r2(:,2),r2(:,3),'g', ...
    r3(:,1),r3(:,2),r3(:,3),'r-.', ...
    r4(:,1),r4(:,2),r4(:,3),'r', ...
    r5(:,1),r5(:,2),r5(:,3),'b-.', ...
    r6(:,1),r6(:,2),r6(:,3),'b', ...
    0,0,0,'ko')
grid on

xlim([-5e8 5e8])
ylim([-5e8 5e8])
zlim([-5e8 5e8])

xlabel("X [km]")
ylabel("Y [km]")
zlabel("Z [km]")
legend(["Initial Orbit" "Spiral Outwards" "After Spiral" "Cranking"...
    "After Cranking" "Spiral Inwards" "Final Orbit"],'Location','best')


%% Orbital Elements Plot

figure(2)

% semi-major axis
subplot(2,3,1)
plot(t0/86400,a0/AU,'k', ...
    t1/86400,a1/AU,'b-.', ...
    t2/86400,a2/AU,'b', ...
    t3/86400,a3/AU,'r-.', ...
    t4/86400,a4/AU,'r')
grid on
title("Semi-major Axis")
xlabel("Time [days]")
ylabel("Semi-Major Axis [AU]")

% eccentricity
subplot(2,3,2)
plot(t0/86400,e0,'k', ...
    t1/86400,e1,'b-.', ...
    t2/86400,e2,'b', ...
    t3/86400,e3,'r-.', ...
    t4/86400,e4,'r')
grid on
title("Eccentricity")
xlabel("Time [days]")
ylabel("Eccentricity [-]")

% argument of periapsis
subplot(2,3,3)
plot(t0/86400,rad2deg(omega0),'k', ...
    t1/86400,rad2deg(omega1),'b-.', ...
    t2/86400,rad2deg(omega2),'b', ...
    t3/86400,rad2deg(omega3),'r-.', ...
    t4/86400,rad2deg(omega4),'r')
grid on
title("Argument of Periapsis")
xlabel("Time [days]")
ylabel("Argument of Periapsis [deg]")
legend(["Initial Orbit" "Spiral" "After Spiral" "Cranking" "After Cranking"],'Location','best')

% longitude of periapsis
subplot(2,3,4)
plot(t0/86400,rad2deg(OMEGA0),'k', ...
    t1/86400,rad2deg(OMEGA1),'b-.', ...
    t2/86400,rad2deg(OMEGA2),'b', ...
    t3/86400,rad2deg(OMEGA3),'r-.', ...
    t4/86400,rad2deg(OMEGA4),'r')
grid on
title("Longitude of Periapsis")
xlabel("Time [days]")
ylabel("Longitude of Periapsis [deg]")

% inclination
subplot(2,3,5)
plot(t0/86400,rad2deg(i0),'k', ...
    t1/86400,rad2deg(i1),'b-.', ...
    t2/86400,rad2deg(i2),'b', ...
    t3/86400,rad2deg(i3),'r-.', ...
    t4/86400,rad2deg(i4),'r')
grid on
title("Inclination")
xlabel("Time [days]")
ylabel("Inclination [deg]")

% mean anomaly
subplot(2,3,6)
plot(t0/86400,rad2deg(M0),'k', ...
    t1/86400,rad2deg(M1),'b-.', ...
    t2/86400,rad2deg(M2),'b', ...
    t3/86400,rad2deg(M3),'r-.', ...
    t4/86400,rad2deg(M4),'r')
grid on
title("Mean Anomaly")
xlabel("Time [days]")
ylabel("Mean Anomaly [deg]")


%% clear

clear
clc

%% input

% initial date
t0 = juliandate("01 Jan 2022","DD mmm YYYY");

% final date
tf = juliandate("01 Jan 2023","DD mmm YYYY");

% number of samples
N = 10000; % [-]

%% compute

% time span
tspan = linspace(t0,tf,N)';

% ephemeris
[r,v] = planetEphemeris(tspan,"Sun","Earth"); % (Nx3) [km], (Nx3) [km/s]

mu = 1.32712440018e+11; % [km^3/s^2]

% initialize
a = zeros(size(tspan));
e = zeros(size(tspan));
omega = zeros(size(tspan));
OMEGA = zeros(size(tspan));
i = zeros(size(tspan));
M = zeros(size(tspan));

for n = 1:N

    [a(n),e(n),omega(n),OMEGA(n),i(n),M(n)] = svec2keps(r(n,:)',v(n,:)',mu);

end % for

%% plot

t = (tspan-2.4515445e+06)/365.25 + 2000;

figure(1)

% semi-major axis
subplot(2,3,1)
plot(t,a,'b')
grid on
xlabel("Year")
ylabel("Semi-Major Axis [km]")
title("Semi-Major Axis")

% eccentricity
subplot(2,3,2)
plot(t,e,'b')
grid on
xlabel("Year")
ylabel("Eccentricity [-]")
title("Eccentricity")

% argument of periapsis
subplot(2,3,3)
plot(t,rad2deg(omega),'b')
grid on
xlabel("Year")
ylabel("Argument of Periapsis [deg]")
title("Argument of Periapsis")

% longitude of periapsis
subplot(2,3,4)
plot(t,rad2deg(OMEGA),'b')
grid on
xlabel("Year")
ylabel("Longitude of Periapsis [deg]")
title("Longitude of Periapsis")

% inclination
subplot(2,3,5)
plot(t,rad2deg(i),'b')
grid on
xlabel("Year")
ylabel("Inclination [deg]")
title("Inclination")

% mean anomaly
subplot(2,3,6)
plot(t,rad2deg(M),'b')
grid on
xlabel("Year")
ylabel("Mean Anomaly [deg]")
title("Mean Anomaly")

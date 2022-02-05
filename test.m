%% Clear
clear
clc

%%

launch = juliandate(datetime(2028,0,0)); % [JD]

[Re_tru,Ve_tru] = planetEphemeris(launch,'Sun','Earth');

[x,y,z,vx,vy,vz] = ephemeris_state(launch,true);
Re = [x(3) y(3) z(3)]*149597870.691;
Ve = [vx(3) vy(3) vz(3)]*149597870.691;

abs(Re_tru-Re)./Re_tru*100;

Ve_tru
Ve

%%

N = 400;
start = juliandate(datetime(2000,0,0));
stop = juliandate(datetime(2050,0,0));

tspan = linspace(start,stop,N);

tic
[R,V] = planetEphemeris(tspan','Sun','Mercury');
toc

%%

figure(1)
plot3(R(:,1),R(:,2),R(:,3))

figure(2)

t = linspace(start,start + 100,N^2);
r = interp1(tspan,R,t);

plot3(r(:,1),r(:,2),r(:,3))


%%

clear
clc

tic
[R,V] = planetEphemeris(2451545.0,'Sun','Mercury','405','AU');
toc

[x,y,z,vx,vy,vz] = ephemeris(2451545.0,true);

r = [x(1) y(1) z(1)];
v = [vx(1) vy(1) vz(1)];
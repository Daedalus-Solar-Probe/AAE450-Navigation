%% Clear

clear
clc

%% Input

% standard gravitational parameter (sun)
mu = 1.3271e+11; % [km^3/s^2]

% departure body    
bd = 'Earth';

% arrival body
ba = 'Venus';

% launch window
launch_0 = juliandate(datetime(2028,0,0)); % [JD]
launch_f = juliandate(datetime(2032,0,0)); % [JD]

% time of flight (Hohman transfer time of ~260 days)
tof_0 = 25; % [days]
tof_f = 600; % [days]

% mesh size
M = 500;
N = 500;

cutoff = 10; % [km/s]

% time spans
launch = linspace(launch_0,launch_f,M);
tof = linspace(tof_0,tof_f,N);

% evaluate subset for interpolation
tspan = linspace(launch(1),launch(end)+tof(end),N)';

% get ephemeris data
tic
[rd,vd] = planetEphemeris(tspan,'Sun',bd); % [km,km/s]
[ra,va] = planetEphemeris(tspan,'Sun',ba); % [km,km/s]
toc

Rd = interp1(tspan,rd,launch,"cubic"); % [km]
Vd = interp1(tspan,vd,launch,"cubic"); % [km/s]

Ra = interp1(tspan,ra,launch' + tof,"cubic"); % [km]
Va = interp1(tspan,va,launch' + tof,"cubic"); % [km/s]

Ra = permute(Ra,[3,1,2]); % [km]
Va = permute(Va,[3,1,2]); % [km/s]

%% solve targeting problem

% delta-V
dV1 = zeros(M,N);
dV2 = zeros(M,N);

tic
for m = 1:M
    for n = 1:N

        t0 = launch(m); % [JD]
        TOF = tof(n); % [days]
        tf = t0 + TOF; % [JD]
    
        [V01,Vf1] = lambert(Rd(m,:),Ra(:,m,n)',TOF,0,mu); % [km/s]
        [V02,Vf2] = lambert(Rd(m,:),Ra(:,m,n)',-TOF,0,mu); % [km/s]

        dV1(m,n) = min( norm( V01 - Vd(m,:) ), norm( V02 - Vd(m,:) ) ); % [km/s]
        dV2(m,n) = min( norm( Va(:,m,n)' - Vf1 ), norm( Va(:,m,n)' - Vf2) ); % [km/s]
    
    end
end
toc

dV = dV1 + dV2;

%% plot departure v_inf

tmp = dV1;
tmp(tmp > cutoff) = cutoff;

[TOF,T0] = meshgrid(tof,launch);

T0 = (T0 - T0(1,1))/365.25 + 2028;

figure(1)
contourf(T0,TOF,tmp,'ShowText','on')
grid on
colormap(flipud(winter))
h = colorbar;

title(sprintf("%s-%s Departure V-inf Porkchop Plot",bd,ba),'FontSize',16);
xlabel("Launch Date [years]",'FontSize',14);
ylabel("Time of Flight [days]",'FontSize',14);

h.Label.String = "Delta-V [km/s]";
h.Label.FontSize = 14;
h.Label.Rotation = 270;
h.Label.Position(1) = 3;

%% plot arrival v_inf

tmp = dV2;
tmp(tmp > cutoff) = cutoff;

figure(2)
contourf(T0,TOF,tmp,'ShowText','on')
grid on
colormap(flipud(winter))
h = colorbar;

title(sprintf("%s-%s Arrival V-inf Porkchop Plot",bd,ba),'FontSize',16);
xlabel("Launch Date [years]",'FontSize',14);
ylabel("Time of Flight [days]",'FontSize',14);

h.Label.String = "Delta-V [km/s]";
h.Label.FontSize = 14;
h.Label.Rotation = 270;
h.Label.Position(1) = 3;

%% plot total delta-v

tmp = dV;
tmp(tmp > 2*cutoff) = 2*cutoff;

figure(3)
contourf(T0,TOF,tmp,'ShowText','on')
grid on
colormap(flipud(winter))
h = colorbar;

title(sprintf("%s-%s Total Delta-V Porkchop Plot",bd,ba),'FontSize',16);
xlabel("Launch Date [years]",'FontSize',14);
ylabel("Time of Flight [days]",'FontSize',14);

h.Label.String = "Delta-V [km/s]";
h.Label.FontSize = 14;
h.Label.Rotation = 270;
h.Label.Position(1) = 3;


%% max arrival v-inf for min departure v-inf

tmp = dV2./dV1;
tmp(tmp < 1) = 1;

figure(4)
contourf(T0,TOF,tmp,'ShowText','on')
grid on
colormap("winter")
h = colorbar;
%% Clear

clear
clc

%% Inputs

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1-Mercury 2-Venus 3-Earth 4-Mars 5-Jupiter 6-Saturn 7-Uranus 8-Neptune %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bd = 3; % [1] departure body   
ba = 2; % [1] arrival body

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Launch Window: dd-mmm-yyyy HH:MM:SS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lw_start = '01-Jan-2028 00:00:00'; % [UTC] launch window start
lw_stop = '31-Dec-2032 23:59:59'; % [UTC] launch window stop

%%%%%%%%%%%%%%%%%%
% Time of Flight %
%%%%%%%%%%%%%%%%%%

tof_min = 50; % [days] minimum time of flight
tof_max = 300; % [days] maximum time of flight

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mesh size, M: Launch Window Cells, N: Time of Flight Cells %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M = 500; % [1] x-axis (launch window) cells
N = 500; % [1] y-axis (time of flight) cells

%%%%%%%%%%%%%%%%
% V_inf Cutoff %
%%%%%%%%%%%%%%%%

cutoff = 20; % [km/s]

%% Constants

% standard gravitational parameter (sun)
mu = 1.3271e+11; % [km^3/s^2]

% conversion between au and km
AU = 149597870.691; % [km/au]

% conversion between days and seconds
day = 86400; % [s/day]

%% Temporal mesh generation (lmao)

% time spans
launch = linspace(juliandate(lw_start,'dd-mmm-yyyy HH:MM:SS'), ...
    juliandate(lw_stop,'dd-mmm-yyyy HH:MM:SS'),M); % [JD]

tof = linspace(tof_min,tof_max,N); % [days]

% departure time mesh
t0 = repmat(launch',1,N); % [JD]

% arrival time mesh
tf = launch' + tof; % [JD]

%%

tic
% departure body
[Rd,Vd] = ephemeris(t0(:)',true); % [au,au/day]
Rd = reshape(Rd(:,:,bd),3,M,N)*AU; % [km]
Vd = reshape(Vd(:,:,bd),3,M,N)*AU/day; % [km/s]
toc

tic
% arrival body
[Ra,Va] = ephemeris(tf(:)',true); % [au,au/day]
Ra = reshape(Ra(:,:,ba),3,M,N)*AU; % [km]
Va = reshape(Va(:,:,ba),3,M,N)*AU/day; % [km/s]
toc

%% solve targeting problem

% delta-V
dV1 = zeros(M,N);
dV2 = zeros(M,N);

vv = zeros(3,M,N);

tic
for m = 1:M
    for n = 1:N

        [V01,Vf1] = lambert(Rd(:,m,n)',Ra(:,m,n)',tof(n),0,mu); % [km/s]
        [V02,Vf2] = lambert(Rd(:,m,n)',Ra(:,m,n)',-tof(n),0,mu); % [km/s]

        dV1(m,n) = min(norm(V01-Vd(:,m,n)'),norm(V02-Vd(:,m,n)')); % [km/s]
        dV2(m,n) = min(norm(Va(:,m,n)'-Vf1),norm(Va(:,m,n)'-Vf2)); % [km/s]

        vv(:,m,n) = V02-Vd(:,m,n)';
    
    end
end
toc

% total delta-V
dV = dV1 + dV2; % [km/s]

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

title("Departure V-inf Porkchop Plot",'FontSize',16);
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

title("Arrival V-inf Porkchop Plot",'FontSize',16);
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

title("Total Delta-V Porkchop Plot",'FontSize',16);
xlabel("Launch Date [years]",'FontSize',14);
ylabel("Time of Flight [days]",'FontSize',14);

h.Label.String = "Delta-V [km/s]";
h.Label.FontSize = 14;
h.Label.Rotation = 270;
h.Label.Position(1) = 3;
function [R,V] = ephemeris(T_eph,useICRF)
% EPHEMERIS Provides cartesian positions of all 8 planets relative to the
% sun in either J2000 Ecliptic or ICRF reference frames. Implemented to be
% vectorized as much as possible and therefore very fast as compared to the
% MATLAB Aerospace Dataset function "planetEphemeris" but is less accurate.
%
% Inputs:
%   T_eph: row vector containing ephemeris times to be evaluated at in
%       Julian Days (see MATLAB Aerosapce toolbox function "juliandate"
%
%   useICRF: boolean controlling the output frame of the data. By default
%       it returns everything in the J2000 Ecliptic system but when useICRF
%       is true, it uses the J2000 ICRF.
%
% Outputs:
%   R: Matrix containing the positions of each planet at the desired
%       ephemeris times. Order is ascending from sun outwards so Mercury,
%       Venus, Earth... The shape of the vector is [3xNx8] where the first
%       dimension contains the position vector, the second is the index of
%       the ephemeris, and the last in the planetary index. Unites are in
%       AU
%
%   V: Same as R but with units of AU/day
%
% Information:
%   Author: Matthew Mader
%   Contact: maderm@purdue.edu
%   Date: Feb 5, 2022
%
% Revision History:
%   Rev: IR
%   Date: Feb 2, 2022
%   Notes: Initial release, conducted limited tests compared against to the
%       "planetEphemeris" function in ICRF for inner planets.
%
%   Rev: A
%   Date: Feb 5, 2022
%   Notes: Added planetary velocities and changed the shape of the output
%       variables. Velocities slightly diverge in the z-axis of the ICRF
%       from limited testing on the inner planets.
%
% Data Source:
%   https://ssd.jpl.nasa.gov/planets/approx_pos.html

if nargin < 2
    useICRF = false;
end

% objects
objects = [
    "Mercury"
    "Venus"
    "EM Bary"
    "Mars"
    "Jupiter"
    "Saturn"
    "Uranus"
    "Neptune"
    ];

% semi-major axis [au]
a0 = [
    0.38709927
    0.72333566
    1.00000261
    1.52371034
    5.20288700
    9.53667594
    19.18916464
    30.06992276
    ];

% semi-major axis rate [au/Cy]
a_dot = [
    0.00000037
    0.00000390
    0.00000562
    0.00001847
   -0.00011607
   -0.00125060
   -0.00196176
    0.00026291
    ];

% eccentricity [1]
e0 = [
    0.20563593
    0.00677672
    0.01671123
    0.09339410
    0.04838624
    0.05386179
    0.04725744
    0.00859048
    ];

% eccentricity rate [1/Cy]
e_dot = [
     0.00001906
    -0.00004107
    -0.00004392
     0.00007882
    -0.00013253
    -0.00050991
    -0.00004397
     0.00005105
    ];

% inclination [rad]
I0 = deg2rad([
     7.00497902
     3.39467605
    -0.00001531
     1.84969142
     1.30439695
     2.48599187
     0.77263783
     1.77004347
    ]);

% inclination rate [rad/Cy]
I_dot = deg2rad([
    -0.00594749
    -0.00078890
    -0.01294668
    -0.00813131
    -0.00183714
     0.00193609
    -0.00242939
     0.00035372
    ]);

% mean longitude [rad]
L0 = deg2rad([
    252.25032350
    181.97909950
    100.46457166
     -4.55343205
     34.39644051
     49.95424423
    313.23810451
    -55.12002969
    ]);

% mean longitude rate [rad/Cy]
L_dot = deg2rad([
    149472.67411175
     58517.81538729
     35999.37244981
     19140.30268499
      3034.74612775
      1222.49362201
       428.48202785
       218.45945325
    ]);

% longitude of perihelion [rad]
long_peri0 = deg2rad([
     77.45779628
    131.60246718
    102.93768193
    -23.94362959
     14.72847983
     92.59887831
    170.95427630
     44.96476227
    ]);

% longitude of perihelion rate [rad/Cy]
long_peri_dot = deg2rad([
     0.16047689
     0.00268329
     0.32327364
     0.44441088
     0.21252668
    -0.41897216
     0.40805281
    -0.32241464
    ]);

% longitude of the ascending node [rad]
long_node0 = deg2rad([
     48.33076593
     76.67984255
      0.0
     49.55953891
    100.47390909
    113.66242448
     74.01692503
    131.78422574
    ]);

% longitude of the ascending node rate [rad/Cy]
long_node_dot = deg2rad([
    -0.12534081
    -0.27769418
     0.0
    -0.29257343
     0.20469106
    -0.28867794
     0.04240589
    -0.00508664
    ]);

% additional terms for mean anomaly
b = deg2rad([
     0.0
     0.0
     0.0
     0.0
    -0.00012452
     0.00025899
     0.00058331
    -0.00041348
    ]);

c = deg2rad([
     0.0
     0.0
     0.0
     0.0
     0.06064060
    -0.13434469
    -0.97731848
     0.68346318
    ]);

s = deg2rad([
     0.0
     0.0
     0.0
     0.0
    -0.35635438
     0.87320147
     0.17689245
    -0.10162547
    ]);

f = deg2rad([
    0.0
    0.0
    0.0
    0.0
    38.35125000
    38.35125000
    7.67025000
    7.67025000
    ]);

% centuries past J2000.0
T = (T_eph - 2451545.0)/36525; % [cY]


% current semi-major axis
a = a0 + a_dot*T; % [au]

% current eccentricity
e = e0 + e_dot*T; % [1]

% current incination
I = I0 + I_dot*T; % [rad]

% current mean longitude
L = L0 + L_dot*T; % [rad]

% current longitude of perihelion
omega_bar = long_peri0 + long_peri_dot*T; % [rad]

% current longitude of the ascending node
OMEGA = long_node0 + long_node_dot*T; % [rad]


% argument of perihelion
omega = omega_bar - OMEGA; % [rad]

% mean anomaly
M = L - omega_bar + b*T.^2 + c.*cos(f*T) + s.*sin(f*T); % [rad]

% mean anomaly can only be, -pi <= M <= pi
M = mod(M+pi,2*pi) - pi; % [rad]


% Newton-Raphson to find eccentric anomaly
E = M; % [rad] initial guesses

% iterate until convergence is reached
tol = 1e-6;
for i = 1:size(objects,1)
    for j = 1:length(T_eph)
        while abs(E(i,j)-M(i,j)-e(i,j)*sin(E(i,j))) > tol
            E(i,j) = E(i,j)-(E(i,j)-M(i,j)-e(i,j)*sin(E(i,j)))/(1-e(i,j)*cos(E(i,j))); % [rad]
        end
    end
end

% heliocentric coordinates in orbit plane
x_orbit = a.*(cos(E)-e); % [au]
y_orbit = a.*sqrt(1-e.^2).*sin(E); % [au]


% heliocentric coodinates in J2000 Ecliptic Plane
x = (cos(omega).*cos(OMEGA)-sin(omega).*sin(OMEGA).*cos(I)).*x_orbit + ...
    (-sin(omega).*cos(OMEGA)-cos(omega).*sin(OMEGA).*cos(I)).*y_orbit;

y = (cos(omega).*sin(OMEGA)+sin(omega).*cos(OMEGA).*cos(I)).*x_orbit + ...
    (-sin(omega).*sin(OMEGA)+cos(omega).*cos(OMEGA).*cos(I)).*y_orbit;

z = sin(omega).*sin(I).*x_orbit + cos(omega).*sin(I).*y_orbit;

% velocity
mu = 2.959e-04; % [au^3/day^2]
Rc = a.*(1-e.*cos(E)); % [au]

vx_orbit = -sqrt(mu*a)./Rc.*sin(E); % [au/day]
vy_orbit = sqrt(mu*a)./Rc.*sqrt(1-e.^2).*cos(E); % [au/day]

vx = (cos(omega).*cos(OMEGA)-sin(omega).*sin(OMEGA).*cos(I)).*vx_orbit + ...
    (-sin(omega).*cos(OMEGA)-cos(omega).*sin(OMEGA).*cos(I)).*vy_orbit;

vy = (cos(omega).*sin(OMEGA)+sin(omega).*cos(OMEGA).*cos(I)).*vx_orbit + ...
    (-sin(omega).*sin(OMEGA)+cos(omega).*cos(OMEGA).*cos(I)).*vy_orbit;

vz = sin(omega).*sin(I).*vx_orbit + cos(omega).*sin(I).*vy_orbit;

% if user want the coordinates in J2000 ICRF
if useICRF
    epsilon = deg2rad(23.43928); % [deg] obliquity at J2000
    y = cos(epsilon)*y - sin(epsilon)*z; % [au]
    z = sin(epsilon)*y + cos(epsilon)*z; % [au]
    
    vy = cos(epsilon)*vy - sin(epsilon)*vz; % [au/day]
    vz = sin(epsilon)*vy + cos(epsilon)*vz; % [au/day]
end

% combine x,y,z vectors
R = [x; y; z]; % [au]
V = [vx; vy; vz]; % [au/day]

% restructure the data into proper vectors
R = reshape(R,8,3,[]); % [au]
V = reshape(V,8,3,[]); % [au/day]

% reorder the indices
R = permute(R,[2,3,1]); % [au]
V = permute(V,[2,3,1]); % [au/day]
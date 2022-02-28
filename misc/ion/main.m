%% Clear

clear
clc

%% Spacecraft Definition

sc.Isp = 4200; % [s] specific impulse
sc.g0 = 9.81; % [m/s^2] standard gravity
sc.m_dry = 150; % [kg] dry mass
sc.m_wet = 1000; % [kg] wet mass
sc.F = 0.6; % [N] thrust
sc.mu = 1.327e11; % [km^3/s^2]
sc.a = 1.5*1.496e+8; % [km]

%%

y0 = [1.496e+8; 0; 0; 0; sqrt(sc.mu/1.496e+8); 0; sc.m_wet];


%% Spiral ODE Function
function dydt = spiralode(~,y,sc)

% initialize output
dydt = zeros(size(y));

% unpack
r = y(1:3); % [km]
v = y(2:4); % [km/s]
m = y(5); % [kg]

% semi-major axis
a = 1/(2/norm(r)-norm(v)^2/mu); % [km]

% gravitational acceleration
a_grav = -sc.mu*r/norm(r)^3; % [km/s^2]

% propulsion direction
if a < sc.a
    n_prop = v/norm(v); % [-]
else
    n_prop = -v/norm(v); % [-]
end % if/else

% propulsion acceleration
a_prop = sc.F/m/1000*n_prop; % [km/s^2]

% mass flow
m_dot = -sc.F/(sc.Isp*sc.g0); % [kg/s]

% derivatives
dydt(1:3) = v; % [km/s] r_dot = v
dydt(4:6) = a_grav + a_prop; % [km/s^2] v_dot = a
dydt(5) = m_dot; % [kg/s] m_dot = m_dot

end % function spiralode
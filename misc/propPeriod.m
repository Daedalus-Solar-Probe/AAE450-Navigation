function [t,r,v] = propPeriod(t0,r0,v0,mu)

% input validation
arguments

    t0 (1,1) {mustBeNumeric}
    r0 (3,1) {mustBeNumeric}
    v0 (3,1) {mustBeNumeric}
    mu (1,1) {mustBeNumeric}

end % arguments

% initial conditions
y0 = [r0; v0];

% semi-major axis
a = rv2coes(r0,v0,mu); % [km]

% orbital period
T = 2*pi*sqrt(a^3/mu); % [s]

% simulation time span
tspan = t0 + [0 T]; % [s]

% simulation options
opts = odeset("RelTol",1e-9,"AbsTol",1e-9);

% simulation 
[t,y] = ode89(@(t,y)odefun(t,y,mu),tspan,y0,opts);

% position
r = y(:,1:3);

% velocity
v = y(:,4:6);

end % function

%% ODE Function
function dydt = odefun(~,y,mu)

% initialize output
dydt = zeros(size(y));

% position
r = y(1:3); % [km]

% velocity
v = y(4:6); % [km/s]

% acceleration due to gravity
a = -mu*r/norm(r)^3; % [km/s^2]

% derivatives
dydt(1:3) = v; % [km/s] r_dot = v
dydt(4:6) = a; % [km/s^2] v_dot = a

end % function
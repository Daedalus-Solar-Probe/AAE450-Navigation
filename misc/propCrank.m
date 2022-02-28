function [t,r,v] = propCrank(t0,r0,v0,P)
%

% initial conditions
y0 = [r0; v0];

% simulation time span
tspan = t0 + [0 1e9]; % [s]

% simulation options
opts = odeset("RelTol",1e-6,"AbsTol",1e-6,"Events",@(t,y)eventsfun(t,y,P));

% simulation 
[t,y] = ode89(@(t,y)odefun(t,y,P),tspan,y0,opts);

% position
r = y(:,1:3);

% velocity
v = y(:,4:6);

end % function

%% ODE Function
function dydt = odefun(~,y,P)

% initialize output
dydt = zeros(size(y));

% unpack values
mu = P.mu; % [km^3/s^2]
accel = P.accel; % [km/s^2]

% position
r = y(1:3); % [km]

% velocity
v = y(4:6); % [km/s]

% specific angular momentum
h = cross(r,v); % [km^2/s]

% acceleration due to gravity
a_grav = -mu*r/norm(r)^3; % [km/s^2]

[~,~,~,~,~,M] = rv2coes(r,v,mu);

% propulsion acceleration
if abs(M) > pi/2
    a_prop = accel*h/norm(h); % [km/s^2]
else
    a_prop = -accel*h/norm(h); % [km/s^2]
end % if/else

% derivatives
dydt(1:3) = v; % [km/s] r_dot = v
dydt(4:6) = a_grav + a_prop; % [km/s^2] v_dot = a

end % function

%% ODE Events Function
function [val,isterm,dir] = eventsfun(~,y,P)

r = y(1:3);
v = y(4:6);
h = cross(r,v);

val = P.inc - acos(h(3)/norm(h));
isterm = 1;
dir = [];

end % function
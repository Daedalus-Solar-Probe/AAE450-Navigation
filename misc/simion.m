function [t,r,v] = simion(r0,v0,mu,a0,a1,a2,di,accel)


[t1,r1,v1] = propCrank(t0,r0,v0,P);

end % function

%% Crank ODE
function dydt = crankode(~,y,mu,accel)

% initialize output
dydt = zeros(size(y));

% position
r = y(1:3); % [km]

% velocity
v = y(4:6); % [km/s]

% specific angular momentum
h = cross(r,v); % [km^2/s]

% acceleration due to gravity
a_grav = -mu*r/norm(r)^3; % [km/s^2]

% propulsion acceleration
a_prop = accel*h/norm(h); % [km/s^2]

% derivatives
dydt(1:3) = v; % [km/s] r_dot = v
dydt(4:6) = a_grav + a_prop; % [km/s^2] v_dot = a

end % function
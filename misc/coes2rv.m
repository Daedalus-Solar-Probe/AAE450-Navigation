function [r,v] = coes2rv(a,e,omega,OMEGA,i,M0,dt,mu)
% SVEC2KEPS Converts classical orbital elements into three dimensional 
% position and velocity vectors. Units are in terms of length L and time T.
%
% Inputs:
%   a: semi-major axis (1x1) [L]
%   e: eccentricity (1x1) [-]
%   omega: argument of periapsis (1x1) [rad]
%   OMEGA: longitude of periapsis (1x1) [rad]
%   i: inclination (1x1) [rad]
%   M0: mean anomaly at t0 (1x1) [rad]
%   dt: change in time from t0 (1x1) [JD]
%   mu: standard gravitational parameter of central body (1x1) [L^3/T^2]
%
% Ouputs:
%   r: position vector (3x1) [L]
%   v: velocity vector (3x1) [L/T]
%
% Information:
%   Author: Matthew Mader
%   Contact: maderm@purdue.edu
%   Date: 27 Feb 2022
%
% Notes:
%   https://downloads.rene-schwarz.com/download/M001-Keplerian_Orbit_Elements_to_Cartesian_State_Vectors.pdf

% input validation
arguments

    a (1,1) {mustBeNumeric}
    e (1,1) {mustBeNumeric}
    omega (1,1) {mustBeNumeric}
    OMEGA (1,1) {mustBeNumeric}
    i (1,1) {mustBeNumeric}
    M0 (1,1) {mustBeNumeric}
    dt (1,1) {mustBeNumeric}
    mu (1,1) {mustBeNumeric}

end % arguments

% current mean anomaly
M = M0 + dt*86400*sqrt(mu/a^3); % (1x1) [rad]
M = mod(M,2*pi); % (1x1) [rad] on domain [0,2*pi)

% solve Kepler's equation for eccentric anomaly
E = M; % (1x1) [rad] initial guess
E_prev = inf; % (1x1) [rad] ensure it runs at least once

epsilon = 1e-6; % (1x1) [rad] absolute convergence criteria

% Newton-Raphson iterations
while abs(E-E_prev) > epsilon
    E_prev = E; % (1x1) [rad]
    E = E_prev - (E_prev-e*sin(E_prev)-M)/(1-e*cos(E_prev)); % (1x1) [rad]
end % while

% true anomaly
nu = 2*atan2(sqrt(1+e)*sin(E/2),sqrt(1-e)*cos(E/2)); % (1x1) [rad]

% distance to central body
r_c = a*(1-e*cos(E)); % (1x1) [L]

% position in the orbital frame
r_o = r_c*[cos(nu); sin(nu); 0]; % (3x1) [L]

% velocity in the orbital frame
v_o = sqrt(mu*a)/r_c*[-sin(E); sqrt(1-e^2)*cos(E); 0]; % (3x1) [L/T]

% DCM from orbital to intertial frame
C = angle2dcm(-omega,-i,-OMEGA,"ZXZ"); % (3x3) [-]

% position in the inertial frame
r = C*r_o; % (3x1) [L]

% velocity in the inertial frame
v = C*v_o; % (3x1) [L/T]

end % function
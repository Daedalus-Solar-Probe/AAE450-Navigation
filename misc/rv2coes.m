function [a,e,omega,OMEGA,i,M] = rv2coes(r,v,mu)
% SVEC2KEPS Converts three dimensional position and velocity vectors to
% classical orbital elements. Units are in terms of length L and time T.
%
% Inputs:
%   r: position vector (3x1) [L]
%   v: velocity vector (3x1) [L/T]
%   mu: standard gravitational parameter of central body (1x1) [L^3/T^2]
%
% Ouputs:
%   a: semi-major axis (1x1) [L]
%   e: eccentricity (1x1) [-]
%   omega: argument of periapsis (1x1) [rad]
%   OMEGA: longitude of periapsis (1x1) [rad]
%   i: inclination (1x1) [rad]
%   M: mean anomaly (1x1) [rad]
%
% Information:
%   Author: Matthew Mader
%   Contact: maderm@purdue.edu
%   Date: 27 Feb 2022
%
% Notes:
%   https://downloads.rene-schwarz.com/download/M002-Cartesian_State_Vectors_to_Keplerian_Orbit_Elements.pdf

% input validation
arguments

    r (3,1) {mustBeNumeric}
    v (3,1) {mustBeNumeric}
    mu (1,1) {mustBeNumeric}

end % arguments

% orbital momentum vector
h = cross(r,v); % (3x1) [km^2/s]

% eccentricity vector
e_vec = cross(v,h)/mu - r/norm(r); % (3x1) [-]

% prevent errors with circular orbits
if norm(e_vec) < 1e-9
    e_vec = zeros(3,1);
end

% orbit eccentricity
e = norm(e_vec); % (1x1) [-]

% vector pointing towards ascending node
n = cross([0,0,1]',h); % (3x1) [km^2/s]

% true anomaly
nu = acos(dot(e_vec,r)/e/norm(r)); % (1x1) [rad]
if dot(r,v) < 0
    nu = 2*pi - nu; % (1x1) [rad]
end % if

% inclination
i = acos(h(3)/norm(h)); % (1x1) [rad]

% eccentric anomaly
E = 2*atan(tan(nu/2)/sqrt((1+e)/(1-e))); % (1x1) [-]

% longitude of the ascending node
OMEGA = acos(n(1)/norm(n)); % (1x1) [rad]
if n(2) < 0
    OMEGA = 2*pi - OMEGA; % (1x1) [rad]
end % if

% argument of periapsis
omega = acos(dot(n,e_vec)/norm(n)/e); % (1x1) [rad]
if e_vec(3) < 0
    omega = 2*pi - omega; % (1x1) [rad]
end % if

% mean anomaly
M = E - e*sin(E); % (1x1) [rad]

% semi-major axis
a = 1/(2/norm(r)-norm(v)^2/mu); % (1x1) [km]

end % function
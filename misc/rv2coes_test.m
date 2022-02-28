function tests = rv2coes_test

tests = functiontests(localfunctions);

end % function svec2keps_test


function test_earth(testCase)

% maximum allowable relative error (not inclusive)
epsilon = 0.0001; % [-] (0.1% error)

% Earth position vector
r = [-1.374082700569961
    0.507824530304229
    0.220150442320484]*1e+08; % [km]

% Earth velocity vector
v = [-11.626000517565023
    -25.460431197730347
    -11.036129072339655]; % [km/s]

% Solar standard gravitational parameter
mu = 1.32712440018e+11; % [km^3/s^2]


% calculated values
[a_calc,e_calc,omega_calc,OMEGA_calc,i_calc,M_calc] = svec2keps(r,v,mu);


% expected Earth semi-major axis
a_exp = 1.496932452160402e+08; % [km]

% expected Earth eccentricity
e_exp = 1.729521265157181e-02; % [-]

% expected Earth argument of periapsis
omega_exp = 1.818853222633770; % [rad]

% expected Earth longitude of ascending node
OMEGA_exp = 3.994572288058343e-05; % [rad]

% expected Earth inclination
i_exp = 0.409030126003880; % [rad]

% expected Earth mean anomaly
M_exp = 0.912062481257041; % [rad]


% verify semi-major axis
verifyTrue(testCase, abs(a_exp-a_calc)/a_exp < epsilon)

% verify eccentricity
verifyTrue(testCase, abs(e_exp-e_calc)/e_calc < epsilon)

% verify argument of periapsis
verifyTrue(testCase, abs(omega_exp-omega_calc)/omega_exp < epsilon)

% verify longitude of ascending node
verifyTrue(testCase, abs(OMEGA_exp-OMEGA_calc)/OMEGA_calc < epsilon)

% verify inclination
verifyTrue(testCase, abs(i_exp-i_calc)/i_exp < epsilon)

% verify mean anomaly
verifyTrue(testCase, abs(M_exp-M_calc)/M_exp < epsilon)


end % function
%% clear

clear
clc

%%

% UTC gregorian date
date = '06-Jan-2023-12'; % [dd-mmm-yyyy-HH]

% UTC JD 
JD = juliandate(date,'dd-mmm-yyyy-HH'); % [JD]

% UTC MJD
MJD = mjuliandate(date,'dd-mmm-yyyy-HH'); % [MJD]

% GMAT offset
offset = mjuliandate('05-Jan-1941-12','dd-mmm-yyyy-HH'); % [MJD]

% GMAT MJD
mMJD = JD2mMJD(JD);

JD2date(JD)

mMJD2date(mMJD)




% function MJD = mjd1(S,F)
% 
% % calculate standard modified julian date
% tmp = mjuliandate(S,F);
% 
% % GMAT offset
% offset = mjuliandate('05-Jan-1941-12','dd-mmm-yyyy-HH');
% 
% % GMAT modified julian date
% MJD = tmp - offset;
% 
% end
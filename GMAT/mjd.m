function MJD = mjd(S,F)

% calculate standard modified julian date
tmp = mjuliandate(S,F);

% GMAT offset
offset = mjuliandate('05-Jan-1941-12','dd-mmm-yyyy-HH');

% GMAT modified julian date
MJD = tmp - offset;

end
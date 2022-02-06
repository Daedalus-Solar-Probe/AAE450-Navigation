function mMJD = date2mMJD(date,F)
% MMJD2DATE converts a datetime string following format F if provided to a
%   Modified Julian Date (with an offset at 05 Jan 1941 12:00:00.000) that
%   GMAT uses (UTCModJulian). The default formatting is
%   dd mmm yyyy HH:MM:SS.FFF, like 01 Jan 2022 12:00:00.000
%
% Inputs:
%   date - datetime string (using formatting F if available)
%
% Outputs:
%   mMJD - Modified Julian Date (as used by GMAT's UTCModJulian)
%
% Information:
%   Author: Matthew Mader
%   Contact: maderm@purdue.edu
%   Date: 6 Feb. 2022
%
% Revision History:
%   Rev: IR
%   Date: 6 Feb. 2022
%   Notes: Initial release, tested on dates between 2000 and 2050.

% default formatting
if nargin < 2
    F = 'dd mmm yyyy HH:MM:SS.FFF';
end

% convert to Julian Date
JD = juliandate(date,F);

% convert to Modified Julian Date
mMJD = JD2mMJD(JD);

end
function date = mMJD2date(mMJD,F)
% MMJD2DATE converts a Modified Julian Date (with an offset at 05 Jan 1941
%   12:00:00.000) that GMAT uses (UTCModJulian) to a datetime string
%   following format F if provided. Default formatting is
%   dd mmm yyyy HH:MM:SS.FFF, like 01 Jan 2022 12:00:00.000
%
% Inputs:
%   mMJD - Modified Julian Date (as used by GMAT's UTCModJulian)
%
% Outputs:
%   date - datetime string (using formatting F if available)
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

% convert Modified Julian Date to Julian Date
JD = mMJD2JD(mMJD);

% conversion
date = JD2date(JD,F);

end
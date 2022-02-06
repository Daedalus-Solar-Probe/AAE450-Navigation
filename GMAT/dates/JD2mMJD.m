function mMJD = JD2mMJD(JD)
% JD2MMJD Converts a Julian Date (from MATLAB's juliandate() function) into
%   the Modified Julian Date (with an offset at 05 Jan 1941 12:00:00.000)
%   that GMAT uses (UTCModJulian). Think of it as modified-modified julian.
%
% Inputs:
%   JD - Julian Date (as output from MATLAB's juliandate() function)
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

% data source:
% http://gmat.sourceforge.net/docs/nightly/html/SpacecraftEpoch.html

% Modified Julian Date offset
offset = 2430000.0;

% conversion
mMJD = JD - offset;

end
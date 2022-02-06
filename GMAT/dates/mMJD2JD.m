function JD = mMJD2JD(mMJD)
% MMJD2JD Converts the Modified Julian Date (with an offset at 05 Jan 1941
%   12:00:00.000) that GMAT uses (UTCModJulian) to a Julian Date (from
%   MATLAB's juliandate() function).
%
% Inputs:
%   mMJD - Modified Julian Date (as used by GMAT's UTCModJulian)
%
% Outputs:
%   JD - Julian Date (as output from MATLAB's juliandate() function)
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
JD = mMJD + offset;

end
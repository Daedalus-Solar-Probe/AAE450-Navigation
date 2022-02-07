function date = JD2date(JD,F)
% JD2DATE converts a Julian Date (from MATLAB's juliandate() function) to a
%   datetime string following format F if provided. Default formatting is
%   dd mmm yyyy HH:MM:SS.FFF, like 01 Jan 2022 12:00:00.000
%
% Inputs:
%   JD - Julian Date (as output from MATLAB's juliandate() function)
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

% corresponds to the juliandate at the start of the common era
offset = 1721058.5;

% conversion
date = datestr(JD-offset,F);

end
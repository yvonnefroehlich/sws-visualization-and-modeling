function vers_out = SWS_modeling_check_matlab_version()

% ==========================================================================
%% This function
% ==========================================================================
% checks the MATLAB version available on your system
% --------------------------------------------------------------------------
% is 
% - modified: Yvonne Fröhlich (ORCID 0000-0002-8566-0619)
%   https://github.com/yvonnefroehlich/sws-visualization-and-modeling
%   Fröhlich (2025) Dissertation
%   https://doi.org/10.5445/IR/1000183786
%   Fröhlich, Grund, Ritter (2024) Geophysical Journal International
%   https://doi.org/10.1093/gji/ggae245
%   Ritter, Fröhlich, Sanz Alonso, Grund (2022) Journal of Seismology
%   https://doi.org/10.1007/s10950-022-10112-w
% - created: Michael Grund (ORCID 0000-0001-8759-2018)
%   https://github.com/michaelgrund/sws_tools
%   Grund, Ritter (2020) Geophysical Journal International
%   https://doi.org/10.1093/gji/ggaa388
%   Grund (2019) Dissertation
%   https://doi.org/10.5445/IR/1000091425
% --------------------------------------------------------------------------
% LICENSE
%
% Copyright (C) 2026  Yvonne Fröhlich, Michael Grund (v2.0)
% Copyright (C) 2022  Yvonne Fröhlich, Michael Grund (v1.0)
% https://github.com/yvonnefroehlich/sws-visualization-and-modeling
% Copyright (C) 2020  Michael Grund (sws_tools)
% https://github.com/michaelgrund/sws_tools
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program. If not, see <http://www.gnu.org/licenses/>.
% --------------------------------------------------------------------------
% TERMS OF USE
%
% The modeling routines are provided "as is" and without any warranty.
% The author cannot be held responsible for anything that happens to you
% or your equipment. Use it at your own risk.
% --------------------------------------------------------------------------
% CONTRIBUTING
%
% Feel free to modify/adjust the code for your needs. Submit improvements
% and report bugs by opening a "New issue" in the GitHub repository (:
% ==========================================================================


vers = version('-release');

vers_yyyy = str2double(vers(1:4));

% MATLAB R2020a or higher
if vers_yyyy>2019
    vers_out=1;
else
    vers_out=0;
end


%==========================================================================
end % EOF

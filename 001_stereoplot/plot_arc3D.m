function P = plot_arc3D(a, b, h, k, r, colfill, coledge, facealphafill)

% ==========================================================================
%% This function
% ==========================================================================
% plots a circular arc as a pie wedge
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
% - extended: Michael Grund (ORCID 0000-0001-8759-2018)
%   https://github.com/michaelgrund/sws_tools
%   Grund, Ritter (2020) Geophysical Journal International
%   https://doi.org/10.1093/gji/ggaa388
%   Grund (2019) Dissertation
%   https://doi.org/10.5445/IR/1000091425
% - based on >>> plot_arc.m <<< by Matt Fig
%   https://de.mathworks.com/matlabcentral/answers/6322-drawing-a-segment-of-a-circle
%   (last access 2022 June 19)
% --------------------------------------------------------------------------
% INPUT
%
%   a                start of arc in radians
%   b                end of arc in radians
%   h,k              center
%   r                radius
%   colfill          fill color
%   coledge          outline color
%   facealphafill    transparency of fill color
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
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program. If not, see <http://www.gnu.org/licenses/>.
% --------------------------------------------------------------------------
% TERMS OF USE
%
% The plotting routines are provided "as is" and without any warranty.
% The author cannot be held responsible for anything that happens to you
% or your equipment. Use it at your own risk.
% --------------------------------------------------------------------------
% CONTRIBUTING
%
% Feel free to modify/adjust the code for your needs. Submit improvements
% and report bugs by opening a "New issue" in the GitHub repository (:
% ==========================================================================



t = linspace(a,b);
x = r*cos(t) + h;
y = r*sin(t) + k;

x = [x h x(1)];
y = [y k y(1)];

P = fill3(x, y, ones(length(x))*1.01, 'r');
% fill(x,y,'r');
set(P, 'facecolor',colfill, 'edgecolor',coledge, 'facealpha',facealphafill)
% axis([h-r-1 h+r+1 k-r-1 k+r+1])
% axis tight;

if ~nargout
    clear P
end


%==========================================================================
end % EOF

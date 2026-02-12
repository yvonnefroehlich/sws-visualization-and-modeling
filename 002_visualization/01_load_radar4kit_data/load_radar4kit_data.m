% ==========================================================================
% Load CSV files with shear wave splitting measurements (SWSMs)
% --------------------------------------------------------------------------
% Datasets are available from RADAR4KIT
% - Upper Rhine Graben area: https://dx.doi.org/10.35097/685
%   related to Fröhlich et al. (2024)
% - Black Forest Observatory: https://dx.doi.org/10.35097/684
%   related to Ritter et al. (2022)
% --------------------------------------------------------------------------
% Related to the publications
% - Fröhlich (2025) Dissertation
%   https://doi.org/10.5445/IR/1000183786
% - Fröhlich, Grund, Ritter (2024) Geophysical Journal International
%   https://doi.org/10.1093/gji/ggae245
% - Ritter, Fröhlich, Sanz Alonso, Grund (2022) Journal of Seismology
%   https://doi.org/10.1007/s10950-022-10112-w
% --------------------------------------------------------------------------
% Created: 2024/07/08
% Author: Yvonne Fröhlich
% ORCID: https://orcid.org/0000-0002-8566-0619
% --------------------------------------------------------------------------
% LICENSE
%
% Copyright (C) 2026  Yvonne Fröhlich (v2.0)
% https://github.com/yvonnefroehlich/sws-visualization-and-modeling
% https://doi.org/10.5281/zenodo.7213156
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
% The loading routines are provided "as is" and without any warranty.
% The author cannot be held responsible for anything that happens to you
% or your equipment. Use it at your own risk.
% --------------------------------------------------------------------------
% CONTRIBUTING
%
% Feel free to modify/adjust the code for your needs. Submit improvements
% and report bugs by opening a "New issue" in the GitHub repository (:
% ==========================================================================



clear all

% Example
root_path = '';
header = 15;  % Number of header lines, given in README
station = 'BFO';
network = 'GR';
obstyp = 'NULLS';

% Load CSV file into table
table_swsm = readtable(...
    [root_path 'splitresults_' obstyp '_goodfair_' network '_' station '.csv'], ...
    "Delimiter",';', ...
    "NumHeaderLines",header ...
);

% Convert column "year_jday" from floating point number to string
table_swsm.year_jday = num2str(table_swsm.year_jday, '%.3f');

% Split column "year_jday" into two new columns "year" and "jday"
table_swsm.year = table_swsm.year_jday(:,1:4);
table_swsm.jday = table_swsm.year_jday(:,6:8);

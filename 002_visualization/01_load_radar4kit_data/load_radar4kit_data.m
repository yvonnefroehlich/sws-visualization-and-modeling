% #############################################################################
% Load CSV files with shear wave splitting measurements (SWSMs)
% -----------------------------------------------------------------------------
% Datasets are available from RADAR4KIT
% - Upper Rhine Graben area: https://dx.doi.org/10.35097/685
%   related to Fröhlich et al. (2024)
% - Blackforest Observatory: https://dx.doi.org/10.35097/684
%   related to Ritter et al. (2022)
% -----------------------------------------------------------------------------
% Related to the publications
% - Fröhlich Y., Grund M. & Ritter J. R. R. (2024).
%   Lateral and vertical variations of seismic anisotropy in the
%   lithosphere-asthenosphere system underneath Central Europe from
%   long-term splitting measurements.
%   Geophysical Journal International, 239(1), 112-135.
%   https://doi.org/10.1093/gji/ggae245.
% - Ritter J. R. R., Fröhlich Y., Sanz Alonso Y. & Grund M. (2022).
%   Short-scale laterally varying SK(K)S shear wave splitting at BFO,
%   Germany – implications for the determination of anisotropic structures.
%   Journal of Seismology, 26, 1137-1156.
%   https://doi.org/10.1007/s10950-022-10112-w.
% -----------------------------------------------------------------------------
% Created: 2024/07/08
% Author: Yvonne Fröhlich
% ORCID: https://orcid.org/0000-0002-8566-0619
% #############################################################################


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

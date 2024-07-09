% #########################################################################
% Load CSV files with shear wave splitting measurements (SWSMs)
% -----------------------------------------------------------------------------
% Data is available from RADAR4KIT
% - Upper Rhine Graben area: https://dx.doi.org/10.35097/685
%   related to Fröhlich et al. (2024)
% - Blackforest Observatory: https://dx.doi.org/10.35097/684
%   related to Ritter et al. (2022)
% -----------------------------------------------------------------------------
% Related to the publications
% - Fröhlich Y., Grund M., Ritter J. R. R. (2024).
%   Lateral and vertical variations of seismic anisotropy in the
%   lithosphere-asthenosphere system underneath Central Europe from
%   long-term splitting measurements. Geophysical Journal International,
%   accepted June 27 2024.
% - Ritter J. R. R., Fröhlich Y., Sanz Alonso Y. & Grund M. (2022).
%   Short-scale laterally varying SK(K)S shear wave splitting at BFO,
%   Germany – implications for the determination of anisotropic structures.
%   Journal of Seismology, 26, 1137-1156.
%   https://doi.org/10.1007/s10950-022-10112-w.
% -------------------------------------------------------------------------
% created: 2024/07/08
% author: Yvonne Fröhlich
% contact: yvonne.froehlich@kit.edu
% #########################################################################

clear all

% Example
station = 'BFO';
network = 'GR';
obstyp = 'NULLS';

% Load CSV file into table
table_swsm = readtable([ ...
    'splitresults_' obstyp '_goodfair_' network '_' station '.csv' ...
]);

% Convert colum  "year_jday" from floating point number to string
table_swsm.year_jday = num2str(table_swsm.year_jday, '%.3f');

% Split column "year_jday" into two new columns "year" and "jday"
table_swsm.year = table_swsm.year_jday(:,1:4);
table_swsm.jday = table_swsm.year_jday(:,6:8);

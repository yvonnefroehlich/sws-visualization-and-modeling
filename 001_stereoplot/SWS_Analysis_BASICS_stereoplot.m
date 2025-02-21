function SWS_Analysis_BASICS_stereoplot(colmap)
% function SWS_Analysis_BASICS_stereoplot(colmap, SL_quality, SL_phase, SL_obs)


% ==========================================================================
%% This function
% ==========================================================================
% reads (single seismological recording station related)
% - single-event analysis (SplitLab, SL) and
% - multi-event analysis (StackSplit) result files, and
% prepares and saves stereoplots as publication-ready pdf, png, eps
% --------------------------------------------------------------------------
% uses the provided MATLAB functions
% - SWS_Analysis_BASICS_check_matlab_version.m
% - SWS_Analysis_BASICS_read_SLresults.m
% - SWS_Analysis_BASICS_read_SSresults.m
% - plot_arc3D.m
% >>> plot_arc3D.m <<< is based on >>> plot_arc.m <<< by Matt Fig
% https://de.mathworks.com/matlabcentral/answers/6322-drawing-a-segment-of-a-circle
% (last access 2022 June 19)
% --------------------------------------------------------------------------
% is
% - based on: >>> stereoplot.m <<< function of SplitLab
%   Wüstefeld et al. (2008) Computers & Geosciences
%   https://doi.org/10.1016/j.cageo.2007.08.002
% - created: Michael Grund (ORCID 0000-0001-8759-2018)
%   https://github.com/michaelgrund/sws_tools
%   Grund (2019) Dissertation
%   https://doi.org/10.5445/IR/1000091425
%   Grund & Ritter (2020) Geophysical Journal International
%   https://doi.org/10.1093/gji/ggaa388
% - modified: Yvonne Fröhlich (ORCID 0000-0002-8566-0619)
%   https://github.com/yvonnefroehlich/sws-visualization-and-modeling
%   Ritter, Fröhlich, Sanz Alonso & Grund (2022) Journal of Seismology
%   https://doi.org/10.1007/s10950-022-10112-w
%   Fröhlich, Grund & Ritter (2024) Geophysical Journal International
%   https://doi.org/10.1093/gji/ggae245
% --------------------------------------------------------------------------
% LICENSE
%
% Copyright (C) 2022  Yvonne Fröhlich & Michael Grund (up on v1.0)
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



%==========================================================================
%% How to use
%==========================================================================
% After processing the whole data of a station with SL (and StackSplit)
%
% 1) Change to the folder with the output files of SL and StackSplit
%
%    required: splitresults_*.txt, splitresultsNULL_*.txt
%    optional: *stackresults.mat
%
% 2) Run this function SWS_Analysis_BASICS_stereoplot()
%
%    If result lists are available, they are loaded and processed
%    completely automatically. Just follow the occurring queries.
%
% If bars should be color-coded with respect to the fast polarization
% direction (phi), pass the colormap of your choice as input, e. g.:
%
%    SWS_Analysis_BASICS_stereoplot('lajolla')
%    SWS_Analysis_BASICS_stereoplot('viridis')
%    SWS_Analysis_BASICS_stereoplot('winter')
%
% If no argument is given, the default colormap parula (flipped) is used:
%
%    SWS_Analysis_BASICS_stereoplot() or SWS_Analysis_BASICS_stereoplot
%
% For plotting without color-coding use:
%
%    SWS_Analysis_BASICS_stereoplot('none')
%--------------------------------------------------------------------------
% Supported colormaps
%
% -> please download them and add the corresponding dictionary to your
% MATLAB search path before running this function
%
% 1) MATLAB colormaps (built-in): parula, winter, summer, copper, ...
%
% 2) MatPlotLib Perceptually Uniform Colormaps
%    - MATLAB: v2.1.3 https://de.mathworks.com/matlabcentral/fileexchange/62729-matplotlib-perceptually-uniform-colormaps
%      (last access 2022 June 26)
%
% 3) Scientific colour maps. F. Crameri (2021) Zenodo.
%    http://doi.org/10.5281/zenodo.1243862
%    http://www.fabiocrameri.ch/colourmaps.php
%    - MATLAB: https://de.mathworks.com/matlabcentral/fileexchange/68546-crameri-perceptually-uniform-scientific-colormaps
%      (last access 2023 April 10)
%      Please note a bug in "crameri.m" of v1.08:
%      Line 97 standardizes all colormap names to be lower-case. As MATLAB
%      is a case-sensitive programming language colormaps containing
%      upper-case letters are not found in the provided MATLAB struct.
%      This should be fixed in v1.09.
%
% 4) cmocean colormaps. Thyng et al. (2016) Oceanography 29(3):9–13.
%    http://dx.doi.org/10.5670/oceanog.2016.66
%    - MATLAB: v2.02 https://de.mathworks.com/matlabcentral/fileexchange/57773-cmocean-perceptually-uniform-colormaps
%      (last access 2022 June 18)
%--------------------------------------------------------------------------
% Automatically appearing queries
%
% -> you can select during running this function
%
% - measuremnt quality (see SWS_Analysis_BASICS_read_SLresults.m)
%    * all
%    * good
%    * good & fair
%    * fair & poor
%    * poor
% - seismological phase (see SWS_Analysis_BASICS_read_SLresults.m)
%    * all
%    * SKS
%    * SKKS
%    * PKS
% - observation type (see SWS_Analysis_BASICS_read_SLresults.m)
%    * all
%    * nulls (no shear wave splitting)
%    * splits (shear wave splitting)
% - shear wave splitting method
%    * RC (rotation-correlation method, Bowman & Ando 1987)
%    * SC (energy minimization method, Silver & Chan 1991)
%    * EV (eigenvalue method, Silver & Chan 1991)
% - multi-event analysis results of StackSplit (if available)
%    * STACK (stacking of error surfaces, Wolfe & Silver 1998)
%    * SIMW (simultaneous inversion of multiple waveforms, Roy et al. 2007)
% - backazimuth sector in white, rest in gray (time intense)
%    * no
%    * [lower_BAZ_limit, upper_BAZ_limit]
% - position of radial axis annotation (incidence angle 5, 10, 15 deg)
%    * no
%    * NE, SE, SW, NW
%--------------------------------------------------------------------------
% Changeable settings (mainly for plotting)
%
% -> please set them before running this function
%
% - plotting of legends and annotations
% - fill color of backazimuth sector or rest as well as background
% - general appearance of symbols
%==========================================================================



% check MATLAB version
vers = SWS_Analysis_BASICS_check_matlab_version;


%==========================================================================
%% Changebale settings
%==========================================================================
% >>> adjust for your needs <<<

%--------------------------------------------------------------------------
% What annotation should be plotted?
status_cb = 'yes'; %% 'yes', 'no'  % colorbar - phi color-coding of bars
status_leg = 'yes'; %% 'yes', 'no'  % legend - null, delay time reference
status_sta_label = 'yes'; %% 'yes', 'no'  % station name - station code
status_sta_marker = 'no'; %% 'yes', 'no'  % station symbol - inverse triangle
status_baz = 'yes'; %% 'yes', 'no'  % angle axis - BAZ - N(orth), E(ast)
filename_add = '';  %% additional string added to the filename


%--------------------------------------------------------------------------
% plot sector
% default is white between 0 and 360 degrees

% selected BAZ sector
white_value = 255.99999999;
% rest
colfill = [219,219,219]./256; % light gray
%colfill = [190,190,190]./256; % darker gray

% Examples:
% 1) BAZ sector between 20 and 120 deg in white, rest is gray
%   - pass the vector when the corresponding query appears [20,120]
%   - set above before running this function colfill = [190,190,190]./256;
% 2) shade full background in light gray
%   - pass the vector when the corresponding query appears [0,0.001]
%   - set above before running this function colfill = [220,220,220]./256;

%--------------------------------------------------------------------------
% symbols
linew = 2.5; % thickness of bars
marks = 7;
linewcirc = 2;
fontsize_baz = 15;

%--------------------------------------------------------------------------
% recording station
color_sta_label = [255 90 0] ./ 256; % orange [0.6350 0.0780 0.1840]; % dark red
color_sta_marker = [255 215 0] ./ 256;

%--------------------------------------------------------------------------
% radial axis & legend & colorbar
col_inc = 'k';
col_leg = 'k';
fontsize_leg = 10;
fontsize_cb = 9;

%--------------------------------------------------------------------------
% no phi color-coding (default colors of MATLAB)
splitcol = [0.3 0.3 0.3]; % dark grey
nullcol = [0.8500 0.3250 0.0980]; % red-orange
multicol_stack = [0 0.4470 0.7410]; % dark blue
multicol_simw = [0.3010 0.7450 0.9330]; % light blue

color_SKS = [216.75 82.875 24.990] ./ 256; % SKS
color_SKKS = [236.895 176.97 31.875] ./ 256; % SKKS
color_PKS = [238 238 0] ./ 2556; % PKS


%==========================================================================
%% NOT changeable settings
%==========================================================================

% horizontal position of legends
startval = 0.254;

% manually adjusted to fit the length to the bars plotted via plotm
lengthbar = 0.0710;



%==========================================================================
%% Colormaps
%==========================================================================


%==========================================================================
% list of supported colormaps

%--------------------------------------------------------------------------
% Scientific colour maps
% - without diverging and multi-sequential colormaps
% - cyclic colormaps bamoO, brocO, corkO, vikO, romaO
crameri_cmap = {
    'batlow', ...
    'devon','lajolla','bamako', ...
    'davos','bilbao','nuuk', ...
    'oslo','grayC','hawaii', ...
    'lapaz','tokyo','buda', ...
    'acton','turku','imola', ...
    'bamO','brocO','corkO','vikO','romaO', ...
};

%--------------------------------------------------------------------------
% MatPlotLib colormaps
mpl_cmap = {'viridis','magma','inferno','plasma'};

%--------------------------------------------------------------------------
% cmocean colormaps
% - without diverging and multi-sequential colormaps
% - cyclic colormap phase
cmocean_cmap = {
    'thermal','haline','solar','ice','gray','deep','dense', ...
    'algae','matter','turbid','speed','amp','tempo','rain', ...
    'phase', ...
};


%==========================================================================
% define colormap

if ~exist('colmap','var') || strcmp(colmap,'parula')
    usecmap = parula(181);
    usecmap = flipud(usecmap);
    colmap = 'parulaflip';
    fast_col = 1;
elseif strcmp(colmap,'none')
    fast_col = 0;
else
    fast_col = 1;

%--------------------------------------------------------------------------
    % check for colormaps on your system
    %......................................................................
    % Scientific colour maps
    if vers~=0  % MATLAB R2016b or higher
        idxpre1 = contains(crameri_cmap,colmap);
        idx1 = sum(  double( strcmp(crameri_cmap(idxpre1),colmap) )  );
    else
        checkcmaps = strfind(crameri_cmap,colmap);
        idx1 = find(~contains('isempty',checkcmaps));
        if isempty(idx1) || length(idx1)>1 || ...
                (isscalar(idx1) && ~strcmp({colmap},crameri_cmap(1)))
           idx1 = 0;
        else
           idx1 = 1;
        end
    end
    %......................................................................
    % MatPlotLib colormaps
    if vers~=0  % MATLAB R2016b or higher
        idxpre2 = contains(mpl_cmap,colmap);
        idx2 = sum(double(strcmp(mpl_cmap(idxpre2),colmap)));
    else
        checkcmaps2 = strfind(mpl_cmap,colmap);
        idx2 = find(~contains('isempty',checkcmaps2));
        if isempty(idx2) || length(idx2)>1 || ...
                (isscalar(idx2) && ~strcmp({colmap},mpl_cmap(1)))
           idx2 = 0;
        else
           idx2 = 1;
        end
    end
    %......................................................................
    % cmocean colormaps
    if vers~=0  % MATLAB R2016b or higher
        idxpre3 = contains(cmocean_cmap,colmap);
        idx3 = sum(double(strcmp(cmocean_cmap(idxpre3),colmap)));
    else
        checkcmaps3 = strfind(cmocean_cmap,colmap);
        idx3 = find(~contains('isempty',checkcmaps3));
        if isempty(idx3) || length(idx3)>1 || ...
                (isscalar(idx3) && strcmp({colmap},cmocean_cmap(1)))
           idx3 = 0;
        else
           idx3 = 1;
        end
    end

%--------------------------------------------------------------------------
    % search for input colormap
    %......................................................................
    % Scientific colour maps
    if idx1==1 && ~isempty(which('crameri.m'))
        usecmap = crameri(colmap,181);
        disp(' ')
        disp('>>> Scientific colour maps found! <<<')
    elseif idx1==1 && isempty(which('crameri.m'))
        warning('Scientific colour maps not found!')
        return
    %......................................................................
    % MatPlotLib colormaps
    elseif idx2==1 && ~isempty(which(colmap))
        usecmap = colormap([colmap '(181)']);
        disp(' ')
        disp('>>> MatPlotLib Colormaps found! <<<')
    elseif idx2==1 && isempty(which(colmap))
        warning('MatPlotLib Colormaps not found!')
        return
    %......................................................................
    % cmocean colormaps
    elseif idx3==1 && ~isempty(which('cmocean.m'))
        usecmap = cmocean(colmap,181);
        disp(' ')
        disp('>>> cmocean colormaps found! <<<')
    elseif idx3==1 && isempty(which('cmocean.m'))
        warning('cmocean colormaps not found!')
        return
    %......................................................................
    % built-in MATLAB colormaps (use: help colormap)
    elseif  idx1==0 && idx2==0 && idx3==0
        if exist(colmap,'file')
            usecmap = colormap([colmap '(181)']);
        else
            error('Colormap not available!')
        end

    end

end

close all



%==========================================================================
%% Read SL single-event-analysis results
%==========================================================================

%--------------------------------------------------------------------------
% Make queries for
% - measurement quality: 0; 1; 2; 3; 4; 5
% - seismological phase: 0; 1; 2; 3
% - observation type: 0; 1; 2

% give numbers directly here, then no queries occur
% [RES_split, RES_nulls, SL_qualtiy, SL_phase, SL_obs] = ...
%     SWS_Analysis_BASICS_read_SLresults(SL_quality, SL_phase, SL_obs);
[RES_split, RES_nulls, SL_qualtiy, SL_phase, SL_obs] = ...
    SWS_Analysis_BASICS_read_SLresults();

% corresponding to numbers in queries before
quality_str = {'all'; 'good'; 'goodfair'; 'fairpoor'; 'fair'; 'poor'};
phase_str = {'XKS'; 'SKS'; 'SKKS'; 'PKS'};
obs_str = {'swsms'; 'nulls'; 'splits'};

single_string = 'single';
if isempty(RES_split) && isempty(RES_nulls)
    single_string = '';
end

%--------------------------------------------------------------------------
% Make query for
% - measurement method: 1; 2; 3

disp(' ')
SL_method = input(['Methode you want to plot [Default is SC]? \n' ...
                   '   [1] SC  [2] RC  [3] EV    | ']);

if ~exist('SL_method','var')==1  % default
    SL_method = 1; % SC
end

% corresponding to numbers in query before
method_str = {'SC'; 'RC'; 'EV'};



%==========================================================================
%% Read multi-event-analysis results if available and make query
%==========================================================================

% default
RES_multi = [];
plot_multi = 0;  % no

% corresponding to number in query
multi_string = {''; 'stack'; 'SIMWNN'; 'stackSIMWNN'};

dir_res_multi = dir('*_stackresults.mat');

if ~isempty(dir_res_multi)
    disp(' ')
    plot_multi = input(['Plot multi-event-analysis results' ...
                        ' (if available)? \n' ...
                        '  [0] no  [1] stack  [2] SIMW(NN)' ...
                        '  [3] stack & SIMW(NN)   | ']);
    if plot_multi>0
        RES_multi = SWS_Analysis_BASICS_read_SSresults(...
                        dir_res_multi, 1, plot_multi);
    end

    if plot_multi==1 && isempty(RES_multi)
        error('No stack results in struct!')
    elseif plot_multi==2 && isempty(RES_multi)
        error('No simw results contained in struct!')
    end

end

if ~exist('plot_multi','var')==1  % default
    plot_multi = 0;  % no multi-event analysis results
end


%==========================================================================
%% Adjust colors for bars based on phases
%==========================================================================

if strcmp(colmap, 'none') && plot_multi == 0
    switch SL_phase
        case 1
            splitcol = color_SKS;
            nullcol = color_SKS;
        case 2
            splitcol = color_SKKS;
            nullcol = color_SKKS;
        case 3
            splitcol = color_PKS;
            nullcol = color_PKS;
    end
end


%==========================================================================
%% Check if data is from one single station
%==========================================================================

station_check = {};

if ~isempty(RES_nulls)
    for i_null=1:1:length(RES_nulls)
        station_check{end+1} = RES_nulls(i_null).staname;
    end
end

if ~isempty(RES_split)
    for i_split=1:1:length(RES_split)
        station_check{end+1} = RES_split(i_split).staname;
    end
end

if ~isempty(RES_multi)
    for i_multi=1:1:length(RES_multi)
        station_check{end+1} = RES_multi(i_multi).staname;
    end
end

station_check = unique(station_check);

% error in case SWSMs are from different stations
if length(station_check) > 1
    error('>>> Input files with data from different stations! <<<')
end

% error in case not SWSMs of the selected qualities are available at this station
if isempty(station_check)
    error(['>>> No shear wave splitting measurements are available ' ...
           'for the selected qualities at this station! <<<'])
end



%==========================================================================
%% Setup variables
%==========================================================================

% (I) single-event analysis
% (i) splits
if ~isempty(RES_split)
    bazi = [RES_split.baz];
    inc = abs([RES_split.inc]);  % problem with negative incidence angle
    switch SL_method
        case 1
            azim_pre = [RES_split.phiSC];
            len = [RES_split.dtSC];
        case 2
            azim_pre = [RES_split.phiRC];
            len = [RES_split.dtRC];
        case 3
            azim_pre = [RES_split.phiEV];
            len = [RES_split.dtEV];
    end
    azim = azim_pre;
    staname = RES_split.staname;

    % second SKKS in SL negative incidence angle
    for n = 1:1:length(RES_split)
        if RES_split(n).inc < 0
            if bazi(n) < 180
                bazi(n) = bazi(n)+180;
                azim(n) = azim(n)-180;
            else
                bazi(n) = bazi(n)-180;
                azim(n) = azim(n)+180;
            end
        end
    end

end

% (ii) null
if ~isempty(RES_nulls)
    bazi_nulls = [RES_nulls.baz];
    inc_nulls = abs([RES_nulls.inc]);
    switch SL_method
        case 1
            azim_nulls_pre = [RES_nulls.phiSC];
            len_nulls = [RES_nulls.dtSC];
        case 2
            azim_nulls_pre = [RES_nulls.phiRC];
            len_nulls = [RES_nulls.dtRC];
        case 3
            azim_nulls_pre = [RES_nulls.phiEV];
            len_nulls = [RES_nulls.dtEV];
    end
    azim_nulls = azim_nulls_pre;
    staname = RES_nulls.staname;

    for n = 1:1:length(RES_nulls)
       if RES_nulls(n).inc < 0
          if bazi_nulls(n) < 180
              bazi_nulls(n) = bazi_nulls(n)+180;
              azim_nulls(n) = azim_nulls(n)-180;
          else
              bazi_nulls(n) = bazi_nulls(n)-180;
              azim_nulls(n) = azim_nulls(n)+180;
          end
       end
    end

end

%--------------------------------------------------------------------------
% (II) multi-event analysis
if ~isempty(RES_multi)

    % calculate mean incidence
    incALLmean = nan(length(RES_multi),1);
    for ii = 1:1:length(RES_multi)
        incALL = nan(length(RES_multi(ii).used_phases),1);

        for jj = 1:1:length(RES_multi(ii).used_phases)
            incALL(jj) = RES_multi(ii).used_phases(jj).results.incline;
        end

        incALLmean(ii) = mean(incALL);
        clear incALL
    end

    bazi_multi = [RES_multi.meanbaz];
    inc_multi = abs(incALLmean);
    azim_multi_pre = [RES_multi.phimulti];
    azim_multi = azim_multi_pre;
    len_multi = [RES_multi.dtmulti];
    staname = RES_multi(1).staname;

    for n = 1:1:length(RES_multi)
       if incALLmean(n) < 0
          if bazi_multi(n) < 180
              bazi_multi(n) = bazi_multi(n)+180;
              azim_multi(n) = azim_multi(n)-180;
          else
              bazi_multi(n) = bazi_multi(n)-180;
              azim_multi(n) = azim_multi(n)+180;
          end
       end
    end

end



%==========================================================================
%% Make query for sector plotting
%==========================================================================
disp(' ')
plotsector = input(['Plot sector in backazimuth range (takes some time)? \n ' ...
                    '   Plot no sector: Press "Enter" [Default is used] \n ' ...
                    '   Plot a sector: Pass a vector, e.g., [0,210]    | ']);

if ~exist('plotsector','var')==1  % default
    plotsector = [];  % BAZ range 0°-360°
end

if ~isempty(plotsector)
    if length(plotsector)==2
        lowlim = plotsector(1);
        upplim = plotsector(2);
    else
        error(['Array length needs to be 2! ' ...
               'Please passe only values between 0 and 360 deg!'])
    end
else
    lowlim = 0;
    upplim = 360;
end



%==========================================================================
%% Make query for location of annotation of radial axis
%==========================================================================
disp(' ')
plotannot = input(['Annotate radial scale [Default is SE]? \n' ...
                   '   [0] no  [1] NE  [2] SE  [3] SW  [4] NW    | ']);

if ~exist('plotannot','var')==1  % default
    plotannot = 2;  % no
end



%==========================================================================
%% Make figure
%==========================================================================

f_stereo = figure();

%{
if ~isempty(RES_split)
    m = max(inc);
elseif ~isempty(RES_nulls)
    m = max(inc_nulls);
elseif ~isempty(RES_multi)
    m = max(incALLmean);
end

m = round(m/10)*10;
lim = [-inf m+7];
%}

lim = [-inf 17];
lim_sector = 0.2610/15*lim(2);

axesm('stereo', 'Frame','on', 'Grid','on' ,'Origin',[90 0], ...
      'MlineLocation',30, 'PlineLocation',5, 'fLatLimit',lim, ...
      'fLineWidth',1, 'GLinestyle','-', 'GLinewidth',0.4, ...
      'Gcolor',[0.8 0.8 0.8]);

axis tight
axis off
view([0 -90])

axes = gca;
axes.SortMethod = 'ChildOrder';  % for right order of layers in eps / pdf

framem('FLinewidth',2)


%==========================================================================
% annotation
L = min(abs(axis));

%--------------------------------------------------------------------------
% North (N) and East (E)
if strcmp(status_baz,'yes')
    text(0, -L-0.005, 'N', ...
         'HorizontalAlignment','Center', 'VerticalAlignment','Base', ...
         'fontsize',fontsize_baz);
    text(L+0.005, 0, 'E', ...
         'HorizontalAlignment','Left', 'verticalAlignment','middle', ...
         'fontsize',fontsize_baz);
end

%--------------------------------------------------------------------------
% station name
if strcmp(status_sta_label,'yes')
    text(-0.2, -L+0.01, staname, ...
         'HorizontalAlignment','Center', 'VerticalAlignment','Base', ...
         'FontWeight','bold', 'fontsize',18, ...
         'color',color_sta_label)
end

%--------------------------------------------------------------------------
% radial axis (inclination angle)
switch plotannot
    case 1  % NE
       text(0.047,-0.047, '5^\circ', 'fontsize',12, 'color',col_inc)
       text(0.105,-0.105, '10^\circ', 'fontsize',12, 'color',col_inc)
       text(0.168,-0.168, '15^\circ', 'fontsize',12, 'color',col_inc)
    case 2  % SE
       text(0.040,0.056, '5^\circ', 'fontsize',12, 'color',col_inc)
       text(0.105,0.105, '10^\circ', 'fontsize',12, 'color',col_inc)
       text(0.175,0.145, '15^\circ', 'fontsize',12, 'color',col_inc)
    case 3  % SW
       text(-0.055,0.055, '5^\circ', 'fontsize',12, 'color',col_inc)
       text(-0.117,0.117, '10^\circ', 'fontsize',12, 'color',col_inc)
       text(-0.180,0.180, '15^\circ', 'fontsize',12, 'color',col_inc)
    case 4  % NW
       text(-0.055,-0.055, '5^\circ', 'fontsize',12, 'color',col_inc)
       text(-0.117,-0.117, '10^\circ', 'fontsize',12, 'color',col_inc)
       text(-0.180,-0.180, '15^\circ', 'fontsize',12, 'color',col_inc)
end

%--------------------------------------------------------------------------
% plot recording station
if strcmp(status_sta_marker, 'yes')
    plot(0, 0, 'v', ...
        'MarkerSize',14, ...
        'MarkerEdgeColor','k', ...
        'MarkerFaceColor', ...
        color_sta_marker)
end


%==========================================================================
% plot sector
% >>> function < plot_arc3D.m > is required, based on plot_arc.m <<<

if exist('plot_arc3D','file')

    if lowlim<upplim
        if lowlim~=0 || upplim~=360
            % first plot whole BAZ range in gray as bottom layer
            startwedge = 0;
            endwedge = 360;
            plot_arc3D( ...
               deg2rad(startwedge-90), deg2rad(endwedge-90), ...
               0, 0, lim_sector, ...
               colfill, ...
               colfill, 1 ...
            );

            % then plot considered BAZ range again on top in white
            startwedge = lowlim;
            endwedge = upplim;
            plot_arc3D( ...
                deg2rad(startwedge-90), deg2rad(endwedge-90), ...
                0, 0, lim_sector, ...
               [white_value white_value white_value]./256, ...
               [white_value white_value white_value]./256, 1 ...
            );
        end

    % this is the case when the higher value is slightly higher than 0
    % and the other one in the SA region
    elseif lowlim>upplim
        if lowlim~=0 || upplim~=360
            % first plot whole BAZ range in gray as bottom layer
            startwedge = 0;
            endwedge = 360;
            plot_arc3D( ...
               deg2rad(startwedge-90), deg2rad(endwedge-90), ...
               0, 0, lim_sector, ...
               colfill, colfill, 1 ...
            );

            % then plot considered BAZ range again on top in white
            % in two steps
            startwedge = 0;
            endwedge = upplim;
            plot_arc3D( ...
               deg2rad(startwedge-90), deg2rad(endwedge-90), ...
               0, 0, lim_sector, ...
               [white_value white_value white_value]./256, ...
               [white_value white_value white_value]./256, 1 ...
            );

            startwedge = lowlim;
            endwedge = 360;
            plot_arc3D( ...
               deg2rad(startwedge-90), deg2rad(endwedge-90), ...
               0, 0, lim_sector, ...
               [white_value white_value white_value]./256, ...
               [white_value white_value white_value]./256, 1 ...
            );
        end
    end

end


%==========================================================================
% nulls
for KK = 1:1:length(RES_nulls)
    if ~isempty(RES_nulls) && fast_col==0
        plotm(90-inc_nulls(KK), bazi_nulls(KK), 'o', ...
              'color',nullcol, 'MarkerSize',marks, ...
              'linewidth',linewcirc, 'markerFacecolor','w');
    elseif ~isempty(RES_nulls) && fast_col==1
        plotm(90-inc_nulls(KK), bazi_nulls(KK), 'o', ...
              'color','k', 'MarkerSize',marks, ...
              'linewidth',linewcirc, 'markerFacecolor','w');
        cmap = usecmap;
        colormap(cmap);
    end
end


%==========================================================================
if ~isempty(RES_split)
    NNull = 1:length(RES_split);

    bazi = bazi(:);
    inc  = inc(:);
    len  = len(:);
    azim = azim(:);

    bazi = [bazi(NNull)  bazi(NNull)]';
    inc  = [inc(NNull)   inc(NNull)]';
    len  = [-len(NNull)  len(NNull)]';
    azim = (bazi-[azim(NNull) azim(NNull)]');

elseif ~isempty(RES_nulls)
    NNull = 1:length(RES_nulls);

    bazi = bazi_nulls(:);
    inc  = inc_nulls(:);
    len  = len_nulls(:);
    azim = azim_nulls(:);

    bazi = [bazi(NNull)  bazi(NNull)]';
    inc  = [inc(NNull)   inc(NNull)]';

    len  = [-len(NNull)  len(NNull)]';
    azim = (bazi-[azim(NNull) azim(NNull)]');
end

%--------------------------------------------------------------------------
if ~isempty(RES_multi) && plot_multi>0
    NNull = 1:length(RES_multi);

    bazi_multi = bazi_multi(:);
    inc_multi  = inc_multi(:);
    len_multi  = len_multi(:);
    azim_multi = azim_multi(:);

    bazi_multi = [bazi_multi(NNull)  bazi_multi(NNull)]';
    inc_multi  = [inc_multi(NNull)   inc_multi(NNull)]';
    len_multi  = [-len_multi(NNull)  len_multi(NNull)]';
    azim_multi = (bazi_multi-[azim_multi(NNull) azim_multi(NNull)]');
end


%==========================================================================
if ~isempty(RES_split)

    % scale marker to output size
    len = len*2; % one second==4 deg (2 deg in both directions)
    if ~isempty(RES_multi)
        len_multi = len_multi*2;
    end

    % singles
    [latout, lonout] = reckon(90-inc, bazi, len, azim, 'degrees');
    hndl = plotm(latout, lonout, 'Linewidth',linew);

    % multi
    if ~isempty(RES_multi) && plot_multi>0
        [latout_multi, lonout_multi] = reckon( ...
            90-inc_multi, bazi_multi, len_multi, azim_multi, 'degrees' ...
        );
        hndl_multi = plotm(latout_multi, lonout_multi, 'Linewidth',linew);
    end

%--------------------------------------------------------------------------
    % color-coding based on phi
    if fast_col==1

        cmap = usecmap;
        step_phi = -90:1:90;
        colormap(usecmap);

        % single
        for ii = 1:1:length(hndl)
            if ~isempty(RES_split)
                azim_rounded = floor(azim_pre(1,ii));
                if isnan(azim_pre(1,ii)) % NaNs in azim_pre producig error
                    azim_rounded = 90;
                end
            else
                azim_rounded = floor(azim_nulls_pre(1,ii));
            end

            index = step_phi==azim_rounded;
            set(hndl(ii), 'color',cmap(index,:))
            set(hndl(ii), 'linewidth',linew)
        end

        % multi
        if ~isempty(RES_multi) && plot_multi>0
            for ii = 1:1:length(hndl_multi)

                azim_rounded = floor(azim_multi_pre(1,ii));

                index = step_phi==azim_rounded;
                set(hndl_multi(ii), 'color',cmap(index,:))
                set(hndl_multi(ii), 'linewidth',linew)
            end
        end

    % no color-coding based on phi
    else
        % single
        set(hndl, 'color',splitcol, 'linewidth',linew)
        % multi
        if ~isempty(RES_multi) && plot_multi>0
            for ii = 1:1:length(RES_multi)
                if strcmp(RES_multi(ii).stack_meth,'SIMW')
                    multicol = multicol_simw;
                else
                    multicol = multicol_stack;
                end
                set(hndl_multi(ii), 'color',multicol, 'linewidth',linew)
            end
        end
    end

end


%==========================================================================
% plot legend (null, delay time reference)

if strcmp(status_leg,'yes')

%--------------------------------------------------------------------------
    % null

    if fast_col==0 && SL_phase==0
        col_leg_null = nullcol;
    elseif fast_col==1 || fast_col==0 && SL_phase~=0
        col_leg_null = col_leg;
    end

    plot([startval startval], ...
         [0.220 0.220], 'o', 'linewidth',linewcirc, 'markersize',marks, ...
         'markerfacecolor','w', 'markeredgecolor',col_leg_null)
    text(startval, 0.195, 'null', ...
         'HorizontalAlignment','center', ...
         'fontsize',fontsize_leg, 'color',col_leg)

%--------------------------------------------------------------------------
    % delay time of splits

    plot([startval-lengthbar/2 startval+lengthbar/2], ...
         [0.260 0.260], '-', 'linewidth',linew, 'color',col_leg)
    text(startval, 0.243, '1 s', ...
         'HorizontalAlignment','center', ...
         'fontsize',fontsize_leg, 'color',col_leg)

    plot([startval-lengthbar*2/2 startval+lengthbar*2/2], ...
         [0.295 0.295], '-', 'linewidth',linew, 'color',col_leg)
    text(startval, 0.278, '2 s', ...
         'HorizontalAlignment','center', ...
         'fontsize',fontsize_leg, 'color',col_leg)

end


%==========================================================================
% colorbar for phi OR legend for split single, stack split, simw split

if strcmp(status_cb,'yes')

%--------------------------------------------------------------------------
    % colorbar for phi
    % phi color-coding and for null stations

    % if ~isempty(RES_split) && fast_col==1 || isempty(RES_split) && fast_col==1
    if fast_col==1

        cb = colorbar('location','north');
        zlab = get(cb,'xlabel');
        set(zlab,'String','        \phi_a / N\circE');

        clim([-90 90])
        set(cb,'xtick',-60:30:60);
        set(cb,'fontsize',fontsize_cb)
        set(cb,'TickDirection','out')
        set(cb,'TickLength',0.0200)

        axpos = get(gca,'position');
        cbpos = get(cb,'position');

        cbpos(1) = 0.38+cbpos(1);
        cbpos(2) = 0.10+cbpos(2);
        cbpos(2) = 0.10+cbpos(2);

        cbpos(3) = 0.4*cbpos(3);
        cbpos(4) = 0.4*cbpos(4);
        set(cb,'position',cbpos)
        set(gca,'position',axpos)

        % use predefined size values
        set(cb,'position',[0.6150 0.94100 0.2200 0.0200])
        set(gca,'position',[0.1300 0.1100 0.7745 0.8150])

%--------------------------------------------------------------------------
    % legend for split single, stack split, simw split
    % no phi color-coding and multi

    elseif fast_col==0 && plot_multi~=0

        pos_line_top = [-L -L];
        pos_text_top = -L-0.02;
        pos_line_mid = [-L+0.035 -L+0.035];
        pos_text_mid = -L+0.015;
        pos_line_bot = [-L+0.07 -L+0.07];
        pos_text_bot = -L+0.05;

        % single
        plot([startval-lengthbar/2 startval+lengthbar/2], ...
             pos_line_top, '-', ...
             'linewidth',linew, 'color',splitcol)
        text(startval, pos_text_top, 'split single', ...
             'HorizontalAlignment','center', ...
             'fontsize',fontsize_leg, 'color',col_leg)

        if plot_multi~=1 % -> simw
            plot([startval-lengthbar/2 startval+lengthbar/2], ...
                 pos_line_mid, '-', ...
                 'linewidth',linew, 'color',multicol_simw)
            text(startval, pos_text_mid, 'split simw', ...
                 'HorizontalAlignment','center', ...
                 'fontsize',fontsize_leg, 'color',col_leg)
        end

        if plot_multi~=2 % -> stack
            if plot_multi==1 % only stack
                pos_line_stack = pos_line_mid;
                pos_text_stack = pos_text_mid;
            elseif plot_multi==3 % stack and simw
                pos_line_stack = pos_line_bot;
                pos_text_stack = pos_text_bot;
            end
            plot([startval-lengthbar/2 startval+lengthbar/2], ...
                 pos_line_stack, '-', ...
                 'linewidth',linew, 'color',multicol_stack)
            text(startval, pos_text_stack, 'split stack', ...
                 'HorizontalAlignment','center', ...
                 'fontsize',fontsize_leg, 'color',col_leg)
        end

    end

end



%==========================================================================
%% Save plot
%==========================================================================

%--------------------------------------------------------------------------
% predefined and optimized to give nearly identical sizes on most systems
set(f_stereo, 'PaperPosition',[-3.4655 -1.1548 20.3046 15.2284]);
set(f_stereo, 'PaperSize',[14 14]); % set paper size

% If saved figures do not fulfill your expectations you can adjust
% the size and scaling here by uncomment the following lines
% and modify the individual values:
% figsize = get(f_stereo, 'PaperPosition');
% set(f_stereo, 'PaperPosition',[figsize(1)-4.1 figsize(2)-7.5 figsize(3) figsize(4)]);

%--------------------------------------------------------------------------
file_path = [];
file_name = [ ...
    'Stereo_' staname '_' ...
     quality_str{SL_qualtiy+1} '_' ...
     method_str{SL_method} '_' ...
     phase_str{SL_phase+1} '_' ...
     obs_str{SL_obs+1} '_' ...
     single_string multi_string{plot_multi+1} '_' ...
     'BAZ' num2str(lowlim) 'to' num2str(upplim) '_' ...
     colmap ...
     filename_add ...
 ];

% MATLAB built-in function "exportgraphics" requires MATLAB 2020a+
% format svg not supported by MATLAB built-in function "exportgraphics",
% transparency not supported for format png by MATLAB built-in function
% "exportgraphics", you can try MATLAB file exchange function "export_fig"
% Yair Altman (2021). Retrieved August 25, 2021.
% MATLAB file exchange: https://de.mathworks.com/matlabcentral/fileexchange/23629-export_fig
% GitHub: https://github.com/altmany/export_fig/releases/tag/v3.15

if vers==2 % MATLAB R2020a and higher
    exportgraphics(f_stereo, [file_path file_name '.png'], 'Resolution',360)
    exportgraphics(f_stereo, [file_path file_name '.eps'], 'ContentType','vector')
    exportgraphics(f_stereo, [file_path file_name '.pdf'], 'ContentType','vector')
else
    saveas(f_stereo, [file_path file_name '.png'])
    saveas(f_stereo, [file_path file_name '.eps'])
    saveas(f_stereo, [file_path file_name '.pdf'])
end

% alternatively on Linux systems print as eps and then convert to pdf
% print ('-depsc', '-painters','-r600', [file_path file_name '.eps'])
% dir_eps_file = dir([file_path file_name '.eps']);
% [status,cmdout] = system(['epstopdf ' dir_eps_file.name]);

disp(' ')
disp(['>>> File "' file_name '" saved! :D <<<'])


%==========================================================================
end  % EOF

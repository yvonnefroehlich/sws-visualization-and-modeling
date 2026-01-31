function modsall_sort = SWS_modeling_calc_misfit( ...
    modelsin, ...
    modrange_low, modrange_upp, ...
    datasplit, datanull, ...
    datastack, datasimw, ...
    domper ...
)

% ==========================================================================
%% This function
% ==========================================================================
% fits synthetic structural anisotropy models to observed splitting
% parameters (fast polarization direction phi and delay time dt)
% - one layer with horizontal symmetry axis (H1)
% - two layers with horizontal symmetry axes (H2)
% - one layer with tilted symmetry axis (T1)
% based on the minimum root mean square error (RMSE)
% --------------------------------------------------------------------------
% uses the provided MATLAB functions
% - SWS_modeling_read_data.m
% - SWS_modeling_plot_results.m
% - SWS_modeling_plot_stereo_synthetic.m
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
% INPUT
%
% - modelsin | str
%   filename of the file which consists all re-computed synthetic models
%   (e. g. modelsin = 'sws_modout_domper8s.mat';)
%
% >>> To firstly create the synthetic models please use the function
% SWS_modeling_precomp_models_main.m <<<
%
% - modrange_low | int
%   lower backazimuth limit to model in degree (0-360)
% - modrange_upp | int
%   upper backazimuth limit to model in degree (0-360)
%
% - datasplit | str
%   filename of split data ('splitresults_*.txt')
% - datanull | str
%   filename of null data ('splitresultsNULL_*.txt')
% - datastack | str
%   filename of stack data ('splitresultsSTACK_*.txt')
% - datasimw | str
%   filename of simw data ('splitresultsSIMW_*.txt')
%
% If a file is not available, pass "" (e. g. datastack = "";)
%
% >>> currently EITHER stack OR simw can be used at a time <<<
%
% >>> all files (datasplit, datanull, datastack, datasimw) need to be in
% standard SplitLab or StackSplit output format! <<<
%
% >>> Be sure to exclude discrepant pairs (SKS-SKKS) from your dataset
% before running this function <<<
%
% - domper | int
%   dominant period in sec used in the forward calculation of the splitting
%   parameter for the synthetic anisotropy models (e. g. domper = 8;)
%
% --------------------------------------------------------------------------
% EXAMPLE using the test data set provided via the download package
%
% 1) change to directory .../testdata
%
% 2) define input variables
%{
    modelsin = 'sws_modout_domper8s.mat';
    modrange_low = 30;
    modrange_upp = 110;
    datasplit = 'splitresults_GR_BFO_TEST.txt';
    datanull = 'splitresultsNULL_GR_BFO_TEST.txt';
    datastack = 'splitresultsSTACK_GR_BFO_TEST.txt';
    datasimw = "" % 'splitresultsSIMW_GR_BFO_TEST.txt';
    domper = 8;
%}
% 3) run misfit routine
%{
    modsall_sort = SWS_modeling_calc_misfit( ...
                       modelsin, modrange_low, modrange_upp, ...
                       datasplit, datanull, datastack, datasimw, ...
                       domper );
%}
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



%==========================================================================
% loading preprocessed synthetic models
disp('Loading model file')
disp(modelsin)

% struct with field splitmods
models = load(modelsin);

% struct with fields phi_eff, dt_eff, mod_paras, type
model_out = models.splitmods;


%==========================================================================
%% INITIAL SETTINGS
%==========================================================================
% >>> adjust for your needs <<<
% mainly for later plotting etc.

% maximum XX best models for finale figures
plot_mod_max = 15;  % 20; % keep only XX best models of ALL model types
keep_mods = 100;  % 500; % keep only XX best models of ALL model types
keep_mods_sep = plot_mod_max; % keep only XX best models of EACH model type

%--------------------------------------------------------------------------
% MODEL a specific BAZ range of the data (baz1:baz2)
modrange_col = [219,219,219]./256;
modrange_edcol = modrange_col;

% PLOT a specific BAZ range
BAZ = 0:1:360;

%--------------------------------------------------------------------------
% MODELS plotting
colmod_bf_1 = [255 90 0]./256; % dark orange best (plot_mod_max = 1)
colmod_bf_2max = [175 175 175]./256; % gray rest (plot_mod_max > 1)

%--------------------------------------------------------------------------
% SPLITS plotting

fontsize = 11;
marksize = 7;
fs_RMSE = 8;
lw_symb = 1; % thickness outline of symbols
lw_mod = 1.2;

% single splits
colorsedge(1,:) = [0 0 0];
colorsfill(1,:) = [66 91 169]./256; % blue

% stack splits
colorsedge(2,:) = [0 0 0];
colorsfill(2,:) = [79 197 104]./256; % green

% simw splits
colorsedge(3,:) = [0 0 0];
colorsfill(3,:) = [0.9290 0.6940 0.1250]; % light orange

%--------------------------------------------------------------------------
% NULLS plotting
ms_null = 6.5;
lw_symb_null = lw_symb; % thickness of outline of symbols
col_edge_null = 'k';
col_face_null = 'w';

%--------------------------------------------------------------------------
% model parameter plotting
mymarkersize_symbols = 7;
mylinewidth_symbols = 1.2;



%==========================================================================
%% check MATLAB struct
%==========================================================================
if ~isfield(model_out,'phi_eff') && ~isfield(model_out,'dt_eff')
    error(['Required fields >phi_eff< & >dt_eff< do not exist \n' ...
           'in variable >model_out< (synthetic models)! \n ' ...
           'Check input struct!'])
end



%==========================================================================
%% color-coding of synthetic anisotropy models based on RMSE
%==========================================================================
% >>> based on the 20 best-fit models <<<

% >>> colormap grayC is part of the Scientific colour maps <<<
% F. Crameri (2021) Zenodo.
% https://doi.org/10.5281/zenodo.1243862
% https://www.fabiocrameri.ch/colourmaps.php
%    - MATLAB: v1.08 https://de.mathworks.com/matlabcentral/fileexchange/68546-crameri-perceptually-uniform-scientific-colormaps
%    (last access 2022 June 25)
%    >>> v1.08: bug in >crameri.m<. To fix it uncomment line 97 of this
%    function. This line standardizes all colormap names to lower-case.
%    As MATLAB is a case-sensitive programming language colormaps
%    containing upper-case letters are not found in the provided MATLAB
%    struct anymore. <<<

%--------------------------------------------------------------------------
% make query use phi color-coding and which colormap
disp(' ')
cmap_rms_ind = input(['Color-coding based on RMSE [Default is no]?: \n' ...
                      '    [0] no  [1] gray  [2] grayC    | ']);

% warning if colormap name contains upper-case letters
if cmap_rms_ind==2 % grayC
    disp(' ')
    warning(['v1.08: bug in >crameri.m<. '  ...
             'To fix it uncomment line 97 of this function. ' ...
             'This line standardizes all colormap names to ' ...
             'lower-case. As MATLAB is a case-sensitive ' ...
             'programming language, colormap names containing ' ...
             'upper-case letters are not longer found in the ' ...
             'provided MATLAB struct.'])
end

% select RMSE colormap and set string for file name
if cmap_rms_ind==0
    cmap_rms_str = 'rmsno';
elseif cmap_rms_ind==1
    % built-in MATLAB
    cmap_rms = flipud( gray(plot_mod_max + 10) );
    cmap_rms_str = 'rmsgray';
elseif cmap_rms_ind==2
    % grayC of Scientific Colormaps version 7.0.1 by Crameri 2021
    % Flip colormap by using '-' before the colormap name
    cmap_rms = crameri('-grayC', plot_mod_max + 10);
    cmap_rms_str = 'rmsgrayC';
end

if cmap_rms_ind~=0
    cmap_rms_sel = cmap_rms(9:plot_mod_max+8,:);
else
    cmap_rms_sel = "";
end

%--------------------------------------------------------------------------
% default values
if isempty(cmap_rms_ind)
    cmap_rms_sel = "";
    cmap_rms_ind = 0; % no
    cmap_rms_str = 'rmsno';
end



%==========================================================================
%% color-coding of symbols in backazimuthal plot based on phi
%==========================================================================
% >>> colormap phase it part of the cmocean colormaps <<<
% Thyng et al. (2016) Oceanography 29(3):9–13.
% https://dx.doi.org/10.5670/oceanog.2016.66
% - MATLAB: v2.02 https://de.mathworks.com/matlabcentral/fileexchange/57773-cmocean-perceptually-uniform-colormaps
%   (last access 2022 June 18)

%--------------------------------------------------------------------------
% make query use phi color-coding and which colormap
disp(' ')
cmap_phi_ind = input(['Color-coding based on phi [Default is no]?: \n' ...
                      '    [0] no  [1] parula(flipped)  [2] phase    | ']);

% select phi colormap and set string for file name
if cmap_phi_ind==0
    cmap_phi = "";
    cmap_phi_str = 'phino';
elseif cmap_phi_ind==1
    % built-in MATLAB
    cmap_phi = flipud( parula(181) );
    cmap_phi_str = 'phiparulaflip';
elseif cmap_phi_ind==2
    % phase of cmocean colormaps by Thyng et al. 2016
    cmap_phi = cmocean('phase', 181);
    cmap_phi_str = 'phiphasemap';
end

%--------------------------------------------------------------------------
% make query plot phi colorbar
disp(' ')
cbar_phi_ind = input(['Add colorbar for phi [Default is no]?: \n' ...
                      '    [0] no  [1] yes    | ']);

% set phi colorbar string for file name
if cmap_phi_ind~=0 & cbar_phi_ind==1
    cbar_phi_str = '_cb_';
else
    cbar_phi_str = '_';
end

%--------------------------------------------------------------------------
% default values
if isempty(cmap_phi_ind)
    cmap_phi_ind = 0; % no
    cmap_phi = "";
    cmap_phi_str = 'phino';
    cbar_phi_ind = 0; % no
    cbar_phi_str = '_';
end



%==========================================================================
%% fitting method
%==========================================================================

%--------------------------------------------------------------------------
% make query which fitting method

disp(' ')
whichfit = input(['Fitting method [Default is joint]?: \n' ...
                  '    [1] only phi  [2] joint phi & dt   | ']);

%--------------------------------------------------------------------------
% default value
if isempty(whichfit)
    whichfit = 2; % joint
end



%==========================================================================
%% observed data input
%==========================================================================

%--------------------------------------------------------------------------
% check that EITHER stack OR simw are provided (NOT both at a time)
if strcmp(datastack,"")==0 && strcmp(datasimw,"")==0
    error("Currently EITHER stack OR simw data can be used!")
    return
end

%--------------------------------------------------------------------------
% read SplitLab and StackSplit results
% >>> all txt files must be in the standard output format <<<

dir_res_split = dir(datasplit);
dir_res_null = dir(datanull);
dir_res_stack = dir(datastack);
dir_res_simw = dir(datasimw);

% use only good & fair qualities
% no query from function >>> SWS_modeling_read_data <<< appears
use_QUAL = 2; % adjust if you want other qualities

[RES_split, RES_null, RES_stack, RES_simw_split, RES_simw_null] = ...
    SWS_modeling_read_data( dir_res_split, dir_res_null, ...
                            dir_res_stack, dir_res_simw, ...
                            use_QUAL );

staname_split = RES_split(1).staname;



%==========================================================================
%% select phases that should be used
%==========================================================================

phaselist_split = {RES_split.phase};
phaselist_null = {RES_null.phase};
phaselist_all = unique( horzcat(phaselist_split, phaselist_null) );

%--------------------------------------------------------------------------
% make query which phases available
disp(' ')
disp('Available phases: ')

for ii = 1:1:length(phaselist_all)+1
    if ii==1
        disp('    [0] ALL')
    else
        disp(['    [' num2str(ii-1) '] only ' phaselist_all{ii-1}])
    end
end

%--------------------------------------------------------------------------
% make query which phases to use
wphases = input(['Which phases to use ' ...
                 '(pass a number from the list)?    | ']);

%--------------------------------------------------------------------------
% default value
if isempty(wphases)
    wphases = 0; % all available phases
elseif ~isempty(wphases) && ~ismember(wphases,0:1:length(phaselist_all))
    error('>>> You selected a non-available entry. Try again!')
end

if wphases~=0
    phases2use=phaselist_all(wphases);
    find_index_split = strcmp(phaselist_split,phases2use);
    find_index_null = strcmp(phaselist_null,phases2use);

    RES_split = RES_split(find_index_split);
    RES_null = RES_null(find_index_null);
end



%==========================================================================
%% prepare variables
%==========================================================================


%==========================================================================
% (I) single-event analysis
data_used_str = 'singleSC';

%--------------------------------------------------------------------------
% (i) splits
% create vector with [lower err phi; phi; upper err phi]
meas_phiSC = [RES_split.phiSC_err_min; ...
              RES_split.phiSC; ...
              RES_split.phiSC_err_max]';
% 1 in 4th column means single splits for later coloring
meas_phiSC(:,4) = 1;
% create vector with [lower err dt; dt; upper err dt]
meas_dtSC = [RES_split.dtSC_err_min; ...
             RES_split.dtSC; ...
             RES_split.dtSC_err_max]';

meas_BAZ = [RES_split.baz]';
% >>> round corresponding BAZs of measured splitting parameters
% (no big difference), % since theoretical values are only available
% as integers, otherwise no comparison possible <<<
meas_BAZ_floor = floor(meas_BAZ);

%--------------------------------------------------------------------------
% (ii) nulls
meas_phiSC_null = [RES_null.phiSC_err_min; ...
                   RES_null.phiSC; ...
                   RES_null.phiSC_err_max]';
meas_dtSC_null = [RES_null.dtSC_err_min; ...
                  RES_null.dtSC; ...
                  RES_null.dtSC_err_max]';
meas_BAZ_null = [RES_null.baz]';
meas_BAZ_floor_null = floor(meas_BAZ_null);


%==========================================================================
% (II) multi-event analysis
% >>> under development <<<
% >>> at the moment EITHER simw OR stack can be used <<<

%--------------------------------------------------------------------------
% (i) simw
% simw splits
if ~isempty(RES_simw_split)

    data_used_str = 'multiSIMW';

    meas_phiSC_simw_split = [RES_simw_split.phiSC_err_min; ...
                             RES_simw_split.phiSC; ...
                             RES_simw_split.phiSC_err_max]';
    % 3 in 4th column means stacked splits for later coloring
    meas_phiSC_simw_split(:,4) = 3;
    meas_dtSC_simw_split = [RES_simw_split.dtSC_err_min; ...
                            RES_simw_split.dtSC; ...
                            RES_simw_split.dtSC_err_max]';
    meas_BAZ_simw_split = [RES_simw_split.meanbaz]';
    meas_BAZ_floor_simw_split = floor(meas_BAZ_simw_split);

    % merge single splits with simw splits
    merged_phiSC = vertcat(meas_phiSC, meas_phiSC_simw_split);
    merged_dtSC = vertcat(meas_dtSC, meas_dtSC_simw_split);
    merged_BAZ_floor = vertcat(meas_BAZ_floor, meas_BAZ_floor_simw_split);

    % rename original variables
    meas_phiSC = merged_phiSC;
    meas_dtSC = merged_dtSC;
    meas_BAZ_floor = merged_BAZ_floor;

    % for plotting only
    meas_phiSC4plot = merged_phiSC;
    meas_dtSC4plot = merged_dtSC;
    meas_BAZ_floor4plot = merged_BAZ_floor;
else
    meas_phiSC4plot = meas_phiSC;
    meas_dtSC4plot = meas_dtSC;
    meas_BAZ_floor4plot = meas_BAZ_floor;
end

% simw nulls
if ~isempty(RES_simw_null)
    meas_phiSC_simw_null = [RES_simw_null.phiSC_err_min; ...
                            RES_simw_null.phiSC; ...
                            RES_simw_null.phiSC_err_max]';
    meas_dtSC_simw_null = [RES_simw_null.dtSC_err_min; ...
                            RES_simw_null.dtSC; ...
                            RES_simw_null.dtSC_err_max]';
    meas_BAZ_simw_null = [RES_simw_null.meanbaz]';
    meas_BAZ_floor_simw_null = floor(meas_BAZ_simw_null);
else
    meas_phiSC_simw_null = [];
    meas_dtSC_simw_null = [];
    meas_BAZ_floor_simw_null = [];
end

% >>> at the moment EITHER simw OR stack can be used <<<

%--------------------------------------------------------------------------
% (ii) stack
if ~isempty(RES_stack)

    data_used_str = 'multiWS';

    meas_phiSC_stack = [RES_stack.phiSTACK_err_min; ...
                        RES_stack.phiSTACK; ...
                        RES_stack.phiSTACK_err_max]';
    % 2 in 4th column means stacked splits for later coloring
    meas_phiSC_stack(:,4) = 2;
    meas_dtSC_stack = [RES_stack.dtSTACK_err_min; ...
                       RES_stack.dtSTACK; ...
                       RES_stack.dtSTACK_err_max]';
    meas_BAZ_stack = [RES_stack.meanbaz]';
    meas_BAZ_floor_stack = floor(meas_BAZ_stack);

    %......................................................................
    % no double counting of for stack results used single results
    kk = 1;
    index = 0;
    for ii = 1:1:length(meas_BAZ)
        baz = 0;
        while baz<359
            s = 0;
            % s=0 -> no stack in this BAZ interval
            % s=1 -> stack in this BAZ interval
            for jj = 1:1:length(meas_BAZ_stack)
                if (meas_BAZ_stack(jj) >= baz) && ...
                    (meas_BAZ_stack(jj) < baz+5)
                    s = 1;
                    break
                end
            end
            if s
                if (meas_BAZ(ii) >= baz) && (meas_BAZ(ii) < baz+5)
                    break
                end
            end
            baz = baz + 5;
        end
        if baz==360
            % if while-loop was not broken
            % -> if res_split(ii).baz is not in stacked intervals
            index(kk) = ii;
            kk = kk+1;
        end
    end

    if index==0 % when only stack results, and no single results remaining
        reduce_phiSC = [];
        reduce_dtSC = [];
        reduce_BAZ_floor = [];
    else
        reduce_phiSC = meas_phiSC(index,:);
        reduce_dtSC = meas_dtSC(index,:);
        reduce_BAZ_floor = meas_BAZ_floor(index);
    end

    %......................................................................
    % merge single splits with stacked splits
    merged_phiSC = vertcat(reduce_phiSC, meas_phiSC_stack);
    merged_dtSC = vertcat(reduce_dtSC, meas_dtSC_stack);
    merged_BAZ_floor = vertcat(reduce_BAZ_floor, meas_BAZ_floor_stack);

    %......................................................................
    % UNCOMMENT if you want to use all single spits and all stack splits

%     merged_phiSC = vertcat(meas_phiSC, meas_phiSC_stack);
%     merged_dtSC = vertcat(meas_dtSC, meas_dtSC_stack);
%     merged_BAZ_floor = vertcat(meas_BAZ_floor, meas_BAZ_floor_stack);

    %......................................................................
    % rename original variables
    meas_phiSC = merged_phiSC;
    meas_dtSC = merged_dtSC;
    meas_BAZ_floor = merged_BAZ_floor;

    % for plotting only
    meas_phiSC4plot = merged_phiSC;
    meas_dtSC4plot = merged_dtSC;
    meas_BAZ_floor4plot = merged_BAZ_floor;

else
    meas_phiSC4plot = meas_phiSC;
    meas_dtSC4plot = meas_dtSC;
    meas_BAZ_floor4plot = meas_BAZ_floor;
end



%==========================================================================
%% sort to model BAZ range
%==========================================================================

findvals = find(meas_BAZ_floor > modrange_low & ...
                meas_BAZ_floor < modrange_upp);

meas_BAZ_floor = meas_BAZ_floor(findvals);
meas_phiSC = meas_phiSC(findvals,:);
meas_dtSC = meas_dtSC(findvals,:);



%==========================================================================
%% calc deviations of measured SWS parameters from theoretical ones
%==========================================================================

% for each synthetic model, then compute RMSE

count_mods = 1;

res_dt = zeros(length(meas_phiSC),1);
res_phi = res_dt;

for ii = 1:length(model_out)
    curr_mod_phi = model_out(ii).phi_eff;
    curr_mod_dt = model_out(ii).dt_eff;
    curr_mod_type = model_out(ii).type;

    sizeSC = size(meas_phiSC);

    for kk = 1:sizeSC(1)
        % find theoretical value for BAZ of corresponding measured value
        find_theo_phi = curr_mod_phi(meas_BAZ_floor(kk));
        find_theo_dt = curr_mod_dt(meas_BAZ_floor(kk));

        % calc residuum
        res_phi(kk,1) = abs( find_theo_phi - meas_phiSC(kk,2) );

        % if there is a residuum > 90, calc 180-res_phi to consider
        % the flip from -90 deg to +90 deg
        if res_phi(kk,1) > 90
            res_phi(kk,1) = 180-res_phi(kk,1);
        end

        res_dt(kk,1) = abs( find_theo_dt - meas_dtSC(kk,2) );
    end % kk

    if strcmp(curr_mod_type,'two_layers')
        modsall(count_mods).phi_eff = curr_mod_phi;
        modsall(count_mods).dt_eff = curr_mod_dt;
        modsall(count_mods).mod_type = curr_mod_type;
        modsall(count_mods).phi = model_out(ii).mod_paras.phi_in;
        modsall(count_mods).dt = model_out(ii).mod_paras.dt_in;
        modsall(count_mods).modrange_low = modrange_low;
        modsall(count_mods).modrange_upp = modrange_upp;
        modsall(count_mods).RMSE_phi = sqrt(sum(res_phi.^2)/length(meas_phiSC));
        modsall(count_mods).RMSE_dt = sqrt(sum(res_dt.^2)/length(meas_phiSC));
        modsall(count_mods).staname = staname_split;
        modsall(count_mods).azi4plot = [];
        modsall(count_mods).fast4plot = [];
        modsall(count_mods).dt4plot = [];
    elseif strcmp(curr_mod_type,'single_layer')
        modsall(count_mods).phi_eff = curr_mod_phi;
        modsall(count_mods).dt_eff = curr_mod_dt;
        modsall(count_mods).mod_type = curr_mod_type;
        modsall(count_mods).phi = model_out(ii).mod_paras.phi_in;
        modsall(count_mods).dt = model_out(ii).mod_paras.dt_in;
        modsall(count_mods).modrange_low = modrange_low;
        modsall(count_mods).modrange_upp = modrange_upp;
        modsall(count_mods).RMSE_phi = sqrt(sum(res_phi.^2)/length(meas_phiSC));
        modsall(count_mods).RMSE_dt = sqrt(sum(res_dt.^2)/length(meas_phiSC));
        modsall(count_mods).staname = staname_split;
        modsall(count_mods).azi4plot = [];
        modsall(count_mods).fast4plot =[];
        modsall(count_mods).dt4plot = [];
    elseif strcmp(curr_mod_type,'dipping')
        modsall(count_mods).phi_eff = curr_mod_phi;
        modsall(count_mods).dt_eff = curr_mod_dt;
        modsall(count_mods).mod_type = curr_mod_type;
        modsall(count_mods).downdipdir = model_out(ii).mod_paras.downdipdir;
        modsall(count_mods).dip = model_out(ii).mod_paras.dip;
        modsall(count_mods).thick = model_out(ii).mod_paras.thick;
        modsall(count_mods).modrange_low = modrange_low;
        modsall(count_mods).modrange_upp = modrange_upp;
        modsall(count_mods).RMSE_phi = sqrt(sum(res_phi.^2)/length(meas_phiSC));
        modsall(count_mods).RMSE_dt = sqrt(sum(res_dt.^2)/length(meas_phiSC));
        modsall(count_mods).staname = staname_split;
        modsall(count_mods).azi4plot = model_out(ii).mod_paras.azi4plot;
        modsall(count_mods).fast4plot = model_out(ii).mod_paras.fast4plot;
        modsall(count_mods).dt4plot = model_out(ii).mod_paras.dt4plot;
    end

    if whichfit==2
    % use phi and dt for joint fitting
        modsall(ii).RMSE = modsall(ii).RMSE_phi/90 + modsall(ii).RMSE_dt/4;
        fit_str = 'phidt';
    elseif whichfit==1
    % only use phi for fitting
        modsall(ii).RMSE = modsall(ii).RMSE_phi/90;
        fit_str = 'phi';
    end

    count_mods = count_mods+1;
end % ii



%==========================================================================
%% sort and select models
%==========================================================================

%--------------------------------------------------------------------------
% all model types together

% sort models based on total RMSE
[~,index] = sort([modsall.RMSE]);

% entry 1 corresponds to minimum total RMSE
modsall_sort_all = modsall(index);

% keep only XX best models
modsall_sort = modsall_sort_all(1:keep_mods);

%--------------------------------------------------------------------------
% separate models based on model type
countH1 = 1; % 1 to <=keep_mods_sep
countH2 = 1;
countT1 = 1;

for ii = 1:1:length(modsall_sort_all)
    curr_mod_type = modsall_sort_all(ii).mod_type;

    if strcmp(curr_mod_type,'single_layer') & countH1<=keep_mods_sep
        modsH1(countH1) = modsall_sort_all(ii);
        countH1 = countH1 + 1;
    elseif strcmp(curr_mod_type,'two_layers') & countH2<=keep_mods_sep
        modsH2(countH2) = modsall_sort_all(ii);
        countH2 = countH2 + 1;
    elseif strcmp(curr_mod_type,'dipping') & countT1<=keep_mods_sep
        modsT1(countT1) = modsall_sort_all(ii);
        countT1 = countT1 + 1;
    end

    if countH1>keep_mods_sep && countH2>keep_mods_sep && ...
        countT1>keep_mods_sep
        break
    end

end

phi_H1 = [modsH1.phi];
dt_H1 = [modsH1.dt];
RMSE_H1 = [modsH1.RMSE];
merge_H1 = [phi_H1' dt_H1' RMSE_H1'];

phi_H2 = [modsH2.phi]; % lower layer and upper layer
dt_H2 = [modsH2.dt];
RMSE_H2 = [modsH2.RMSE];
merge_H2 = [phi_H2(1,:)' phi_H2(2,:)' dt_H2(1,:)' dt_H2(2,:)' RMSE_H2'];

downdipdir_T1 = [modsT1.downdipdir];
dip_T1 = [modsT1.dip];
thick_T1 = [modsT1.thick];
RMSE_T1 = [modsT1.RMSE];
merge_T1 = [downdipdir_T1' dip_T1' thick_T1' RMSE_T1'];



%==========================================================================
%% write 20 best-fit models for each model type to txt file
%==========================================================================

fname_start = [modsall_sort(1).staname '_anisomodel_' ...
               data_used_str '_' fit_str '_'];
fname_end = ['_BAZ' num2str(modrange_low) 'to' num2str(modrange_upp) ...
             '_domper' num2str(domper) 's.txt'];

fname_H1 = [fname_start 'H1' fname_end];
filename_H2 = [fname_start 'H2' fname_end];
filename_T1 = [fname_start 'T1' fname_end];

writematrix(merge_H1, fname_H1, 'Delimiter','tab')
writematrix(merge_H2, filename_H2, 'Delimiter','tab')
writematrix(merge_T1, filename_T1, 'Delimiter','tab')



%==========================================================================
%% call plotting functions
%==========================================================================


%--------------------------------------------------------------------------
% backazimuth variation and model type distribution

SWS_modeling_plot_results( ...
    BAZ, modsall_sort, plot_mod_max, ...
    meas_BAZ_floor_null, meas_phiSC_null, meas_dtSC_null, ...
    meas_BAZ_floor_simw_null, meas_phiSC_simw_null, meas_dtSC_simw_null, ...
    modrange_low, modrange_upp, colmod_bf_1, colmod_bf_2max ,lw_mod, ...
    colorsfill, colorsedge, marksize, ms_null, lw_symb, lw_symb_null, ...
    col_face_null, col_edge_null, ...
    fontsize, fs_RMSE, modrange_col, modrange_edcol, ...
    meas_BAZ_floor4plot, meas_phiSC4plot, meas_dtSC4plot, staname_split, fit_str, ...
    domper, keep_mods, keep_mods_sep, data_used_str, ...
    modsH1, phi_H1, dt_H1, RMSE_H1, ...
    modsH2, phi_H2, dt_H2, RMSE_H2, ...
    modsT1, downdipdir_T1, dip_T1, thick_T1, RMSE_T1, ...
    cmap_rms_sel, cmap_rms_ind, cmap_rms_str, ...
    cmap_phi, cmap_phi_ind, cmap_phi_str, cbar_phi_ind, cbar_phi_str, ...
    mymarkersize_symbols, mylinewidth_symbols ...
)

%--------------------------------------------------------------------------
% stereoplot of synthetic splitting parameter of the best fit models

plotnum = 1; % adjust manually

% all model types together
for ii = 1:1:plotnum
    SWS_modeling_plot_stereo_synthetic( ...
        modsall_sort, 'all', ii, ...
        modrange_low, modrange_upp, ...
        fit_str, domper, data_used_str, ...
        cmap_phi, cmap_phi_str ...
    )
end
% one layer with horizontal symmetry axis (H1)
for ii = 1:1:plotnum
    SWS_modeling_plot_stereo_synthetic( ...
        modsH1, 'sep', ii, ...
        modrange_low, modrange_upp, ...
        fit_str, domper, data_used_str, ...
        cmap_phi, cmap_phi_str ...
    )
end
% two layers with horizontal symmetry axes (H2)
for ii = 1:1:plotnum
    SWS_modeling_plot_stereo_synthetic( ...
        modsH2, 'sep', ii, ...
        modrange_low, modrange_upp, ...
        fit_str, domper, data_used_str, ...
        cmap_phi, cmap_phi_str ...
    )
end
% one layer with tilted symmetry axis (T1)
for ii = 1:1:plotnum
    SWS_modeling_plot_stereo_synthetic( ...
        modsT1, 'sep', ii, ...
        modrange_low, modrange_upp, ...
        fit_str, domper, data_used_str, ...
        cmap_phi, cmap_phi_str ...
    )
end


%==========================================================================
end % EOF

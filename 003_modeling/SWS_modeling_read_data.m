function [RES_split, RES_null, RES_stack, RES_simw_split, RES_simw_null] = ...
    SWS_modeling_read_data( ...
        dir_res_split, dir_res_null, ...
        dir_res_stack, dir_res_simw, ...
        varargin ...
    )

% ==========================================================================
%% This function
% ==========================================================================
% reads txt files with shear wave splitting (SWS) measurements (SWSM)
% related to a single seismological recording station
% - single-event analysis (SplitLab, SL)
% - multi-event analysis (StackSplit)
% outputs structs for
% - splits and nulls (SL)
% - stacks and simw (StackSplit)
% based on the selected qualities
% --------------------------------------------------------------------------
% is
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
%% main function
%==========================================================================


%==========================================================================
% scale the delay time, e. g. for GMT, with a factor (same for all results)
% >>> adjust manually for your needs <<<

scalingfac = 1;


%==========================================================================
% check SL version

%--------------------------------------------------------------------------
% search for original SL folder

curr_dir = pwd;

[folderSL, ~, ~] = fileparts(which('install_SplitLab.m'));

%--------------------------------------------------------------------------
% check if SL version by Rob Porritt is available

[folderSLRP, ~, ~] = fileparts(which('SL_swap_QT_components.m'));
if ~isempty(folderSL) && isempty(folderSLRP)
    cd(folderSL)
    SL_version = 1;
    disp(' ')
    disp('>>> SplitLab version 1.0.5 found! <<<')
    disp(' ')
elseif ~isempty(folderSL) && ~isempty(folderSLRP) && ...
        strcmp(folderSL,folderSLRP)
    cd(folderSL)
    SL_version = 2;
    disp(' ')
    disp('>>> SplitLab version 1.2.1 (by Rob Porritt) found! <<<')
    disp(' ')
else
    errordlg('>>> No SplitLab version found! <<<')
    return
end

cd(curr_dir)


%==========================================================================
% search for quality as input argument and make query

if isempty(varargin)
    select_qual = input(['Qualities (splits, nulls) you want to use ' ...
                         '(default all)? \n',...
                         '   [0] all \n',...'
                         '   [1] good \n',...'
                         '   [2] good & fair \n',...'
                         '   [3] fair & poor \n',...'
                         '   [4] fair \n',...'
                         '   [5] poor    | ']);
    if isempty(select_qual) % default
        select_qual = 0; % all
    end
else
   select_qual = varargin{1};
end


%==========================================================================
% check for SWS measurements available

%--------------------------------------------------------------------------
% check for files available

% single splits
if isempty(dir_res_split)
   disp('   >>> No file with SINGLE SPLIT results found!')
   RES_split = [];
else
   scaling_factor = scalingfac;
   RES_split = read_results(dir_res_split, select_qual, scaling_factor, ...
                            SL_version);
end

% single nulls
if isempty(dir_res_null)
   disp('   >>> No file with SINGLE NULL results found!')
   RES_null = [];
else
   scaling_factor = scalingfac;
   RES_null = read_results(dir_res_null, select_qual, scaling_factor, ...
                            SL_version);
end

% stack splits
if isempty(dir_res_stack)
   disp('   >>> No file with STACK SPLIT results found!')
   RES_stack = [];
else
   scaling_factor = scalingfac;
   RES_stack = read_results_STACK(dir_res_stack, scaling_factor);
end

% simw splits and nulls
if isempty(dir_res_simw)
   disp('   >>> No file with SIMW (split or null) results found!')
   RES_simw_split = [];
   RES_simw_null = [];
else
   scaling_factor = scalingfac;
   RES_simw = read_results_SIMW(dir_res_simw, scaling_factor);
   % only good or fair SIMW results
   idx_gf_NN = strcmp({RES_simw.quality_manual},'good') & ...
                    strcmp({RES_simw.NULL},'No') | ...
               strcmp({RES_simw.quality_manual},'fair') & ...
                    strcmp({RES_simw.NULL},'No');
   idx_gf_N = strcmp({RES_simw.quality_manual},'good') & ...
                    strcmp({RES_simw.NULL},'Yes') | ...
              strcmp({RES_simw.quality_manual},'fair') & ...
                    strcmp({RES_simw.NULL},'Yes');
   RES_simw_split = RES_simw(idx_gf_NN); % split
   RES_simw_null = RES_simw(idx_gf_N); % null
end

%--------------------------------------------------------------------------
% check for qualities available

stanamecheck={};

if isempty(RES_split)
   disp('   >>> No match for selected SINGLE SPLIT qualities!')
else
   disp('   >>> Found SINGLE SPLIT results!')
   for i_split=1:1:length(RES_split)
       stanamecheck{end+1} = RES_split(i_split).staname;
   end
end

if isempty(RES_null)
   disp('   >>> No match for selected SINGLE NULL qualities!')
else
   disp('   >>> Found SINGLE NULL results!')
   for i_null=1:1:length(RES_null)
       stanamecheck{end+1} = RES_null(i_null).staname;
   end
end

if isempty(RES_stack)
   disp('   >>> No STACK SPLIT results!')
else
   disp('   >>> Found STACK SPLIT results!')
   for i_stack=1:1:length(RES_stack)
       stanamecheck{end+1} = RES_stack(i_stack).staname;
   end
end

if isempty(RES_simw_split)
   disp('   >>> No SIMW SPLIT results!')
else
   disp('   >>> Found SIMW SPLIT results!')
   for i_simw_split=1:1:length(RES_simw_split)
       stanamecheck{end+1} = RES_simw_split(i_simw_split).staname;
   end
end

if isempty(RES_simw_null)
   disp('   >>> No SIMW NULL results!')
else
   disp('   >>> Found SIMW NULL results!')
   for i_simw_null=1:1:length(RES_simw_null)
       stanamecheck{end+1} = RES_simw_null(i_simw_null).staname;
   end
end

%--------------------------------------------------------------------------
% check for station uniform

stanamecheck = unique(stanamecheck);

if length(stanamecheck) > 1
    error('>>> Input files are from different stations! <<<')
end

%==========================================================================
end % EOmF



%==========================================================================
%% subfunction - single results
%==========================================================================

function RES_out = read_results(dir_res, select_qual, scaling_factor, ...
                                SL_version)

% allocate fields for speed
res_split.date_doy = 'NaN';
res_split.staname = 'NaN';
res_split.sta_lon = NaN; % still empty
res_split.sta_lat = NaN; % still empty
res_split.phase = 'NaN';
res_split.baz = NaN;
res_split.inc = NaN;
res_split.filter = [NaN NaN];
res_split.phiRC = NaN;
res_split.dtRC = NaN;
res_split.phiSC_err_min = NaN;
res_split.phiSC = NaN;
res_split.phiSC_err_max = NaN;
res_split.dtSC_err_min = NaN;
res_split.dtSC = NaN;
res_split.dtSC_err_max = NaN;
res_split.phiEV = NaN;
res_split.dtEV = NaN;
res_split.SNRSC = NaN;
res_split.quality_manual = 'NaN';
res_split.NULL = 'NaN';

%--------------------------------------------------------------------------
% SL original version
if SL_version==1

    curr_folder = pwd;
    cd(dir_res.folder)

    fid=fopen(dir_res.name);
    C = textscan(fid,'%s %s %s %f %f %s %f %f %s %f %f %f %s %f %s %f %f %s %f %s %f %f %f %f %s %s','headerlines',3);
    fclose(fid);

    for k=1:length(C{1})

        res_split(k).date_doy = C{1,1}{k,1};
        res_split(k).staname = C{1,2}{k,1};
        res_split(k).sta_lon = NaN;       % still empty
        res_split(k).sta_lat = NaN;       % still empty
        res_split(k).phase = C{1,3}{k,1};
        res_split(k).baz = C{1,4}(k);
        res_split(k).inc = C{1,5}(k);
        res_split(k).filter = [C{1,7}(k) C{1,8}(k)];
        res_split(k).phiRC = C{1,10}(k);
        res_split(k).dtRC = C{1,11}(k);
        res_split(k).phiSC_err_min = C{1,12}(k);
        res_split(k).phiSC = C{1,14}(k);
        res_split(k).phiSC_err_max = C{1,16}(k);
        res_split(k).dtSC_err_min = C{1,17}(k);
        res_split(k).dtSC = C{1,19}(k)/scaling_factor;
        res_split(k).dtSC_err_max = C{1,21}(k);
        res_split(k).phiEV = C{1,22}(k);
        res_split(k).dtEV = C{1,23}(k);
        res_split(k).SNRSC = C{1,24}(k);
        res_split(k).quality_manual = C{1,25}{k,1};
        res_split(k).NULL = C{1,26}{k,1};

        if strcmp(res_split(k).phase,'SKS')
            res_split(k).phase_col = [0    0.4470    0.7410];
        elseif strcmp(res_split(k).phase,'SKKS')
            res_split(k).phase_col = [0.8500    0.3250    0.0980];
        elseif strcmp(res_split(k).phase,'PKS')
            res_split(k).phase_col = [0.4660    0.6740    0.1880];
        else
            res_split(k).phase_col = [0.9290    0.6940    0.1250];
        end

    end

%--------------------------------------------------------------------------
% SL version by Rob Porritt
elseif SL_version==2

    curr_folder = pwd;
    cd(dir_res.folder)

    fid=fopen(dir_res.name);
    C = textscan(fid,'%s %s %s %f %f [%f %f] %f %f %f<%f <%f %f< %f <%f %f %f %s %s','headerlines',3);
    fclose(fid);

    for k = 1:length(C{1})

        res_split(k).date_doy = C{1,1}{k,1};
        res_split(k).staname = C{1,2}{k,1};
        res_split(k).sta_lon = NaN;       % still empty
        res_split(k).sta_lat = NaN;       % still empty
        res_split(k).phase = C{1,3}{k,1};
        res_split(k).baz = C{1,4}(k);
        res_split(k).inc = C{1,5}(k);
        res_split(k).filter = [C{1,6}(k) C{1,7}(k)];
        res_split(k).phiRC = C{1,8}(k);
        res_split(k).dtRC = C{1,9}(k);
        res_split(k).phiSC_err_min = C{1,10}(k);
        res_split(k).phiSC = C{1,11}(k);
        res_split(k).phiSC_err_max = C{1,12}(k);
        res_split(k).dtSC_err_min = C{1,13}(k);
        res_split(k).dtSC = C{1,14}(k)/scaling_factor;
        res_split(k).dtSC_err_max = C{1,15}(k);
        res_split(k).phiEV = C{1,16}(k);
        res_split(k).dtEV = C{1,17}(k);
        res_split(k).SNRSC = [];
        res_split(k).quality_manual = C{1,18}{k,1};
        res_split(k).NULL = C{1,19}{k,1};

        if strcmp(res_split(k).phase,'SKS')
            res_split(k).phase_col = [0    0.4470    0.7410];
        elseif strcmp(res_split(k).phase,'SKKS')
            res_split(k).phase_col = [0.8500    0.3250    0.0980];
        elseif strcmp(res_split(k).phase,'PKS')
            res_split(k).phase_col = [0.4660    0.6740    0.1880];
        else
            res_split(k).phase_col = [0.9290    0.6940    0.1250];
        end

    end

end % SL_version


%==========================================================================
% sort quality

sel_ev = [];

if select_qual==0 % all
   sel_ev = res_split;
elseif select_qual==1 % good
   find_ev = strcmp({res_split.quality_manual},'good');
   sel_ev = res_split(find_ev);
elseif select_qual==2 % good & fair
   find_ev = strcmp({res_split.quality_manual},'good') | ...
             strcmp({res_split.quality_manual},'fair');
   sel_ev = res_split(find_ev);
elseif select_qual==3 % fair & poor
   find_ev = strcmp({res_split.quality_manual},'fair') | ...
             strcmp({res_split.quality_manual},'poor');
   sel_ev = res_split(find_ev);
elseif select_qual==4 % fair
   find_ev = strcmp({res_split.quality_manual},'fair');
   sel_ev = res_split(find_ev);
elseif select_qual==5 % poor
   find_ev = strcmp({res_split.quality_manual},'poor');
   sel_ev = res_split(find_ev);
end

RES_out = sel_ev;

cd(curr_folder)

%==========================================================================
end % EOmF



%==========================================================================
%% subfunction - multi results - STACK
%==========================================================================

function RES_out = read_results_STACK(dir_res_stack, scaling_factor)

curr_folder = pwd;
cd(dir_res_stack.folder)

fid=fopen(dir_res_stack.name);
C=textscan(fid,'%s %d %d %f %f %f %f %f %f %f %s %f %s %f %f %s %f %s %f %s %s %s %*[^\n]','headerlines',3);
fclose(fid);

for k = 1:1:length(C{1})

        res_split(k).staname = C{1,1}{k,1};
        res_split(k).sta_lon = NaN;
        res_split(k).sta_lat = NaN;
        res_split(k).stack_meth = C{1,20}(k);
        res_split(k).nsurf = C{1,2}(k);
        res_split(k).ndf = C{1,3}(k);
        res_split(k).minbaz = C{1,4}(k);
        res_split(k).maxbaz = C{1,5}(k);
        res_split(k).meanbaz = C{1,6}(k);
        res_split(k).mindis = C{1,7}(k);
        res_split(k).maxdis = C{1,8}(k);
        res_split(k).meandis = C{1,9}(k);
        res_split(k).phiSTACK_err_min = C{1,10}(k);
        res_split(k).phiSTACK = C{1,12}(k);
        res_split(k).phiSTACK_err_max = C{1,14}(k);
        res_split(k).dtSTACK_err_min = C{1,15}(k)/scaling_factor;
        res_split(k).dtSTACK = C{1,17}(k)/scaling_factor;
        res_split(k).dtSTACK_err_max = C{1,19}(k)/scaling_factor;

end

RES_out = res_split;

cd(curr_folder)

%==========================================================================
end % EOsF



%==========================================================================
%% subfunction - multi results - SIMW
%==========================================================================

function RES_out = read_results_SIMW(dir_res_simw, scaling_factor)

curr_folder = pwd;
cd(dir_res_simw.folder)

fid=fopen(dir_res_simw.name);
C=textscan(fid,'%s %d %f %f %f %f %f %f %d %f %f %f %s %f %s %f %f %s %f %s %f %f %f %s %s %*[^\n]','headerlines',3);
fclose(fid);

for k=1:length(C{1})

    res_split(k).staname = C{1,1}{k,1};
    res_split(k).nwave = C{1,2}(k);
    res_split(k).minbaz = C{1,3}(k);
    res_split(k).maxbaz = C{1,4}(k);
    res_split(k).meanbaz = C{1,5}(k);
    res_split(k).mindis = C{1,6}(k);
    res_split(k).maxdis = C{1,7}(k);
    res_split(k).meandis = C{1,8}(k);
    res_split(k).taper = C{1,9}(k);
    res_split(k).phiRC = C{1,10}(k);
    res_split(k).dtRC = C{1,11}(k);
    res_split(k).phiSC_err_min = C{1,12}(k);
    res_split(k).phiSC = C{1,14}(k);
    res_split(k).phiSC_err_max = C{1,16}(k);
    res_split(k).dtSC_err_min = C{1,17}(k)/scaling_factor;
    res_split(k).dtSC = C{1,19}(k)/scaling_factor;
    res_split(k).dtSC_err_max = C{1,21}(k)/scaling_factor;
    res_split(k).phiEV = C{1,22}(k);
    res_split(k).dtEV = C{1,23}(k);
    res_split(k).quality_manual = C{1,24}{k,1};
    res_split(k).NULL = C{1,25}{k,1};

end

RES_out = res_split;

cd(curr_folder)


%==========================================================================
end % EOsF

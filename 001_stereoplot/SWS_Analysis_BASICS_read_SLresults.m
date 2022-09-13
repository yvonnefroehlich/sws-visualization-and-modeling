function [RES_split, RES_nulls, SL_quality] = ...
    SWS_Analysis_BASICS_read_SLresults(varargin)

%==========================================================================
%% This function
%==========================================================================
% reads shear wave splitting single-event analysis results of SplitLab (SL)
% (txt files for nulls and splits)
% outputs structs for splits and nulls based on the selected qualities
%--------------------------------------------------------------------------
% is
% - created and mainly written: Michael Grund (ORCID 0000-0001-8759-2018)
%   https://github.com/michaelgrund/sws_tools
%   Grund PhD (2019)
%   https://doi.org/10.5445/IR/1000091425
%   Grund & Ritter (2020) Geophysical Journal International
%   https://doi.org/10.1093/gji/ggaa388
% - modified: Yvonne Fröhlich (ORCID 0000-0002-8566-0619)
%   https://github.com/yvonnefroehlich/sws-visualization-and-modeling
%   Ritter, Fröhlich, Sanz Alonso & Grund (2022) Journal of Seismology
%   https://doi.org/10.1007/s10950-022-10112-w
%--------------------------------------------------------------------------
% TERMS OF USE
%
% The plotting routines are provided "as is" and without any warranty.
% The author cannot be held responsible for anything that happens to you
% or your equipment. Use it at your own risk.
%==========================================================================



%==========================================================================
%% main function
%==========================================================================


%==========================================================================
% check SplitLab version

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
elseif ~isempty(folderSL) && ~isempty(folderSLRP) && ...
        strcmp(folderSL,folderSLRP)
    cd(folderSL)
    SL_version = 2;
    disp(' ')
    disp('>>> SplitLab version 1.2.1 (by Rob Porritt) found! <<<')
else
    errordlg('>>> No SplitLab version found! <<<')
    return
end

cd(curr_dir)


%==========================================================================
% search for quality and make query

disp(' ')
if isempty(varargin)
   SL_quality = input(['Qualities you want to plot (default is all)? \n', ...
                    '   [0] all \n', ...'
                    '   [1] good \n', ...'
                    '   [2] good & fair \n', ...'
                    '   [3] fair & poor \n', ...'
                    '   [4] fair \n', ...'
                    '   [5] poor    | ']);
    if isempty(SL_quality) % default
        SL_quality = 0; % all
    end
else
   SL_quality = varargin{1};
end


%==========================================================================
% search for input files

if ~isempty(varargin)
    if length(varargin) > 1
        dir_res_split = varargin{2};
        dir_res_nulls = varargin{3};
    else
        dir_res_split = dir('splitresults_*.txt');
        dir_res_nulls = dir('splitresultsNULL_*.txt');
    end
else
   dir_res_split = dir('splitresults_*.txt');
   dir_res_nulls = dir('splitresultsNULL_*.txt');
end

%--------------------------------------------------------------------------
% read input files considering choosen quality
if isempty(dir_res_split)
   disp('>>> No file with SPLIT results found! <<<')
   RES_split = [];
else
   % length of delay time vector, default 1, in general a good choice
   scaling_factor = 1;
   RES_split = read_results(dir_res_split, SL_quality, scaling_factor, ...
                            SL_version);
end

if isempty(dir_res_nulls)
   disp('>>> No file with NULL results found! <<<')
   RES_nulls = [];
else
   scaling_factor = 1; % length of delay time vector, default 1
   RES_nulls = read_results(dir_res_nulls, SL_quality, scaling_factor, ...
                            SL_version);
end

%--------------------------------------------------------------------------
% display station and found SWSMs with choosen quality
if ~isempty(RES_split)
    disp(' ')
    disp(['>>> This is station ' RES_split(1).staname ' <<<'])
elseif ~isempty(RES_nulls)
    disp(' ')
    disp(['>>> This is station ' RES_nulls(1).staname ' <<<'])
else
    RES_split = [];
    RES_nulls = [];
    return
end

if isempty(RES_split)
   disp(' ')
   disp('   >>> No match for your selected SPLIT qualities! <<<')
   pause(2)
end

if isempty(RES_nulls)
   disp(' ')
   disp('   >>> No match for your selected NULL qualities! <<<')
   pause(2)
end


end % EOmF



%==========================================================================
%% subfunction
%==========================================================================

function RES_out = read_results(dir_res, SL_quality, scaling_factor, ...
                                SL_version)

% allocate fields for speed
res_split.date_doy = 'NaN';
res_split.staname = 'NaN';
res_split.sta_lon = NaN;
res_split.sta_lat = NaN;
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
res_split.remark = 'NaN';

%--------------------------------------------------------------------------
% SL original version
if SL_version==1

    fid = fopen(dir_res.name);
    C = textscan(fid,'%s %s %s %f %f %s %f %f %s %f %f %f %s %f %s %f %f %s %f %s %f %f %f %f %s %s %s','headerlines',3);
    fclose(fid);

    for k = 1:1:length(C{1})
        res_split(k).date_doy = C{1,1}{k,1};
        res_split(k).staname = C{1,2}{k,1};
        res_split(k).sta_lon = NaN; % still empty
        res_split(k).sta_lat = NaN; % still empty
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
        try
            res_split(k).remark = C{1,20}{k,1};
        catch
            text = 'no remark'; %disp('no remark')
        end
    end

%--------------------------------------------------------------------------
% SL version by Rob Porritt
elseif SL_version==2

    fid = fopen(dir_res.name);
    C = textscan(fid,'%s %s %s %f %f [%f %f] %f %f %f<%f <%f %f< %f <%f %f %f %s %s %s','headerlines',3);
    fclose(fid);

    for k = 1:1:length(C{1})
        res_split(k).date_doy = C{1,1}{k,1};
        res_split(k).staname = C{1,2}{k,1};
        res_split(k).sta_lon = NaN; % still empty
        res_split(k).sta_lat = NaN; % still empty
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
        try
            res_split(k).remark = C{1,20}{k,1};
        catch
            text = 'no remark'; %disp('no remark')
        end
    end

end % SL_version


res_split_depth = res_split;


%==========================================================================
% further selection

%--------------------------------------------------------------------------
% phase

% sel_ev_phase = [];
% find_mes = strcmp({res_split.phase},'SKKS'); % adjust for your needs
% sel_ev_phase = res_split_depth(find_mes);
% res_split_depth = sel_ev_phase;

%--------------------------------------------------------------------------
% observation typ
% 'No' equals split or non-null, 'Yes' equals null

% sel_ev_phase = [];
% find_mes = strcmp({res_split.NULL},'No'); % adjust for your needs
% sel_ev_null = res_split_depth(find_mes);
% res_split_depth = sel_ev_null;


%==========================================================================
% sort quality

sel_ev_quality = [];

if SL_quality==0 % all
   sel_ev_quality = res_split_depth;
elseif SL_quality==1 % good
   find_ev = strcmp({res_split_depth.quality_manual},'good');
   sel_ev_quality = res_split_depth(find_ev);
elseif SL_quality==2 % good & fair
   find_ev = strcmp({res_split_depth.quality_manual},'good') | ...
        strcmp({res_split_depth.quality_manual},'fair');
   sel_ev_quality = res_split_depth(find_ev);
elseif SL_quality==3 % fair & poor
   find_ev = strcmp({res_split_depth.quality_manual},'fair') | ...
        strcmp({res_split_depth.quality_manual},'poor');
   sel_ev_quality = res_split_depth(find_ev);
elseif SL_quality==4 % fair
   find_ev = strcmp({res_split_depth.quality_manual},'fair');
   sel_ev_quality = res_split_depth(find_ev);
elseif SL_quality==5 % poor
   find_ev = strcmp({res_split_depth.quality_manual},'poor');
   sel_ev_quality = res_split_depth(find_ev);
end

RES_out = sel_ev_quality; % -> RES_split AND RES_null


%==========================================================================
end % EOsF
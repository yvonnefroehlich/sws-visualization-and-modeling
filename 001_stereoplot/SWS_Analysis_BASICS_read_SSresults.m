function RES_out = SWS_Analysis_BASICS_read_SSresults( ...
    dir_res_multi, scaling_factor, SSmethod ...
)

%==========================================================================
%% This function
%==========================================================================
% reads shear wave splitting multi-event analysis results of StackSplit
% (MATLAB structs)
% - stacking of error surfaces (STACK)
% - simultaneous inversion of multiple waveforms (SIMW) (Roy et al. 2017)
% outputs MATLAB structs
%--------------------------------------------------------------------------
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
%--------------------------------------------------------------------------
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
%--------------------------------------------------------------------------
% TERMS OF USE
%
% The plotting routines are provided "as is" and without any warranty.
% The author cannot be held responsible for anything that happens to you
% or your equipment. Use it at your own risk.
%--------------------------------------------------------------------------
% CONTRIBUTING
%
% Feel free to modify/adjust the code for your needs. Submit improvements
% and report bugs by opening a "New issue" in the GitHub repository (:
%==========================================================================



predat = load(dir_res_multi.name);
pre_RES_multi = predat.eqstack;
res_multi = [];

k = 1;

for kk = 1:1:length(pre_RES_multi)

%--------------------------------------------------------------------------
    % stacking of error surfaces
    if (~strcmp(pre_RES_multi(kk).results.stack_meth,'SIMW')) && ...
       (SSmethod==1 || SSmethod==3)

        res_multi(k).staname = pre_RES_multi(kk).results.stnname;
        res_multi(k).sta_lon = pre_RES_multi(kk).results.slat;
        res_multi(k).sta_lat = pre_RES_multi(kk).results.slong;
        res_multi(k).stack_meth = pre_RES_multi(kk).results.stack_meth;
        res_multi(k).nsurf = pre_RES_multi(kk).results.nsurf;
        res_multi(k).minbaz = pre_RES_multi(kk).results.bazi_min;
        res_multi(k).maxbaz = pre_RES_multi(kk).results.bazi_max;
        res_multi(k).meanbaz = pre_RES_multi(kk).results.bazi_mean;
        res_multi(k).mindis = pre_RES_multi(kk).results.dist_min;
        res_multi(k).maxdis = pre_RES_multi(kk).results.dist_max;
        res_multi(k).meandis = pre_RES_multi(kk).results.dist_mean;
        res_multi(k).phimulti_err_min = pre_RES_multi(kk).results.phi_stack(1);
        res_multi(k).phimulti = pre_RES_multi(kk).results.phi_stack(2);
        res_multi(k).phimulti_err_max = pre_RES_multi(kk).results.phi_stack(3);
        res_multi(k).dtmulti_err_min = pre_RES_multi(kk).results.dt_stack(1)/scaling_factor;
        res_multi(k).dtmulti = pre_RES_multi(kk).results.dt_stack(2)/scaling_factor;
        res_multi(k).dtmulti_err_max = pre_RES_multi(kk).results.dt_stack(3)/scaling_factor;
        res_multi(k).remark = pre_RES_multi(kk).results.remark;
        res_multi(k).used_phases = pre_RES_multi(kk).results.events_in;

        k = k+1;

%--------------------------------------------------------------------------
    % simultaneous inversion of multiple waveforms; ONLY splits
    elseif ( strcmp(pre_RES_multi(kk).results.stack_meth,'SIMW') && ...
            strcmp(pre_RES_multi(kk).results.Null_simw,'No ') ) && ...
           (SSmethod==2 || SSmethod==3)

        res_multi(k).staname = pre_RES_multi(kk).results.stnname;
        res_multi(k).sta_lon = pre_RES_multi(kk).results.slong;
        res_multi(k).sta_lat = pre_RES_multi(kk).results.slat;
        res_multi(k).stack_meth = pre_RES_multi(kk).results.stack_meth;
        res_multi(k).nsurf = pre_RES_multi(kk).results.nwave;
        res_multi(k).minbaz = pre_RES_multi(kk).results.bazi_min;
        res_multi(k).maxbaz = pre_RES_multi(kk).results.bazi_max;
        res_multi(k).meanbaz = pre_RES_multi(kk).results.bazi_mean;
        res_multi(k).mindis = pre_RES_multi(kk).results.dist_min;
        res_multi(k).maxdis = pre_RES_multi(kk).results.dist_max;
        res_multi(k).meandis = pre_RES_multi(kk).results.dist_mean;
        res_multi(k).phimulti_err_min = pre_RES_multi(kk).results.phiSC_simw(1);
        res_multi(k).phimulti = pre_RES_multi(kk).results.phiSC_simw(2);
        res_multi(k).phimulti_err_max = pre_RES_multi(kk).results.phiSC_simw(3);
        res_multi(k).dtmulti_err_min = pre_RES_multi(kk).results.dtSC_simw(1)/scaling_factor;
        res_multi(k).dtmulti = pre_RES_multi(kk).results.dtSC_simw(2)/scaling_factor;
        res_multi(k).dtmulti_err_max = pre_RES_multi(kk).results.dtSC_simw(3)/scaling_factor;
        res_multi(k).remark = pre_RES_multi(kk).results.remark;
        res_multi(k).used_phases = pre_RES_multi(kk).results.events_in;

        k = k+1;

    end

end % kk

RES_out = res_multi;


%===============================================================================
end % EOF

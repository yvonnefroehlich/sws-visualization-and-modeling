function modout = SWS_modeling_precomp_dippinglayer( ...
    dfreq, stepdddir, stepdips, stepthick ...
)

% ==========================================================================
%% This function
% ==========================================================================
% generates synthetic apparent (effective) splitting parameters for
%   dipping one-layer models
% outputs corresponding MATLAB struct
% >>> for small step sizes < stepdddir >, < stepdips >, < stepthick >
% computation time and struct size increases significantly <<<
% --------------------------------------------------------------------------
% uses the provided MATLAB function
% - SWS_modeling_calc_dipping
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



%==========================================================================
%% model set up
%==========================================================================

%--------------------------------------------------------------------------
% define incidence angle that should be used
inc = 10; % deg

%--------------------------------------------------------------------------
% set up parameters
downdipdir = 0:stepdddir:(360-stepdddir);
dips = stepdips:stepdips:75; % avoid horizontal layer
thickness = stepthick:stepthick:250; % avoid layer of zero thickness

%--------------------------------------------------------------------------
% get all possible combinations
comb_vecs = combvec(downdipdir, dips, thickness);
N = length(comb_vecs);

disp(' ')
disp('Total number of dipping-layer models to generate: ' + num2str(N))
disp('Generate models...')



%==========================================================================
%% generate models
%==========================================================================

modout = repmat(struct('phi_eff',zeros(1,360), ...
                       'dt_eff',zeros(1,360), ...
                       'mod_paras',struct('downdipdir',0, ...
                                          'dip',0, ...
                                          'thick',0, ...
                                          'azi4plot',zeros(1,360)), ...
                       'type',zeros(1,1)), ...
                N, 1);


% if problems occur, replace parfor by standard for loop
% parfor ii = 1:N
for ii = 1:N

    currmod = comb_vecs(:,ii);
    downdipdir = currmod(1,:);
    dips = currmod(2,:);
    thickness = currmod(3,:);

    [fast_eff, tlag_eff, azi4plot, fast4plot, tlag4plot] = ...
        SWS_modeling_calc_dipping(inc, dips, downdipdir, thickness, dfreq);

     modout(ii).phi_eff = fast_eff;
     modout(ii).dt_eff = tlag_eff;
     modout(ii).mod_paras.downdipdir = downdipdir;
     modout(ii).mod_paras.dip = dips;
     modout(ii).mod_paras.thick = thickness;
     modout(ii).mod_paras.azi4plot = azi4plot;
     modout(ii).mod_paras.fast4plot = fast4plot;
     modout(ii).mod_paras.dt4plot = tlag4plot;
     modout(ii).type = 'dipping';

    if rem(ii/1000, 1) == 0 % whole number
        disp([num2str(ii) ' / ' num2str(N) ' models done.'])
    end

end

% save(['sws_modout_domper_dipping_' ...
%    num2str(1/dfreq) 's_' ...
%    num2str(stepdddir) 'deg_' num2str(stepdips) 'deg_' num2str(stepthick) 'km.mat'], ...
%    'modout', '-v7.3')

disp('Dipping-layer models done!')


%==========================================================================
end % EOF

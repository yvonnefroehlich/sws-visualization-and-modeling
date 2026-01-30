function modout = SWS_modeling_precomp_singlelayer(stepphi, stepdt)

% ==========================================================================
%% This function
% ==========================================================================
% generates synthetic apparent (or effective) splitting parameters for
%   horizontal single-layer models
% outputs corresponding MATLAB struct
% --------------------------------------------------------------------------
% is
% - modified: Yvonne Fröhlich (ORCID 0000-0002-8566-0619)
%   https://github.com/yvonnefroehlich/sws-visualization-and-modeling
%   Fröhlich (2025) Dissertation
%   https://doi.org/10.5445/IR/1000183786
%   Fröhlich, Grund, Ritter (2024) Geophysical Journal International
%   https://doi.org/10.1093/gji/ggae245
%   Ritter, Fröhlich, Sanz Alonso & Grund (2022) Journal of Seismology
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
% Copyright (C) 2026  Yvonne Fröhlich & Michael Grund (v2.0)
% Copyright (C) 2022  Yvonne Fröhlich & Michael Grund (v1.0)
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
%% model set up
%==========================================================================

%--------------------------------------------------------------------------
% one layer
phi = (-90+stepphi):stepphi:90; % in deg; -90 equivalent to 90
dt = stepdt:stepdt:4; % in sec; avoid dt=0

%--------------------------------------------------------------------------
% get all possible combinations
comb_vecs = combvec(phi, dt);
N = length(comb_vecs);

disp(' ')
disp(['Total number of single-layer models to generate: ' num2str(N)])
disp('Generate models...')



%==========================================================================
%% pre-calculate models
%==========================================================================

% pre-allocate structs
modout = repmat(struct('phi_eff',zeros(1,360), ...
                       'dt_eff',zeros(1,360), ...
                       'mod_paras', struct('phi_in',zeros(1,1), ...
                                           'dt_in',zeros(1,1), ...
                                           'counter',zeros(1,1)), ...
                       'type',zeros(1,1)), ...
                N, 1);


for ii = 1:1:N

    currmod = comb_vecs(:,ii);
    modphis = currmod(1:length(currmod)/2,:);
    moddts = currmod(length(currmod)/2+1:end,:);

    modout(ii).phi_eff = ones(1,length(0:1:360)) * modphis;
    modout(ii).dt_eff = ones(1,length(0:1:360)) * moddts;
    modout(ii).mod_paras.phi_in = modphis;
    modout(ii).mod_paras.dt_in = moddts;
    modout(ii).mod_paras.counter = 1;
    modout(ii).type = 'single_layer';

    if rem(ii/1000,1)==0 % whole number
     disp([num2str(ii) ' models done.'])
    end

end

% save(['sws_modout_domper_single_layer_' ...
%    num2str(stepphi) 'deg_' num2str(stepdt) 's.mat'], ...
%    'modout', '-v7.3')

disp('Single-layer models done!')


%==========================================================================
end % EOF

function modout = SWS_modeling_precomp_twolayers(dfreq, stepphi, stepdt)

%==========================================================================
%% This function
%==========================================================================
% generates synthetic apparent (or effective) splitting parameters for
%   horizontal two-layer models
% outputs corresponding sturctur
% >>> for small step sizes < stepphi >, < stepdt > computation time and
% structure size increases significantly <<<
%--------------------------------------------------------------------------
% >>> requires the MSAT package by Walker & Wookey (2012) <<<
% please download it from (last access 2022 June 21)
% - https://www1.gly.bris.ac.uk/MSAT/
% - https://github.com/andreww/MSAT
%--------------------------------------------------------------------------
% is
% - created and mainly written: Michael Grund (ORCID 0000-0001-8759-2018)
%   https://github.com/michaelgrund/sws_tools
%   Grund PhD (2019)
%   https://doi.org/10.5445/IR/1000091425
%   Grund & Ritter (2020) Geophysical Journal International
%   https://doi.org/10.1093/gji/ggaa388
% - lightly modified: Yvonne Fröhlich (ORCID 0000-0002-8566-0619)
%   https://github.com/yvonnefroehlich/sws-visualization-and-modeling
%   Ritter, Fröhlich, Sanz Alonso & Grund (2022) Journal of Seismology
%   https://doi.org/10.1007/s10950-022-10112-w
%--------------------------------------------------------------------------
% LICENSE
%
% Copyright (C) 2020  Michael Grund
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
%--------------------------------------------------------------------------
% TERMS OF USE
%
% The modeling routines are provided "as is" and without any warranty.
% The author cannot be held responsible for anything that happens to you
% or your equipment. Use it at your own risk.
%==========================================================================



%==========================================================================
%% model set up
%==========================================================================

%--------------------------------------------------------------------------
% BAZ range to calculate theoretical SWS parameters
BAZ = 0:1:360;

%--------------------------------------------------------------------------
% two layers
phi1 = (-90+stepphi):stepphi:90; % in deg; -90 equivalent to 90
dt1 = stepdt:stepdt:4; % in sec; avoid dt=0 sec
phi2 = (-90+stepphi):stepphi:90; % deg; -90 equivalent to 90
dt2 = stepdt:stepdt:4; % in sec; avoid dt=0 sec

%--------------------------------------------------------------------------
% get all possible combinations
comb_vecs = combvec(phi1,phi2,dt1,dt2);

disp(' ')
disp(['Total number of two-layer models to generate: ' ...
      num2str(length(comb_vecs))])
disp('Generate models...')


%==========================================================================
%% generate models
%==========================================================================

N = length(comb_vecs);

% preallocate structs
modout = repmat(struct('phi_eff',zeros(1,360), ...
                       'dt_eff',zeros(1,360), ...
                       'mod_paras', struct('phi_in',zeros(1,2), ...
                                             'dt_in',zeros(1,2), ...
                                            'counter',zeros(1,1)), ...
                       'type', zeros(1,1)), ...
                N, 1);

parfor ii = 1:N

    currmod_in = comb_vecs(:,ii);
    modphis = currmod_in(1:length(currmod_in)/2,:);
    moddts = currmod_in(length(currmod_in)/2+1:end,:);
    [phi_eff_out, dt_eff_out] = nest(dfreq, BAZ, modphis, moddts);

    modout(ii).phi_eff = phi_eff_out;
    modout(ii).dt_eff = dt_eff_out;
    modout(ii).mod_paras.phi_in = modphis;
    modout(ii).mod_paras.dt_in = moddts;
    modout(ii).mod_paras.counter = 1;
    modout(ii).type = 'two_layers';

    if rem(ii/5000,1)==0 % whole number
        disp([num2str(ii) ' models done.'])
    end

end

disp('Two-layer models done!')


%==========================================================================
end % EOmF



%==========================================================================
%% subfunction
%==========================================================================

% outsource second parfor loop

function [phi_eff_out, dt_eff_out] = nest(dfreq, BAZ, modphis, moddts)

    parfor (jj = 1:length(BAZ), 8) % parallel computing to speed up calc.
    %for jj = 1:length(BAZ) % uncomment, if problems occur with parfor
        [phi_eff_out(jj), dt_eff_out(jj)] = ...
         MS_effective_splitting_N( ...
            dfreq, BAZ(jj), modphis', moddts', 'mode','S&S' );
    end

%==========================================================================
end % EOsF
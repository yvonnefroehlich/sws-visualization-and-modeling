function SWS_modeling_precomp_models_main()

% ==========================================================================
%% This function
% ==========================================================================
% pre-computes synthetic apparent (or effective) splitting parameters
% for structural anisotropy models
% - horizontal single-layer models (H1, 'single_layer')
% - horizontal two-layer models (H2, 'two_layers')
% - dipping one-layer models (T1, 'dipping')
% --------------------------------------------------------------------------
% uses the provided MATLAB functions
% - SWS_modeling_precomp_singlelayer.m
% - SWS_modeling_precomp_twolayers.m
% - SWS_modeling_precomp_dippinglayer.m
% - SWS_modeling_calc_dipping.m
% --------------------------------------------------------------------------
% >>> requires the MSAT package by Walker & Wookey (2012) <<<
% please download it from (last access 2022 June 21)
% - https://www1.gly.bris.ac.uk/MSAT/
% - https://github.com/andreww/MSAT
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
% Copyright (C) 2020  Michael Grund (sws_tools)
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
%% set stepsizes
%==========================================================================
% >>> adjust for your needs <<<
% >>> for small stepsizes computation time and struct size increases
% significantly <<<

% dominant period of shear wave
domper = 8; % in sec

% single-layer models
stepphis = 45; % in deg
stepdts = 1; % in sec

% two-layer models
stepphim = 45; % in deg
stepdtm = 1; % in sec

% dipping layer models
stepdddir = 45; % in deg
stepdips = 15; % in deg
stepthick = 100; % in km



%==========================================================================
%% pre-compute models
%==========================================================================

disp(' ')
disp(['Model setup for shear-wave splitting modeling using \n' ...
      'single-layer, two-layer and dipping layer models!'])

modout1 = SWS_modeling_precomp_singlelayer(stepphis, stepdts);
modout2 = SWS_modeling_precomp_twolayers(1/domper, stepphim, stepdtm);
modout3 = SWS_modeling_precomp_dippinglayer(1/domper, stepdddir, ...
                                            stepdips, stepthick);

% merge models
disp(' ')
disp('Merge models and save into file...')
splitmods = vertcat(modout1, modout2, modout3);

% save to mat-file
save(['sws_modout_domper' num2str(domper) 's.mat'], 'splitmods', '-v7.3')

disp(' ')
disp('Model setup done!')


%==========================================================================
end % EOF

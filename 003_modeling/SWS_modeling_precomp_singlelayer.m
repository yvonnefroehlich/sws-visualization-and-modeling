function modout = SWS_modeling_precomp_singlelayer(stepphi, stepdt)

%==========================================================================
%% This function
%==========================================================================
% generates synthetic apparent (or effective) splitting parameters for
%   horizontal single-layer models
% outputs corresponding structur
%--------------------------------------------------------------------------
% is
% - created and mainly written: Michael Grund (ORCID 0000-0001-8759-2018)
%   https://github.com/michaelgrund/sws_tools
%   Grund PhD (2019)
%   https://doi.org/10.5445/IR/1000091425
%   Grund & Ritter (2020) GJI
%   https://doi.org/10.1093/gji/ggaa388
% - lightly modified: Yvonne Fröhlich (ORCID 0000-0002-8566-0619)
%   https://github.com/yvonnefroehlich
%   Ritter, Fröhlich, Sanz Alonso & Grund (2022) Journal of Seismology
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

disp('Single-layer models done!')


%==========================================================================
end % EOF
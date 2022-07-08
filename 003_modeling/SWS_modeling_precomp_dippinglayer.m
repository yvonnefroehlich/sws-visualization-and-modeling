function modout = SWS_modeling_precomp_dippinglayer( ...
    dfreq, stepdddir, stepdips, stepthick ...
    )

%==========================================================================
%% This function
%==========================================================================
% generates synthetic apparent (effective) splitting parameters for
%   dipping one-layer models
% outputs corresponding MATLAB structure
% >>> for small step sizes < stepdddir >, < stepdips >, < stepthick >
% computation time and structure size increases significantly <<<
%--------------------------------------------------------------------------
% uses the provided MATLAB function
% - SWS_modeling_calc_dipping
%--------------------------------------------------------------------------
% is
% - created and mainly written: Michael Grund (ORCID 0000-0001-8759-2018)
%   https://github.com/michaelgrund/sws_tools
%   Grund PhD (2019)
%   https://doi.org/10.5445/IR/1000091425
%   Grund & Ritter (2020) GJI
%   https://doi.org/10.1093/gji/ggaa388
% - lightly modified: Yvonne Fröhlich (ORCID 0000-0002-8566-0619)
%   https://github.com/yvonnefroehlich/sws-visualization-and-modeling
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

% define incidence angle that should be used
inc = 10; % deg

% set up parameters
downdipdir = 0:stepdddir:(360-stepdddir);
dips = stepdips:stepdips:75; % avoid horizontal layer
thickness = stepthick:stepthick:250; % avoid layer of zero thickness

% get all possible combinations
comb_vecs = combvec(downdipdir, dips, thickness);

disp(' ')
disp(['Total number of dipping-layer models to generate: ' ...
      num2str(length(comb_vecs))])
disp('Generate models...')



%==========================================================================
%% generate models
%==========================================================================

N = length(comb_vecs);

modout = repmat(struct('phi_eff',zeros(1,360), ...
                       'dt_eff',zeros(1,360), ...
                       'mod_paras',struct('downdipdir',0, ...
                                          'dip',0, ...
                                          'thick',0, ...
                                          'azi4plot',zeros(1,360)), ...
                       'type',zeros(1,1)), ...
                N, 1);


parfor ii = 1:N % if problems occur, replace parfor by standard for loop

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

    if rem(ii/1000,1)==0 % whole number
     disp([num2str(ii) ' models done.'])
    end

end

disp('Dipping-layer models done!')


%==========================================================================
end % EOF
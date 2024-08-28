function [fast_eff, tlag_eff, azi4plot, fast4plot, tlag4plot] = ...
        SWS_modeling_calc_dipping( incval, dips, downdipdir, thick, dfreq )

%==========================================================================
%% This function
%==========================================================================
% calculates synthetic apparent (effective) splitting parameters for
%   dipping one-layer models
%--------------------------------------------------------------------------
% >>> requires the MSAT package by Walker & Wookey (2012) <<<
% please download it from (last access 2022 June 21)
% - https://www1.gly.bris.ac.uk/MSAT/
% - https://github.com/andreww/MSAT
%--------------------------------------------------------------------------
% is
% - mainly based on the MSAT function >>> split_model.m <<< and
%   modified to fit into the workflow
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
%   Fröhlich, Grund & Ritter (2024) Geophysical Journal International
%   https://doi.org/10.1093/gji/ggae245
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
% handle optional arguments and defaults
e_split_mode = 's&s'; % Default mode for effective splitting calc.
min_azi = 0.0;
max_azi = 360.0;
del_azi = 1.0;

%--------------------------------------------------------------------------
% model parameters
aoi = incval;

%--------------------------------------------------------------------------
% layer parameters
%L1_depth = 0;
L1_thick = thick ;
L1_dip = dips; % layer geometry
L1_aaz = 0; % a axis azimuth (relative to down dip direction)
L1_faln = 0.3; % fraction aligned

%--------------------------------------------------------------------------
% imaging parameters

% 0 is down dip direction (perp. to strike)
% note that this is the seismic backazimuth + 180 deg
% wave is assumed to be polarized in this direction
azi = [min_azi:del_azi:max_azi] -downdipdir-180;

% since here azi=0 deg is downdip direction (fast axis points in direction
% of downdip), for an estimated ~90 deg downdip direction one has to add
% 90 deg (subtract 270 deg) to the azi values to get the corresponding
% curves over BAZ
inc = ones(size(azi)).*90 - aoi; % 90 deg is vertical (set aoi==0)

%--------------------------------------------------------------------------
% calculate distances
[dist1] = distance_in_dipping_layer(L1_dip, aoi, L1_thick, azi);



%==========================================================================
%% anisotropy
%==========================================================================

%--------------------------------------------------------------------------
% load anisotropy, and generate an isotropic version of it

[Cani, rh] = MS_elasticDB('olivine');
[Ciso] = MS_decomp( MS_axes(Cani) );
Cani = MS_rot3(Cani, 90, 0, 0); % orientation for dry upper mantle

%--------------------------------------------------------------------------
% generate layer elasticities

% Voigt-Reuss-Hill average of the appropriately rotated olivine tensor and
% its isotropic equivalent
[L1_C, ~] = MS_VRH( ...
    [L1_faln 1-L1_faln], ...
    MS_rot3(Cani,0,-L1_dip,L1_aaz,'order',[3 2 1]), rh, Ciso, rh) ;

%--------------------------------------------------------------------------
% interrogate elasticities to generate splitting parameters for layer
[pol, ~, vs1, vs2, ~, ~, ~] = MS_phasevels(L1_C, rh, inc, azi);
fast1 = MS_unwind_pm_90( (azi+pol') ); % geog. reference frame
tlag1 = dist1./vs2' - dist1./vs1';



%==========================================================================
%% splitting parameters
%==========================================================================

%--------------------------------------------------------------------------
% calculate the effective splitting for the (dipping) layer
fast_eff = zeros(size(azi));
tlag_eff = fast_eff;

for i = 1:length(azi)
    [fast_eff(i),tlag_eff(i)] = MS_effective_splitting_N( ...
        dfreq,azi(i), fast1(i), tlag1(i), 'mode',e_split_mode ...
        );
end

%--------------------------------------------------------------------------
% variables needed for stereoplot displaying synthetic parameters
fast4plot = fast_eff;
tlag4plot = tlag_eff;
azi4plot = azi-180;

%--------------------------------------------------------------------------
% recalculate azi to bazi for correct plotting over bazi 0 to 360 deg

azi = min_azi:del_azi:max_azi;

% recalculate fast_eff
for ii = 1:1:length(azi)
    if (fast_eff(ii)+downdipdir) > 90
        fast_eff(ii) = fast_eff(ii)+downdipdir-180;
    elseif (fast_eff(ii)+downdipdir) < -90
        fast_eff(ii) = fast_eff(ii)+downdipdir+180;
    else
        fast_eff(ii) = fast_eff(ii)+downdipdir;
    end
end

% now correctly sort between +-90 deg
for ii = 1:1:length(azi)
    if (fast_eff(ii)) > 90
        fast_eff(ii) = fast_eff(ii)-180;
    elseif (fast_eff(ii)) < -90
        fast_eff(ii) = fast_eff(ii)+180;
    else
        fast_eff(ii) = fast_eff(ii);
    end
end

return


%==========================================================================
end % EOmF



%==========================================================================
%% subfunction
%==========================================================================

function [dist] = distance_in_dipping_layer(dip, aoi, thick, azi)

% calculate apparent dip
alp = atand( tand(dip) .* sind(azi-90) );

% calculate distances
gam = 90 - alp - aoi;
bet = 90 + alp;
dist = thick .* sind(bet) ./ sind(gam);

return

%==========================================================================
end % EOsF
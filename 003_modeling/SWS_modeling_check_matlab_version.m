function vers_out = SWS_modeling_check_matlab_version()

%==========================================================================
%% This function
%==========================================================================
% checks the MATLAB version available on your system
%--------------------------------------------------------------------------
% is
% - created: Michael Grund (ORCID 0000-0001-8759-2018)
%   https://github.com/michaelgrund/sws_tools
%   Grund PhD (2019)
%   https://doi.org/10.5445/IR/1000091425
%   Grund & Ritter (2020) Geophysical Journal International
%   https://doi.org/10.1093/gji/ggaa388
% - modified: Yvonne Fröhlich (ORCID 0000-0002-8566-0619)
%   https://github.com/yvonnefroehlich/sws-visualization-and-modeling
%   Ritter, Fröhlich, Sanz Alonso & Grund (2022) Journal of Seismology
%==========================================================================


vers = version('-release');

vers_yyyy = str2double(vers(1:4));

% MATLAB R2020a or higher
if vers_yyyy>2019
    vers_out=1;
else
    vers_out=0;
end


%==========================================================================
end % EOF
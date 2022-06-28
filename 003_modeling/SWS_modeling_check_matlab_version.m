function vers_out = SWS_modeling_check_matlab_version()

%==========================================================================
%% This function
%==========================================================================
% checks the MATLAB version available on your system
%--------------------------------------------------------------------------
% is
% - created: Michael Grund
%   Grund PhD (2019), Grund & Ritter (2020)
%   https://github.com/michaelgrund/sws_tools
% - extended: Yvonne Fröhlich
%   ORCID 0000-0002-8566-0619
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
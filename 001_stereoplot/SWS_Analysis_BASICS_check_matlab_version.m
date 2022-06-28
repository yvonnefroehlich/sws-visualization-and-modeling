function vers_out = SWS_Analysis_BASICS_check_matlab_version()

%==========================================================================
%% This function
%==========================================================================
% checks the MATLAB version available on your system
%--------------------------------------------------------------------------
% is
% - created: Michael Grund
%   https://github.com/michaelgrund/sws_tools
%   Grund PhD (2019)
%   https://doi.org/10.5445/IR/1000091425
%   Grund & Ritter (2020) GJI
%   https://doi.org/10.1093/gji/ggaa388
% - extended and modified: Yvonne Fröhlich
%   https://github.com/yvonnefroehlich
%   Ritter, Fröhlich, Sanz Alonso & Grund (2022) Journal of Seismology
%--------------------------------------------------------------------------
% TERMS OF USE
%
% The plotting routines are provided "as is" and without any warranty.
% The author cannot be held responsible for anything that happens to you
% or your equipment. Use it at your own risk.
%==========================================================================



vers = version('-release');

vers_yyyy = str2double(vers(1:4));
vers_letter = vers(5);

% MATLAB R2016b to R2019b
if (vers_yyyy>2016 && vers_yyyy<2020) || ...
   (vers_yyyy==2016 && strcmp(vers_letter,'b'))
   vers_out = 1;
% MATLAB R2020a or higher
elseif vers_yyyy>2019
   vers_out = 2;
else
   vers_out = 0;
end


%==========================================================================
end % EOF
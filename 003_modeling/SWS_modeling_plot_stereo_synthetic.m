function SWS_modeling_plot_stereo_synthetic( ...
    modsall_sort, sort_str, plotnum, ...
    modrange_low, modrange_upp, ...
    fit_str, domper, data_used, ...
    cmap_phi_in, cmap_phi_str ...
    )

%==========================================================================
%% This function
%==========================================================================
% generates stereoplots displaying the synthetic splitting pattern of the
% best-fitting structural anisotropy model
% based on the minimum root mean square error (RSME)
%--------------------------------------------------------------------------
% is
% - based on: >>> stereoplot.m <<< function of SplitLab
%   Wüstefeld et al. (2008)
%   https://doi.org/10.1016/j.cageo.2007.08.002
% - created and mainly written: Michael Grund (ORCID 0000-0001-8759-2018)
%   https://github.com/michaelgrund/sws_tools
%   Grund PhD (2019)
%   https://doi.org/10.5445/IR/1000091425
%   Grund & Ritter (2020) Geophysical Journal International
%   https://doi.org/10.1093/gji/ggaa388
% - modified and extended: Yvonne Fröhlich (ORCID 0000-0002-8566-0619)
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
%--------------------------------------------------------------------------
% CONTRIBUTING
%
% Feel free to modify/adjust the code for your needs. Submit improvements
% and report bugs by opening a "New issue" in the GitHub repository (:
%==========================================================================



fig_stereo = figure('visible','off');

%==========================================================================
%% plot settings
%==========================================================================
% >>> adjust for your needs <<<

%--------------------------------------------------------------------------
% stereoplot
linew = 3;
marks = 7;
linewcirc = 2;

%--------------------------------------------------------------------------
% legend and colorbar
myfontsize = 8;
myfontsize_BAZ = 15;
col_leg = 'k'; % black
lengthbar_leg = 0.0710; % manually adjusted to fit bar length!

%--------------------------------------------------------------------------
% sectors
colfill = [219,219,219] ./256;
col_val_white = 255.9999999999999;
col_wedge = [col_val_white, col_val_white, col_val_white]./256;

%--------------------------------------------------------------------------
% station label
color_sta = [255 90 0]./256; % orange

%--------------------------------------------------------------------------
% schematic sketches
lengthbar_lay = 0.145;

% Ritter et al. (2022) Journal of Seismology
col_lay_upp = [255 165 079] ./256; % orange
col_lay_low = [0.8 0.8 0.8]; % light gray
col_lay_one = [255 165 079] ./256; % orange
col_dip_sta = 'w'; % white
col_dip_cryst = 'b'; % black
col_dip_ray = 'r'; % red

% horizontal layers
y_low = -0.2933;
y_upp = -0.3220;
y_H = 0.008;

% manually adjusted to have nearly the same legend size
% for horizontal layer(s) and dipping layer
x_add = -0.01;



%==========================================================================
%% colormap for color-coding of bars based on fast axis (phi)
%==========================================================================
% >>> using colormap <phasemap> requires cmocean colormaps by Thyng
% et al. (2016) be available on your system <<<

%--------------------------------------------------------------------------
if strcmp(cmap_phi_str,'phino')

    %......................................................................
    % make query use phi color-coding and which colormap
    disp(' ')
    cmap_phi_ind = input(['Color-coding based on phi in stereoplot?: \n' ...
                          '    [1] parula(fliped)  [2] phase    | ']);

    %......................................................................
    % select phi colormap and set string for file name
    if cmap_phi_ind==1
        % build-in MATLAB
        cmap = flipud( parula(181) );
        cmap_phi_str = 'phiparulaflip';
    elseif cmap_phi_ind==2
        % phase of cmocean colormaps by Thyng et al. 2016
        cmap = cmocean('phase', 181);
        cmap_phi_str = 'phiphasemap';
    end

%--------------------------------------------------------------------------
else
    cmap = cmap_phi_in;
end



%==========================================================================
%% select synthetic splitting parameters and model parameter
%==========================================================================
stepsize = 5;

if strcmp(modsall_sort(plotnum).mod_type,'dipping')

    tlag_eff = modsall_sort(plotnum).dt4plot;
    fast_eff = modsall_sort(plotnum).fast4plot;
    azi = modsall_sort(plotnum).azi4plot;
    downdipdir = modsall_sort(plotnum).downdipdir;
    dips = modsall_sort(plotnum).dip;

    phi0 = fast_eff(1:stepsize:360);
    dt0 = tlag_eff(1:stepsize:360);

    phi0 = phi0';
    dt0 = dt0';
    bazi = azi(1:stepsize:180);

    for ii = 1:length(bazi)
        if bazi(ii) > 360
            bazi(ii) = bazi(ii)-360;
        end
    end

    modtyp_str = 'T1';

elseif strcmp(modsall_sort(plotnum).mod_type,'two_layers')

    phi0 = modsall_sort(plotnum).phi_eff(1:stepsize:360);
    dt0 = modsall_sort(plotnum).dt_eff(1:stepsize:360);

    phi0 = phi0';
    dt0 = dt0';
    bazi = 0:stepsize:179;

    phi0 = [phi0; phi0];
    dt0  = [dt0; dt0];
    bazi = [bazi, bazi+180]';

    bazi = bazi';

    modtyp_str = 'H2';

elseif strcmp(modsall_sort(plotnum).mod_type,'single_layer')

    phi0 = modsall_sort(plotnum).phi_eff(1:stepsize:360);
    dt0 = modsall_sort(plotnum).dt_eff(1:stepsize:360);

    phi0 = phi0';
    dt0 = dt0';
    bazi = 0:stepsize:179;

    phi0 = [phi0; phi0];
    dt0  = [dt0; dt0];
    bazi = [bazi, bazi+180]';

    bazi = bazi';

    modtyp_str = 'H1';

end

len = dt0/1; % scale bar lengths if they e.g. range over frame


%==========================================================================
bazi = [bazi, bazi+180]';
inc = ones(size(bazi)) * 10; % default 10 deg inclination
azim = phi0;
azim_pre = phi0';

m = max(inc);
m = round(m/10)*10;
lim = [-inf m+7];
lim_gray = 0.2610/15*lim(2); % for radius of area

axesm('stereo', 'Frame','on', 'Grid','on', 'Origin',[90 0], ...
      'MlineLocation',30, 'PlineLocation',5, 'fLatLimit',lim, ...
      'fLineWidth',1, 'GLinestyle','-', 'GLinewidth',0.4, ...
      'Gcolor',[0.8 0.8 0.8]);

framem('FLinewidth',2)



%==========================================================================
%% mark null area in shaded gray
%==========================================================================

if modsall_sort(1).modrange_low < modsall_sort(1).modrange_upp

    if (modsall_sort(1).modrange_low~=0 || modsall_sort(1).modrange_upp~=360)

        % first plot whole range in gray as bottom layer
        startwedge = 0;
        endwedge = 360;
        plot_wedge3D(deg2rad(startwedge-90), deg2rad(endwedge-90), ...
                     0, 0, lim_gray, colfill);

        % then plot modelled range again on top in white ;)
        startwedge = modsall_sort(1).modrange_low;
        endwedge = modsall_sort(1).modrange_upp;
        plot_wedge3D(deg2rad(startwedge-90), deg2rad(endwedge-90), ...
                     0, 0, lim_gray, col_wedge);
    end

% this is the case when the higher value is slightly higher than 0 and
% the other one in the SA region
elseif modsall_sort(1).modrange_low > modsall_sort(1).modrange_upp

    if modsall_sort(1).modrange_low~=0 || modsall_sort(1).modrange_upp~=360

        % first plot whole range in gray as botom layer
        startwedge = 0;
        endwedge = 360;
        plot_wedge3D(deg2rad(startwedge-90), deg2rad(endwedge-90), ...
                     0, 0, lim_gray, colfill);

        % then plot modelled range again on top in white in two steps
        startwedge = 0;
        endwedge = modsall_sort(1).modrange_upp;
        plot_wedge3D(deg2rad(startwedge-90), deg2rad(endwedge-90), ...
                     0, 0, lim_gray, col_wedge);

        startwedge = modsall_sort(1).modrange_low;
        endwedge = 360;
        plot_wedge3D(deg2rad(startwedge-90), deg2rad(endwedge-90), ...
                     0, 0, lim_gray, col_wedge);
    end

end



%==========================================================================
%% plot colored bars
%==========================================================================
colormap(cmap);

NNull = 1:length(phi0);
bazi = bazi(:);

inc  = inc(:);
len  = len(:);
azim = azim(:);

bazi = [bazi(NNull)  bazi(NNull)]';
inc  = [inc(NNull)   inc(NNull)]';
len  = [-len(NNull)  len(NNull)]';
azim = (bazi - [azim(NNull) azim(NNull)]');

len = len*2; % one second == 4 deg (2 deg in both directions)

if strcmp(modsall_sort(plotnum).mod_type,'dipping')
    bazi = bazi + downdipdir;
end

% plot the splits as bars
[latout, lonout] = reckon(90-inc, bazi, len, azim, 'degrees');
hndl = plotm(latout, lonout, -0.1);


%==========================================================================
step_phi = -90:1:90;

azim_rounded_all = nan(length(hndl),1);

for ii = 1:length(hndl)

    if strcmp(modsall_sort(plotnum).mod_type,'dipping')
        azim_rounded = floor(azim_pre(1,ii)+downdipdir);
    else
        azim_rounded = floor(azim_pre(1,ii));
    end

    if azim_rounded < -90 && azim_rounded >=-270
        azim_rounded = azim_rounded+180;
        index = step_phi==azim_rounded;
    elseif azim_rounded > 90 && azim_rounded <=270
        azim_rounded = azim_rounded-180;
        index = step_phi==azim_rounded;
    elseif azim_rounded < -270
        azim_rounded = azim_rounded+360;
        index = step_phi==azim_rounded;
    elseif azim_rounded > 270
        azim_rounded = azim_rounded-360;
        index = step_phi==azim_rounded;
    else
        index = step_phi==azim_rounded;
    end

    azim_rounded_all(ii) = azim_rounded;
    set(hndl(ii),'color',cmap(index,:))
    set(hndl(ii),'linewidth',linew)

end



%==========================================================================
%% plot potential nulls
%==========================================================================

%--------------------------------------------------------------------------
if strcmp(modsall_sort(plotnum).mod_type,'dipping')
% in direction of down-dip direction and perpendicular to it
% >>> only valid for assuming the fast axis plunges with down-dip dir <<<

    % first plot nulls in downdipdir +-180 deg
    bazi_nulls = [downdipdir downdipdir+180];
    inc_nulls = [10 10];
    for k = 1:2
        plotm(90-inc_nulls(k), bazi_nulls(k), -0.15, ...
              'o', 'color','k', 'MarkerFaceColor','w', ...
              'MarkerSize',marks, 'linewidth',linewcirc);
    end

    % second find positions where bar orientation exact +-90 deg BAZ
    % >>> for downdipdir 90, 270 deg the nulls are not correctly shown <<<
    if downdipdir < 90
        first = azim_rounded_all + 90;
        firstmin = min(first);
        second = azim_rounded_all - 90;
        secondmax = max(second);
    elseif downdipdir > 90 && downdipdir < 180
        first = azim_rounded_all + 90;
        firstmin = max(first);
        second = azim_rounded_all - 90;
        secondmax = min(second);
    elseif downdipdir > 270
        first = azim_rounded_all - 90;
        firstmin = max(first);
        second = azim_rounded_all + 90;
        secondmax = min(second);
    elseif downdipdir > 180 && downdipdir < 270
        first = azim_rounded_all(azim_rounded_all<0) -90; % original, unclear
        %first = azim_rounded_all(azim_rounded_all>0) + 90; changed
        firstmin = max(first);
        second = azim_rounded_all(azim_rounded_all>0) - 90;
        secondmax = min(second);
    else
        first = azim_rounded_all + 90;
        firstmin = max(first);
        second = azim_rounded_all - 90;
        secondmax = min(second);
    end


    if abs(secondmax-firstmin) < 5
        firstmin = firstmin+180;
    end

    if abs(secondmax-firstmin) < 5
        secondmax = secondmax-180;
    end


    bazi_nulls = [firstmin secondmax];
    inc_nulls = [10 10];

    for k = 1:2
        plotm(90-inc_nulls(k), bazi_nulls(k), -0.15, ...
              'o', 'color','k', 'MarkerFaceColor','w', ...
              'MarkerSize',marks, 'linewidth',linewcirc);
    end

%--------------------------------------------------------------------------
elseif strcmp(modsall_sort(plotnum).mod_type,'two_layers')
% >>> under development <<<
% >>> for 90 deg and 270 deg the nulls may not correctly shown <<<

    inc_nulls = 10;

    % use characteristic 90 deg jump of phi in case of two layers
    dd_phi = phi0(2:end) - phi0(1:end-1); % calculate difference
    %disp(dd_phi)
    dd_phi_max = max(dd_phi); % calculate maximum
    %disp(dd_phi_max)
    dd_phi_max_abs = max( abs(dd_phi) ); % calculate absolute maximum
    %disp(dd_phi_max_abs)

    % >>> phi0 as -90 deg to 90 deg (not 0 deg to 180 deg) <<<
    count_phi_neg = 0;
    count_phi_pos = 0;
    for pp = 1:1:length(phi0)
        if phi0(pp)<0
            count_phi_neg = count_phi_neg+1;
        elseif phi0(pp)>0
            count_phi_pos = count_phi_pos+1;
        end
    end

    if count_phi_pos >= count_phi_neg
        if dd_phi_max_abs == dd_phi_max % abs max taken from pos. value
            dd_phi_max_sign = dd_phi_max_abs;
            %disp('positive')
        elseif dd_phi_max_abs ~= dd_phi_max % abs max taken from neg. value
            dd_phi_max_sign = -1*dd_phi_max_abs;
            %disp('negative')
        end
    elseif count_phi_pos < count_phi_neg
        if dd_phi_max_abs == dd_phi_max
            dd_phi_max_sign = dd_phi_max_abs;
            %disp('positive')
        elseif dd_phi_max_abs ~= dd_phi_max
            dd_phi_max_sign = dd_phi_max;
            %disp('negative')
        end
    end

    % index for first occurrence
    dd_phi_max_index = find(dd_phi==dd_phi_max_sign, 1);
    %disp(dd_phi_max_index)
    bazi_phi_max = azim_rounded_all(dd_phi_max_index); % corresponding BAZ
    %disp(bazi_phi_max)

    % use occurrence in 90 deg intervals
    plotm(90-inc_nulls, bazi_phi_max, -0.15, ...
          'o', 'color','k', 'MarkerFaceColor','w', ...
          'MarkerSize',marks, 'linewidth',linewcirc);
    plotm(90-inc_nulls, bazi_phi_max+90, -0.15, ...
          'o', 'color','k', 'MarkerFaceColor','w', ...
          'MarkerSize',marks, 'linewidth',linewcirc);
    plotm(90-inc_nulls, bazi_phi_max+180, -0.15, ...
          'o', 'color','k', 'MarkerFaceColor','w', ...
          'MarkerSize',marks, 'linewidth',linewcirc);
    plotm(90-inc_nulls, bazi_phi_max+270, -0.15, ...
          'o', 'color','k', 'MarkerFaceColor','w', ...
          'MarkerSize',marks, 'linewidth',linewcirc);

%--------------------------------------------------------------------------
elseif strcmp(modsall_sort(plotnum).mod_type,'single_layer')

    inc_nulls = 10;

    plotm(90-inc_nulls, modsall_sort(plotnum).phi, -0.15, ...
          'o', 'color','k', 'MarkerFaceColor','w', ...
          'MarkerSize',marks, 'linewidth',linewcirc);
    plotm(90-inc_nulls, modsall_sort(plotnum).phi+90, -0.15, ...
          'o', 'color','k', 'MarkerFaceColor','w', ...
          'MarkerSize',marks, 'linewidth',linewcirc);
    plotm(90-inc_nulls, modsall_sort(plotnum).phi+180, -0.15, ...
          'o', 'color','k', 'MarkerFaceColor','w', ...
          'MarkerSize',marks, 'linewidth',linewcirc);
    plotm(90-inc_nulls, modsall_sort(plotnum).phi+270, -0.15, ...
          'o', 'color','k', 'MarkerFaceColor','w', ...
          'MarkerSize',marks, 'linewidth',linewcirc);

end



%==========================================================================
%% plot legend (lower right corner)
%==========================================================================

startval = 0.25;

plot([startval startval], ...
     [0.230 0.230], ...
     'o', 'linewidth',linewcirc, 'markersize',marks, 'color',col_leg)
text(startval, 0.205, 'null', ...
     'HorizontalAlignment','center', 'fontsize',myfontsize)

plot([startval-lengthbar_leg/2 startval+lengthbar_leg/2], ...
     [0.265 0.265], ...
     '-', 'linewidth',linew, 'color',col_leg)
text(startval, 0.252, '1 s', ...
     'HorizontalAlignment','center', 'fontsize',myfontsize)

lengthbar_leg = 2*lengthbar_leg;
plot([startval-lengthbar_leg/2 startval+lengthbar_leg/2], ...
     [0.295 0.295], ...
     '-', 'linewidth',linew, 'color',col_leg)
text(startval, 0.282, '2 s', ...
     'HorizontalAlignment','center', 'fontsize',myfontsize)



%==========================================================================
%% plot colorbar (upper right corner)
%==========================================================================

cb = colorbar('location','north', ...
              'TickDirection','out', ...
              'TickLength',0.025);

zlab = get(cb, 'xlabel');
set(zlab,'String','   \phi_a / N\circE');

caxis([-90 90])
set(cb, 'xtick',-60:30:60);
set(cb, 'fontsize',myfontsize)

% [left bottom width height]
if strcmp(modsall_sort(plotnum).mod_type,'dipping')
    set(cb, 'position',[0.615, 0.915 0.220 0.020])
else
    set(cb, 'position',[0.625, 0.920 0.220 0.020])
end



%==========================================================================
%% display station name (lower left corner)
%==========================================================================

text(-0.270, 0.252, modsall_sort(1).staname, ...
     'HorizontalAlignment','center', ...
     'fontsize',myfontsize+8, 'color',color_sta)



%==========================================================================
%% best-fit anisotropy model (upper left corner)
%==========================================================================
% plot schematic sketch of model type and display model parameters

%--------------------------------------------------------------------------
if strcmp(modsall_sort(plotnum).mod_type,'two_layers')

    % plot two horizontal layers as two tick lines
    % lower
    plot([-startval-lengthbar_lay/2+x_add -startval+lengthbar_lay*0.8], ...
         [y_low+y_H y_low+y_H], ...
         '-', 'linewidth',linew+9, 'color',col_lay_low)
    % upper
    plot([-startval-lengthbar_lay/2+x_add -startval+lengthbar_lay*0.8], ...
         [y_upp+y_H y_upp+y_H], ...
         '-', 'linewidth',linew+9, 'color',col_lay_upp)

    % display model parameters of best-fit anisotropy model
    % index 1 - lower layer
    % index 2 - upper layer
    phi1 = modsall_sort(plotnum).phi(1);
    phi2 = modsall_sort(plotnum).phi(2);
    dt1 = modsall_sort(plotnum).dt(1);
    dt2 = modsall_sort(plotnum).dt(2);

    text(-startval-lengthbar_lay/2+x_add/2, y_low+y_H, ...
         ['\phi_1=N' num2str(phi1) '\circE,  ' ...
          '\delta{\itt}_1=' num2str(dt1,'%.2f') 's'], ...
         'HorizontalAlignment','left', 'fontsize',myfontsize, 'color','k')
    text(-startval-lengthbar_lay/2+x_add/2, y_upp+y_H, ...
         ['\phi_2=N' num2str(phi2) '\circE,  ' ...
          '\delta{\itt}_2=' num2str(dt2,'%.2f') 's'], ...
         'HorizontalAlignment','left', 'fontsize',myfontsize, 'color','k')

%--------------------------------------------------------------------------
elseif strcmp(modsall_sort(plotnum).mod_type,'single_layer')

    % plot one horizontal layer as orange thick line
    plot([-startval-lengthbar_lay/2+x_add -startval+lengthbar_lay*0.8], ...
         [y_upp+y_H y_upp+y_H], ...
         '-', 'linewidth',linew+9, 'color',col_lay_one)

    % display model parameters of best-fit anisotropy model
    phi = modsall_sort(plotnum).phi(1);
    dt = modsall_sort(plotnum).dt(1);

    text(-startval-lengthbar_lay/2+x_add/2, y_upp+y_H, ...
         ['\phi=N' num2str(phi) '\circE,  ' ...
          '\delta{\itt}=' num2str(dt,'%.2f') 's'], ...
         'HorizontalAlignment','left', 'fontsize',myfontsize, 'color','k')

%--------------------------------------------------------------------------
elseif strcmp(modsall_sort(plotnum).mod_type,'dipping')

    %......................................................................
    % plot dip angle sector as gray wedge starting from horizontal
    xdirset = 0.0;
    ydirset = -0.3;
    startwedge = 90;
    endwedge = 90 + dips;
    colfill = [210 210 210]./256;

    plot_wedge(deg2rad(startwedge-90), deg2rad(endwedge-90), ...
               -startval+xdirset-lengthbar_lay/2, ydirset, 0.12, colfill);

    %......................................................................
    % surface as black line
    plot([-startval+xdirset-lengthbar_lay/2 -startval+lengthbar_lay/1.65], ...
         [ydirset ydirset], '-', 'linewidth',linew-1, 'color','k')

    %......................................................................
    % plot dipping layer as thick line with aligned crystals as dashed line
    ang = dips;
    leng = 0.13;
    scaleddist = 1;
    beg = [-startval+xdirset-lengthbar_lay/2 ydirset];
    fin = beg + leng*[cosd(ang) sind(ang)];

    if dips==0 % case horizontal layer
        plot([-startval+xdirset-lengthbar_lay/2 fin(1)], ...
             [ydirset fin(2)], ...
             '-', 'linewidth',linew+35, 'color',col_lay_one)
        plot([-startval+xdirset-lengthbar_lay/2 fin(1)], ...
             [ydirset fin(2)], ...
             '--', 'linewidth',linew-1, 'color',col_dip_cryst)
    else
        plot([-startval+xdirset-lengthbar_lay/2 fin(1)], ...
             [ydirset fin(2)], ...
             '-', 'linewidth',linew+20, 'color',col_lay_one)
        plot([-startval+xdirset-lengthbar_lay/2 fin(1)], ...
             [ydirset+0.0*scaleddist fin(2)+0.0*scaleddist], ...
             '--', 'linewidth',linew-1, 'color',col_dip_cryst)
    end

    %......................................................................
    % plot schematic ray path with ~10 deg incidence angle as lines
    xdirset = 0.03;
    leng = 0.140;

    % right line
    ang = 79;
    beg = [-startval+xdirset-lengthbar_lay/2 ydirset];
    fin = beg + leng*[cosd(ang) sind(ang)];
    plot([-startval+xdirset-lengthbar_lay/2 fin(1)], [ydirset fin(2)], ...
         '-', 'linewidth',linew-1, 'color',col_dip_ray)
    % left line
    ang = 101;
    beg = [-startval+xdirset-lengthbar_lay/2 ydirset];
    fin = beg + leng*[cosd(ang) sind(ang)];
    plot([-startval+xdirset-lengthbar_lay/2 fin(1)], [ydirset fin(2)], ...
         '-', 'linewidth',linew-1, 'color',col_dip_ray)

    %......................................................................
    % polishing
    xdirset = 0.0;

    % hide edge of orange dippping layer due to small white line
    plot([-startval+xdirset-lengthbar_lay/2 -startval+lengthbar_lay/1.65], ...
         [ydirset-0.019 ydirset-0.019], '-', 'linewidth',linew+11, ...
         'color','w')
    % plot surface again due to small black line
    plot([-startval+xdirset-lengthbar_lay/2 -startval+lengthbar_lay/1.65], ...
         [ydirset ydirset], '-', 'linewidth',linew-1, ...
         'color','k')

    %......................................................................
    % plot station as inverse triangle
    xdirset = -0.293;
    ydirset = -0.317;

    plot(xdirset, ydirset, 'v', 'markerfacecolor',col_dip_sta, ...
         'markersize',11, 'markeredgecolor','k', 'linewidth',1.5)

    %......................................................................
    % write value of dip angle
    text(-startval+0.075-lengthbar_lay/2, ydirset-0.001, ...
         ['\Psi=' num2str(dips) '\circ'], ...
          'HorizontalAlignment','left', ...
          'fontsize',myfontsize+4, 'color','k')

    %......................................................................
    % plot down-dip direction as arrow

    % adopted from the drawArrow function available here:
    % Matthew Kelly (2020). drawArrow.
    % https://www.mathworks.com/matlabcentral/fileexchange/55181-drawarrow,
    % MATLAB Central File Exchange. Retrieved June 22, 2020.

    % >>> partly the arrow is lightly shifted form the center <<<

    % input values for function
    p0 = [-0.07 0];
    p1 = [ 0.07,0];
    color_arrow = [120 120 120]./256;

    % parameters of arrow layout
    W1 = 0.14; % half width of the arrow head normalized by length of arrow
    W2 = 0.04; % half width of the arrow shaft
    L1 = 0.38; % length of the arrow head normalized by length of arrow
    L2 = 0.33; % length of the arrow inset

    % tail and tip of the arrow
    x0 = p0(1);
    y0 = p0(2);
    x1 = p1(1);
    y1 = p1(2);

    % drawing an arrow from 0 to 1 on the x-axis
    P = [ 0, (1-L2), (1-L1), 1, (1-L1), (1-L2),   0;
         W2,     W2,     W1, 0,    -W1,    -W2, -W2 ];

    % scale, rotate, shift, and plot
    dx = x1-x0;
    dy = y1-y0;
    Length = sqrt(dx*dx + dy*dy);
    Angle = atan2(-dy,dx);

    P = Length*P;
    P = [cos(Angle),sin(Angle); -sin(Angle),cos(Angle)] * P;
    P = p0(:)*ones(1,7) + P;

    arrw = patch( P(1,:),P(2,:),color_arrow );
    arrw.EdgeColor = 'k';
    set(arrw, 'LineWidth',1.5)
    rotate(arrw, [0 0 1], downdipdir-90) % corresponding to syn. aniso model

end



%==========================================================================
%% annotation N(north) and E(east)
%==========================================================================

axis tight
axis off

lval = min(abs(axis));
text(0, -lval-0.005, 'N', ...
    'HorizontalAlignment','Center', 'VerticalAlignment','Base', ...
    'fontsize',myfontsize_BAZ);
text(lval+0.005, 0, 'E', ...
    'HorizontalAlignment','Left', 'VerticalAlignment','middle', ...
    'fontsize',myfontsize_BAZ);

view([0 -90])

axes = gca;
axes.SortMethod = 'ChildOrder'; % for right order of layers in eps / pdf



%==========================================================================
%% save figure
%==========================================================================

%--------------------------------------------------------------------------
% check MATLAB version
vers_out = SWS_modeling_check_matlab_version();

%--------------------------------------------------------------------------
% create file name
filename = ['PLOT_SWSmod_' modsall_sort(1).staname '_' num2str(plotnum) ...
            'synStereo' modtyp_str sort_str '_' fit_str '_' data_used ...
            '_baz' num2str(modrange_low) 'to' num2str(modrange_upp) ...
            '_domper' num2str(domper) 's_' cmap_phi_str];

%--------------------------------------------------------------------------
% save
if vers_out==1 % R2020a or higher
    exportgraphics(fig_stereo, [filename '.png'], 'Resolution',360)
    exportgraphics(fig_stereo, [filename '.eps'], 'ContentType','vector')
    exportgraphics(fig_stereo, [filename '.pdf'], 'ContentType','vector')
else % below R2020a
    saveas(export_fig,[filename '.png'])
    % uncomment if you want plots in pdf format, takes time...
    %saveas(export_fig,[filename '.eps'])
    %saveas(export_fig,[filename '.pdf'])
end


%==========================================================================
end % EOmF



%==========================================================================
%% helper functions
%==========================================================================
% based on >>> plot_arc.m <<< by Matt Fig
% https://de.mathworks.com/matlabcentral/answers/6322-drawing-a-segment-of-a-circle
% (last access 2022 June 19)


%==========================================================================
function hndl = plot_wedge(start_val, end_val, h, k, r, colfill)

t = linspace(start_val,end_val);
x = r*cos(t) + h;
y = r*sin(t) + k;
x = [x h x(1)];
y = [y k y(1)];

hndl = fill(x, y, 'r');
set(hndl, 'facecolor',colfill, 'edgecolor','white', 'facealpha',1)

axis([h-r-1 h+r+1 k-r-1 k+r+1])
axis tight;

end % EOsF


%==========================================================================
function hndl = plot_wedge3D(start_val, end_val, h, k, r, colfill)

t = linspace(start_val,end_val);
x = r*cos(t) + h;
y = r*sin(t) + k;
x = [x h x(1)];
y = [y k y(1)];

hndl = fill3(x, y, ones(length(x))*0.5, 'r');
set(hndl, 'facecolor',colfill, 'edgecolor','white', 'facealpha',1)

axis([h-r-1 h+r+1 k-r-1 k+r+1])
axis tight;

end % EOsF
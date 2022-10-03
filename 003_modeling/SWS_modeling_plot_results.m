function SWS_modeling_plot_results( ...
    BAZ, models_sort, plot_mod_max, ...
    meas_BAZ_floor_null, meas_phiSC_null, meas_dtSC_null, ...
    meas_BAZ_floor_simw_null, meas_phiSC_simw_null, meas_dtSC_simw_null, ...
    modrange_low, modrange_upp, colmod_bf_1, colmod_bf_2max, lw_mod, ...
    colfill, coledge, marksize, marksize_null, lw_symbols, lw_symbols_null, ...
    color_face_null, col_edge_null, ...
    fontsize, fontsize_RMSE, modrange_col, modrange_edcol, ...
    meas_BAZ_floor, meas_phiSC, meas_dtSC, staname_split, fit_str, ...
    domper, keep_mods, keep_mods_sep, data_used, ...
    modsH1, phi_H1, dt_H1, RMSE_H1, ...
    modsH2, phi_H2, dt_H2, RMSE_H2, ...
    modsT1, downdipdir_T1, dip_T1, thick_1, RMSE_T1, ...
    cmap_rms_sel, cmap_rms_ind, cmap_rms_str, ...
    cmap_phi, cmap_phi_ind, cmap_phi_str, cbar_phi_ind, cbar_phi_str, ...
    mymarkersize_symbols, mylinewidth_symbols ...
    )

%==========================================================================
%% This function
%==========================================================================
% plots shear-wave splitting modeling results
%--------------------------------------------------------------------------
% is
% - created and mainly written: Michael Grund (ORCID 0000-0001-8759-2018)
%   https://github.com/michaelgrund/sws_tools
%   Grund PhD (2019)
%   https://doi.org/10.5445/IR/1000091425
%   Grund & Ritter (2020) Geophysical Journal International
%   https://doi.org/10.1093/gji/ggaa388
% - strongly modified and extended: Yvonne Fröhlich (ORCID 0000-0002-8566-0619)
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
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
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



% check your matlab version
vers_out = SWS_modeling_check_matlab_version();


%==========================================================================
%% figure settings
%==========================================================================
% >>> adjust for your needs <<<

myresolution = 360;
status_save = 'yes'; %% 'yes','no'
modtyp_str_all = {'All','H1','H2','T1'};

%--------------------------------------------------------------------------
% BAZ range to plot

% % Ritter et al. 2022 Journal of Seismology
% baz2plot_start = 0;
% baz2plot_end = 100;

baz2plot_start = 0;
baz2plot_end = 360;

%--------------------------------------------------------------------------
% legend symbols
leg_symb_col_split = 'w'; % color-coding based on phi
leg_symb_single_split = 'diamond';
leg_symb_single_null = 'o';
leg_symb_simw_split = 'square';
leg_symb_simw_null = 'square';

% % Ritter et al. 2022 Journal of Seismology
% % -> only stack splits in NE quadrant
% leg_symb_stack_split = 'diamond';

leg_symb_stack_split = 'square';

% TickLabel position of RMSE colorbar
tick_label_pos_str = 'mid'; %% 'low','mid'

% RMSE colorbar label
if strcmp(fit_str,'phidt')
    cb_label_rmse = 'RMSE_{tot}';
elseif strcmp(fit_str,'phi')
    cb_label_rmse = 'RMSE_{phi}';
end


%==========================================================================
% predefined colors for singles and multis
% >>> do not change here <<<
singlcol = [66 91 169]./256; % color for single
simwcol = [0.9290 0.6940 0.1250]; % color for simw


%==========================================================================
%% shift phi -90°-90° -> 0°-180°
%==========================================================================
% >>> adjust stations for your needs <<<

station_shift = {};

%--------------------------------------------------------------------------
% synthetic anisotropy models
if any( strcmp(station_shift,staname_split) )
    for ii = 1:1:plot_mod_max
        for jj = 1:1:length([models_sort(1).phi_eff])
            if models_sort(ii).phi_eff(jj)<0
                models_sort(ii).phi_eff(jj) = models_sort(ii).phi_eff(jj) + 180;
            end
        end
        for jj = 1:1:length([modsH1(1).phi_eff])
             if modsH1(ii).phi_eff(jj)<0
                modsH1(ii).phi_eff(jj) = modsH1(ii).phi_eff(jj) + 180;
            end
            if modsH2(ii).phi_eff(jj)<0
                modsH2(ii).phi_eff(jj) = modsH2(ii).phi_eff(jj) + 180;
            end
            if modsT1(ii).phi_eff(jj)<0
                modsT1(ii).phi_eff(jj) = modsT1(ii).phi_eff(jj) + 180;
            end
        end
    end
end

%--------------------------------------------------------------------------
% single splits and stack splits or simw splits
if any( strcmp(station_shift,staname_split) )
    for ii = 1:1:length(meas_phiSC(:,2))
        if meas_phiSC(ii,1)<0
            meas_phiSC(ii,1) = meas_phiSC(ii,1) + 180;
        end
        if meas_phiSC(ii,2)<0
            meas_phiSC(ii,2) = meas_phiSC(ii,2) + 180;
        end
        if meas_phiSC(ii,3)<0
            meas_phiSC(ii,3) = meas_phiSC(ii,3) + 180;
        end
    end
end

% single nulls
if ~isempty(meas_phiSC_null)
    if any( strcmp(station_shift,staname_split) )
        for ii = 1:1:length(meas_phiSC_null(:,2))
            if meas_phiSC_null(ii,1)<0
                meas_phiSC_null(ii,1) = meas_phiSC_null(ii,1) + 180;
            end
            if meas_phiSC_null(ii,2)<0
                meas_phiSC_null(ii,2) = meas_phiSC_null(ii,2) + 180;
            end
            if meas_phiSC_null(ii,3)<0
                meas_phiSC_null(ii,3) = meas_phiSC_null(ii,3) + 180;
            end
        end
    end
end

% simw nulls
if ~isempty(meas_phiSC_simw_null)
    if any( strcmp(station_shift,staname_split) )
        for ii = 1:1:length(meas_phiSC_simw_null(:,2))
            if meas_phiSC_simw_null(ii,1)<0
                meas_phiSC_simw_null(ii,1) = meas_phiSC_simw_null(ii,1) + 180;
            end
            if meas_phiSC_simw_null(ii,2)<0
                meas_phiSC_simw_null(ii,2) = meas_phiSC_simw_null(ii,2) + 180;
            end
            if meas_phiSC_simw_null(ii,3)<0
                meas_phiSC_simw_null(ii,3) = meas_phiSC_simw_null(ii,3) + 180;
            end
        end
    end
end

%--------------------------------------------------------------------------
% colormap for color-coding based on phi
if cmap_phi_ind~=0 && any( strcmp(station_shift,staname_split) )
    cmap_phi = [cmap_phi(92:181,:); cmap_phi(1:91,:)];
end



%==========================================================================
%% Figure I - BAZ variation of splitting parameters
%==========================================================================
% model types together and model types spereated

for mt = 1:1:length(modtyp_str_all)

    modtyp_str = modtyp_str_all{mt};

    if strcmp(modtyp_str,'All')==1
        models_2plot = models_sort;
    elseif strcmp(modtyp_str,'H1')==1
        models_2plot = modsH1;
    elseif strcmp(modtyp_str,'H2')==1
        models_2plot = modsH2;
    elseif strcmp(modtyp_str,'T1')==1
        models_2plot = modsT1;
    end


    fig_splitpara = figure('visible','off');

%==========================================================================
    % subplot - panel 1 - phi
    s1 = subplot(2,1,1);
    box on, grid on
    hold on

    %......................................................................
    % plot model range only if not full range is used
    if modrange_low~=0 && modrange_upp~=360

        % not considered BAZ range in background gray
        xdir = [0 0 modrange_low modrange_low];
        ydir = [-90 90 90 -90];
        if any( strcmp(station_shift,staname_split) )
            ydir = [0 180 180 0];
        end
        p1 = patch(xdir,ydir,'white');
        set(p1,'facecolor',modrange_col,'edgecolor',modrange_edcol)

        xdir = [modrange_upp modrange_upp 360 360];
        ydir = [-90 90 90 -90];
        if any( strcmp(station_shift,staname_split) )
            ydir = [0 180 180 0];
        end
        p1 = patch(xdir,ydir,'white');
        set(p1,'facecolor',modrange_col,'edgecolor',modrange_edcol)
    end

    %......................................................................
    % plot theoretical BAZ curves
    if cmap_rms_ind~=0 % color-coding based on RMSE
        for ii = plot_mod_max:-1:1
             plot(BAZ, models_2plot(ii).phi_eff, ...
                  'linewidth',lw_mod, 'color',cmap_rms_sel(ii,:))
        end
        % highligth best-fit model i. e. #1
        plot(BAZ, models_2plot(1).phi_eff, ...
             'linewidth',lw_mod, 'color',colmod_bf_1)
    elseif cmap_rms_ind==0
        % plot best-fit models #2 - #max
        for ii = 2:1:plot_mod_max
             plot(BAZ, models_2plot(ii).phi_eff, ...
                  'linewidth',lw_mod, 'color',colmod_bf_2max)
        end
        % highligth best-fit model i. e. #1
        plot(BAZ, models_2plot(1).phi_eff, ...
             'linewidth',lw_mod, 'color',colmod_bf_1)
    end

    %......................................................................
    % workaround to plot "axis" on top
    % after real axis was set to bottom layer
    % to allow symbols plotted on top
    if any( strcmp(station_shift,staname_split) )
        % axis; % gives xmin xmax ymin ymax
        a = [baz2plot_start, baz2plot_end, 0, 180];
    else
        a = [baz2plot_start, baz2plot_end, -90, 90];
    end

    cx = get(gca,'Xcolor'); % color of x axis
    plot([a(1) a(2)], [a(3) a(3)], 'color', cx)
    plot([a(1) a(2)], [a(4) a(4)], 'color', cx)

    cy = get(gca,'Ycolor'); % color of y axis
    plot([a(1) a(1)], [a(3) a(4)], 'color', cy)
    plot([a(2) a(2)], [a(3) a(4)], 'color', cy)

    %......................................................................
    % plot measured values
    % single nulls
    if ~isempty(meas_phiSC_null)
        for ii = 1:1:length(meas_phiSC_null(:,2))
            plot(meas_BAZ_floor_null(ii), meas_phiSC_null(ii,2), ...
                'ok', 'markerfacecolor',color_face_null, ...
                'markersize',marksize_null-2, 'LineWidth',lw_symbols_null, ...
                'markeredgecolor',col_edge_null);
        end
    end
    % simw nulls
    if ~isempty(meas_phiSC_simw_null)
        for ii = 1:1:length(meas_phiSC_simw_null(:,2))
            plot(meas_BAZ_floor_simw_null(ii), meas_phiSC_simw_null(ii,2), ...
                'sk', 'markerfacecolor', color_face_null, ... %o
                'markersize',marksize_null-2, 'LineWidth',lw_symbols_null, ...
                'markeredgecolor',col_edge_null);
        end
    end

    %......................................................................
    % single splits and stack splits or simw splits
    sizeSC = size(meas_phiSC);
    for FF = 1:sizeSC
        if cmap_phi_ind==0 % no color-coding based on phi
            %h1(FF) =
            errorbar(meas_BAZ_floor(FF), meas_phiSC(FF,2), ...
                              abs(meas_phiSC(FF,2)-meas_phiSC(FF,1)), ...
                              abs(meas_phiSC(FF,2)-meas_phiSC(FF,3)), ...
                              'ok', ... % black error bars
                              'markerfacecolor',colfill(meas_phiSC(FF,4),:), ...
                              'markersize',marksize-1, ...
                              'LineWidth',lw_symbols, ...
                              'markeredgecolor',coledge(meas_phiSC(FF,4),:), ...
                              'color',coledge(meas_phiSC(FF,4),:));
        elseif cmap_phi_ind~=0 % color-coding based on phi
            %h1(FF) =
            errorbar(meas_BAZ_floor(FF),meas_phiSC(FF,2), ...
                              abs(meas_phiSC(FF,2)-meas_phiSC(FF,1)), ...
                              abs(meas_phiSC(FF,2)-meas_phiSC(FF,3)), ...
                              'ok', ... % black error bars
                              'markersize',1, ...
                              'markeredgecolor','none');
            if isequal(colfill(meas_phiSC(FF,4),:),singlcol)
                symbol_method = leg_symb_single_split; % single splits
            elseif isequal(colfill(meas_phiSC(FF,4),:),simwcol)
                symbol_method = leg_symb_simw_split; % simw splits
            elseif isequal(colfill(meas_phiSC(FF,4),:),([79 197 104]./256))
                symbol_method = leg_symb_stack_split; % stack splits
            end
            scatter(meas_BAZ_floor(FF),meas_phiSC(FF,2), ...
                    marksize*4.5, ... % MarkerSize
                    meas_phiSC(FF,2),'filled', ... % MarkerFaceColor
                    symbol_method, ...
                    'MarkerEdgeColor',coledge(meas_phiSC(FF,4),:))
            colormap(cmap_phi);
            caxis([-90 90])
            if any( strcmp(station_shift,staname_split) )
                caxis([0 180])
            end
        end
    end % FF

    %......................................................................
    % axis
    ylabel('\phi_a / N\circE','fontsize',fontsize)
    set(gca, 'xticklabel',[])
    set(gca, 'fontsize',fontsize)
    set(gca, 'xtick',0:45:360)
    if baz2plot_end<180
        set(gca, 'xtick',0:30:360)
    else
        set(gca, 'xtick',0:45:360)
    end
    set(gca, 'ytick',-90:45:180)
    set(gca, 'TickLength',[0.01 0.01], 'XMinorTick','on', 'YMinorTick','on')
    set(gca, 'TickDir','out')
    set(gca, 'layer','bottom') % important to plot symbols on top of axis
    %set(gca, 'layer','top') % important to see background grid
    xlim([baz2plot_start baz2plot_end])
    ylim([-90 90])
    if any( strcmp(station_shift,staname_split) )
        ylim([0 180])
    end

    %......................................................................
	% position of label (c); later used for model pararameter plot
	if cmap_rms_ind==0 % no RMSE color-coding
	    pos_label_c = [0.165 0.892 0 0]; % [x_begin y_begin length height]
    elseif cmap_rms_ind~=0 % RMSE color-coding
	    pos_label_c = [0.155 0.892 0 0];
	end

    % plot labels (a), (b)
    annotation('textbox', [0.165 0.895 0 0], ...
               'String','\bf(a)', 'FitBoxToText','on', ...
               'Margin',0, ... % space around text in pixel units
               'HorizontalAlignment','center', ...
               'VerticalAlignment','middle', ...
               'Color','k', 'BackgroundColor','white', 'FaceAlpha',0.7)
    annotation('textbox', [0.165 0.520 0 0], ...
               'String','\bf(b)', 'FitBoxToText','on', ...
               'Margin',0, ...
               'HorizontalAlignment','center', ...
               'VerticalAlignment','middle', ...
               'Color','k', 'BackgroundColor','white', 'FaceAlpha',0.7)

    %......................................................................
    % plot station name
    annotation('textbox', [0.700 0.875 0 0], ...
               'String',staname_split, 'FitBoxToText','on', ...
               'Margin',0, ...
               'HorizontalAlignment','center', ...
               'VerticalAlignment','bottom', ...
               'Color',colmod_bf_1, 'BackgroundColor','white', ...
               'FaceAlpha',0.7)

    %......................................................................
    % plot model type
    annotation('textbox', [0.785 0.875 0 0], ...
               'String',modtyp_str, 'FitBoxToText','on', 'Margin',0, ...
               'HorizontalAlignment','center', ...
               'VerticalAlignment','bottom', ...
               'Color',colmod_bf_1, 'BackgroundColor','white', ...
               'FaceAlpha',0.7)

    %......................................................................
    % plot RMSE values

    % standard for plot backazimuth 0°-360°
    str_rms = ['RMSE_{tot} = ' num2str(models_2plot(1).RMSE,'%4.2f'), ...
               ', RMSE_{\phi} = ' num2str(models_2plot(1).RMSE_phi,'%4.2f') '\circ', ...
               ', RMSE_{\delta{\itt}} = ' num2str(models_2plot(1).RMSE_dt,'%4.2f') ' s'];
    annotation('textbox', [0.520 0.486 0 0], ...
               'String',str_rms, 'FitBoxToText','on', 'Margin',0, ...
               'HorizontalAlignment','center', ...
               'VerticalAlignment','bottom', ...
               'Color',colmod_bf_1, 'BackgroundColor','white', ...
               'FaceAlpha',0.7)

%     % Ritter et al. 2022, Journal of Seismology
%     % plot backazimuth 0°-110°
%     str_rms = {['RMSE_{tot} = ' num2str(models_2plot(1).RMSE,'%4.2f')], ...
%                ['RMSE_{\phi} = ' num2str(models_2plot(1).RMSE_phi,'%4.2f') '\circ'], ...
%                ['RMSE_{\delta{\itt}} = ' num2str(models_2plot(1).RMSE_dt,'%4.2f') ' s']};
%     annotation('textbox', [0.140 0.595 0 0], ...
%                'String',str_rms, 'FitBoxToText','on', 'Margin',1, ...
%                'HorizontalAlignment','left', ...
%                'VerticalAlignment','bottom', ...
%                'Color',colmod_bf_1, 'BackgroundColor','white', ...
%                'FaceAlpha',0.7)


%==========================================================================
    % subplot - panel 2 - dt
    s2 = subplot(2,1,2);
    hold on
    box on, grid on

    %......................................................................
     % plot model range only if not full range is used
    if modrange_low~=0 && modrange_upp~=360

        % not considered BAZ range in background gray
        xdir = [0 0 modrange_low modrange_low];
        ydir = [-90 90 90 -90];
        if any( strcmp(station_shift,staname_split) )
            ydir = [0 180 180 0];
        end
        p1 = patch(xdir,ydir,'white');
        set(p1,'facecolor',modrange_col,'edgecolor',modrange_edcol)

        xdir = [modrange_upp modrange_upp 360 360];
        ydir = [-90 90 90 -90];
        if any( strcmp(station_shift,staname_split) )
            ydir = [0 180 180 0];
        end
        p1 = patch(xdir,ydir,'white');
        set(p1,'facecolor',modrange_col,'edgecolor',modrange_edcol)

    end

    %......................................................................
    % plot theoretical BAZ curves
    if cmap_rms_ind~=0
        for ii = plot_mod_max:-1:2
             plot(BAZ, models_2plot(ii).dt_eff, ...
                  'linewidth',lw_mod, 'color',cmap_rms_sel(ii,:))
        end
        % highligth best-fit model
        plot(BAZ, models_2plot(1).dt_eff, ...
             'linewidth',lw_mod, 'color',colmod_bf_1)
    elseif cmap_rms_ind==0
        for ii = 2:1:plot_mod_max
             plot(BAZ, models_2plot(ii).dt_eff, ...
                  'linewidth',lw_mod, 'color',colmod_bf_2max)
        end
        % highlight best-fit model
        plot(BAZ, models_2plot(1).dt_eff, ...
             'linewidth',lw_mod, 'color',colmod_bf_1)
    end

    %......................................................................
    % plot measured values
    % singel splits and stack splits or simw splits
    for FF = 1:1:sizeSC
        if cmap_phi_ind==0 % no color-coding based on phi
            %h2(FF) =
            errorbar(meas_BAZ_floor(FF),meas_dtSC(FF,2), ...
                              abs( meas_dtSC(FF,2)-meas_dtSC(FF,1) ), ...
                              abs( meas_dtSC(FF,2)-meas_dtSC(FF,3) ), ...
                              'ok', ...
                              'markerfacecolor',colfill(meas_phiSC(FF,4),:), ...
                              'markersize',marksize-1, ...
                              'LineWidth',lw_symbols, ...
                              'markeredgecolor',coledge(meas_phiSC(FF,4),:), ...
                              'color',coledge(meas_phiSC(FF,4),:));
        elseif cmap_phi_ind~=0 % color-coding based on phi
            %h2(FF) =
            errorbar(meas_BAZ_floor(FF),meas_dtSC(FF,2), ...
                              abs(meas_dtSC(FF,2)-meas_dtSC(FF,1)), ...
                              abs(meas_dtSC(FF,2)-meas_dtSC(FF,3)), ...
                              'ok', ...
                              'markersize',1, ...
                              'markeredgecolor','none');
            if isequal(colfill(meas_phiSC(FF,4),:),([66 91 169]./256))
                symbol_method = leg_symb_single_split; % single split
            elseif isequal(colfill(meas_phiSC(FF,4),:),([0.9290 0.6940 0.1250]))
                symbol_method = leg_symb_simw_split; % simw split
            elseif isequal(colfill(meas_phiSC(FF,4),:),([79 197 104]./256))
                symbol_method = leg_symb_stack_split; % stack split
            end
            scatter(meas_BAZ_floor(FF),meas_dtSC(FF,2), ...
                    marksize*4.5, ... % MarkerSize
                    meas_phiSC(FF,2),'filled', ... % MakerFaceColor
                    symbol_method, ...
                    'MarkerEdgeColor',coledge(meas_phiSC(FF,4),:))
            colormap(cmap_phi);
            caxis([-90 90])
            if any( strcmp(station_shift,staname_split) )
                caxis([0 180])
            end
        end
    end

    %......................................................................
    % workaround to plot "axis" on top
    % after real axis was set to bottom layer
    % to allow null symbols / circles for nulls plotted on top ;)
    % axis; % xmin xmax ymin ymax
    a = [baz2plot_start, baz2plot_end, 0, 4];

    cx = get(gca,'Xcolor'); % color of x axis
    plot([a(1) a(2)], [a(3) a(3)], 'color', cx)
    plot([a(1) a(2)], [a(4) a(4)], 'color', cx)

    cy = get(gca,'Ycolor'); % color of y axis
    plot([a(1) a(1)], [a(3) a(4)], 'color', cy)
    plot([a(2) a(2)], [a(3) a(4)], 'color', cy)

    %......................................................................
    % plot measured nulls, delay times are set manually to zero
    % single nulls
    if ~isempty(meas_dtSC_null)
        for ii = 1:1:length(meas_dtSC_null(:,2))
            plot(meas_BAZ_floor_null(ii),0,'o', ...
                 'markerfacecolor',color_face_null, ...
                 'markersize',marksize_null-2, ...
                 'LineWidth',lw_symbols_null, ...
                 'markeredgecolor',col_edge_null)
        end
    end
    % simw nulls
    if ~isempty(meas_dtSC_simw_null)
        for ii = 1:1:length(meas_dtSC_simw_null(:,2))
            plot(meas_BAZ_floor_simw_null(ii), 0, 's', ...
                 'markerfacecolor',color_face_null, ...
                 'markersize',marksize_null-2, ...
                 'LineWidth',lw_symbols_null, ...
                 'markeredgecolor',col_edge_null)
        end
    end

    %......................................................................
    % colorbar for phi color-coding
    if cmap_phi_ind~=0 & cbar_phi_ind==1
        cb = colorbar('Position', [0.85 0.21 0.025 0.715], ...
                      'TickDirection','out'); % [left bottom width height]
        cb.FontSize = 8; % fontsize of TickLabels (not Label)
        cb.Label.String = '\phi_a / N\circE';
        if any( strcmp(station_shift,staname_split) )
            cb_label_pos_y = -15;
        else
            cb_label_pos_y = -105;
        end
        cb.Label.Position = [0.70 cb_label_pos_y]; % below vertical colorbar
        cb.Label.Rotation = 0; % horizontal colorbar label text
        if any( strcmp(station_shift,staname_split) )
            cb.Ticks = 0:30:180;
        else
            cb.Ticks = -90:30:90;
        end
    end

    %......................................................................
    xlabel('backazimuth / \circ', 'fontsize',fontsize)
    ylabel('\delta{\itt}_a / s', 'fontsize',fontsize)
    set(gca,'fontsize',fontsize)
    if baz2plot_end<180
        set(gca,'xtick',0:30:360)
    else
        set(gca,'xtick',0:45:360)
    end
    set(gca, 'ytick',0:1:4)
    set(gca, 'TickLength',[0.01 0.01], 'XMinorTick','on', 'YMinorTick','on')
    set(gca, 'TickDir','out')
    set(gca, 'layer','bottom') % important to plot symbols on top of axis
    %set(gca, 'layer','top') % important to see background grid
    xlim([baz2plot_start baz2plot_end])
    ylim([0 4])


%==========================================================================
    % create legend

    % are single splits and multi splits in struct
    wcols = unique(colfill(meas_phiSC(:,4),:),'rows');
    swcols = size(wcols);

%--------------------------------------------------------------------------
    % no phi color-coding
    if cmap_phi_ind==0

    %......................................................................
        % BOTH, single splits AND multi splits in dataset
        if swcols(1)==2 && isequal(wcols(1,:),singlcol)

            % multi splits = simw splits
            if isequal(wcols(1,:),simwcol) || isequal(wcols(2,:),simwcol)
                if isequal(wcols(1,:),singlcol)
                    l1 = plot(-10, -10, 'ok', ...
                              'markerfacecolor',wcols(1,:), ...
                              'markersize',marksize-1, ...
                              'LineWidth',lw_symbols, ...
                              'markeredgecolor',[0 0 0]);
                    l2 = plot(-10, -10, 'ok', ...
                              'markerfacecolor',wcols(2,:), ...
                              'markersize',marksize-1, ...
                              'LineWidth',lw_symbols, ...
                              'markeredgecolor',[0 0 0]);
                    handvec = [l1(1), l2(1)];
                    strvec = {'single split', 'simw split'};
                else
                    l1 = plot(-10, -10, 'ok', ...
                              'markerfacecolor',wcols(2,:), ...
                              'markersize',marksize-1, ...
                              'LineWidth',lw_symbols, ...
                              'markeredgecolor',[0 0 0]);
                    l2 = plot(-10, -10, 'ok', ...
                              'markerfacecolor',wcols(1,:), ...
                              'markersize',marksize-1, ...
                              'LineWidth',lw_symbols, ...
                              'markeredgecolor',[0 0 0]);
                    handvec = [l2(1), l1(1)];
                    strvec = {'single split', 'simw split'};
                end
            % multi splits = stack splits
            else
                % Ritter et al. 2022 Journal of Seismology
                % no single splits in NE quadrant after stack,
                % so no entry in legend for single splits
                % and only entry for stack splits
                if strcmp(staname_split,'BFO') && isequal(baz2plot_end,100)
                    l1 = plot(-10, -10, 'ok', ...
                              'markerfacecolor',wcols(2,:), ...
                              'markersize',marksize-1, ...
                              'LineWidth',lw_symbols, ...
                              'markeredgecolor',[0 0 0]);
                    handvec = l1(1);
                    strvec = {'stack split'};
                elseif isequal(wcols(1,:),singlcol)
                    l1 = plot(-10, -10, 'ok', ...
                              'markerfacecolor',wcols(1,:), ...
                              'markersize',marksize-1, ...
                              'LineWidth',lw_symbols, ...
                              'markeredgecolor',[0 0 0]);
                    l2 = plot(-10, -10, 'ok', ...
                              'markerfacecolor',wcols(2,:), ...
                              'markersize',marksize-1, ...
                              'LineWidth',lw_symbols, ...
                              'markeredgecolor',[0 0 0]);
                    handvec = [l1(1), l2(1)];
                    strvec = {'single split', 'stack split'};
                else
                    l1 = plot(-10, -10, 'ok', ...
                              'markerfacecolor',wcols(2,:), ...
                              'markersize',marksize-1, ...
                              'LineWidth',lw_symbols, ...
                              'markeredgecolor',[0 0 0]);
                    l2 = plot(-10, -10, 'ok', ...
                              'markerfacecolor',wcols(1,:), ...
                              'markersize',marksize-1, ...
                              'LineWidth',lw_symbols, ...
                              'markeredgecolor',[0 0 0]);
                    handvec = [l2(1), l1(1)];
                    strvec = {'single split', 'stack split'};
                end
            end

    %......................................................................
        % ONLY single splits OR multi splits in dataset
        else
            if isequal(wcols(1,:),singlcol)
                l1 = plot(-10, -10, 'ok', ...
                          'markerfacecolor',wcols(1,:), ...
                          'markersize',marksize-1, ...
                          'LineWidth',lw_symbols, ...
                          'markeredgecolor',[0 0 0]);
                handvec = l1(1);
                strvec = {'single split'};
            elseif isequal(wcols(1,:),simwcol)
                l1 = plot(-10, -10, 'ok', ...
                          'markerfacecolor',wcols(1,:), ...
                          'markersize',marksize-1, ...
                          'LineWidth',lw_symbols, ...
                          'markeredgecolor',[0 0 0]);
                handvec = l1(1);
                strvec = {'simw split'};
            else
                l1 = plot(-10, -10, 'ok', ...
                          'markerfacecolor',wcols(1,:), ...
                          'markersize',marksize-1, ...
                          'LineWidth',lw_symbols, ...
                          'markeredgecolor',[0 0 0]);
                handvec = l1(1);
                strvec = {'stack split'};
            end
        end

%--------------------------------------------------------------------------
    elseif cmap_phi_ind~=0 % phi color-coding

    %......................................................................
        % BOTH, single splits AND multi splits in dataset
        if swcols(1)==2 && isequal(wcols(1,:),singlcol)

            % multi splits = simw splits
            if isequal(wcols(1,:),simwcol) || isequal(wcols(2,:),simwcol)
                if isequal(wcols(1,:),singlcol)
                    l1 = plot(-10, -10, leg_symb_single_split, ...
                              'markerfacecolor', leg_symb_col_split, ...
                              'markersize',marksize-1, ...
                              'LineWidth',lw_symbols, ...
                              'markeredgecolor',[0 0 0]);
                    l2 = plot(-10, -10, leg_symb_simw_split, ...
                              'markerfacecolor', leg_symb_col_split, ...
                              'markersize',marksize-1, ...
                              'LineWidth',lw_symbols, ...
                              'markeredgecolor',[0 0 0]);
                    handvec = [l1(1), l2(1)];
                    strvec = {'single split', 'simw split'};
                else
                    l1 = plot(-10, -10, leg_symb_simw_split, ...
                              'markerfacecolor', leg_symb_col_split, ...
                              'markersize',marksize-1, ...
                              'LineWidth',lw_symbols, ...
                              'markeredgecolor',[0 0 0]);
                    l2 = plot(-10, -10, leg_symb_single_split, ...
                              'markerfacecolor', leg_symb_col_split, ...
                              'markersize',marksize-1, ...
                              'LineWidth',lw_symbols, ...
                              'markeredgecolor',[0 0 0]);
                    handvec = [l2(1), l1(1)];
                    strvec = {'single split', 'simw split'};
                end

            % multi splits = stack splits
            else
                % Ritter et al. 2022 Journal of Seismology
                % no single splits in NE quadrant after stack,
                % so no entry in legend for single splits
                % and only entry for stack splits
                if strcmp(staname_split,'BFO') && isequal(baz2plot_end,100)
                    l1 = plot(-10, -10, leg_symb_stack_split, ...
                              'markerfacecolor',leg_symb_col_split, ...
                              'markersize',marksize-1, ...
                              'LineWidth',lw_symbols, ...
                              'markeredgecolor',[0 0 0]);
                    handvec = l1(1);
                    strvec = {'stack split'};
                elseif isequal(wcols(1,:),singlcol)
                    l1 = plot(-10, -10, leg_symb_single_split, ...
                              'markerfacecolor',leg_symb_col_split, ...
                              'markersize',marksize-1, ...
                              'LineWidth',lw_symbols, ...
                              'markeredgecolor',[0 0 0]);
                    l2 = plot(-10, -10, leg_symb_stack_split, ...
                              'markerfacecolor',leg_symb_col_split, ...
                              'markersize',marksize-1, ...
                              'LineWidth',lw_symbols, ...
                              'markeredgecolor',[0 0 0]);
                    handvec = [l1(1), l2(1)];
                    strvec = {'single split', 'stack split'};
                else
                    l1 = plot(-10, -10, leg_symb_stack_split, ...
                              'markerfacecolor',leg_symb_col_split, ...
                              'markersize',marksize-1, ...
                              'LineWidth',lw_symbols, ...
                              'markeredgecolor',[0 0 0]);
                    l2 = plot(-10, -10, leg_symb_single_split, ...
                              'markerfacecolor',leg_symb_col_split, ...
                              'markersize',marksize-1, ...
                              'LineWidth',lw_symbols, ...
                              'markeredgecolor',[0 0 0]);
                    handvec = [l2(1), l1(1)];
                    strvec = {'single split', 'stack split'};
                end
            end

    %......................................................................
        % ONLY single splits OR multi splits in dataset
        else
            if isequal(wcols(1,:),singlcol)
                l1 = plot(-10, -10, leg_symb_single_split, ...
                          'markerfacecolor',leg_symb_col_split, ...
                          'markersize',marksize-1, 'LineWidth',lw_symbols, ...
                          'markeredgecolor',[0 0 0]);
                handvec = l1(1);
                strvec = {'single split'};
            elseif isequal(wcols(1,:),simwcol)
                l1 = plot(-10, -10, leg_symb_simw_split, ...
                          'markerfacecolor',leg_symb_col_split, ...
                          'markersize',marksize-1, 'LineWidth',lw_symbols, ...
                          'markeredgecolor',[0 0 0]);
                handvec = l1(1);
                strvec = {'simw split'};
            else
                l1 = plot(-10, -10, leg_symb_stack_split, ...
                          'markerfacecolor',leg_symb_col_split, ...
                          'markersize',marksize-1, 'LineWidth',lw_symbols, ...
                          'markeredgecolor',[0 0 0]);
                handvec = l1(1);
                strvec = {'stack split'};
            end
        end
    end

%--------------------------------------------------------------------------
    % add nulls to legend if in dataset

    % single nulls
    if ~isempty(meas_dtSC_null)
        n1 = plot(-10, -10, leg_symb_single_null, ...
                  'markersize',marksize_null-2, 'LineWidth',lw_symbols_null, ...
                  'markerfacecolor',color_face_null, ...
                  'markeredgecolor',col_edge_null);
        handvec(end+1) = n1;
        strvec{end+1} = 'single null';
    end

    % simw nulls
    if ~isempty(meas_dtSC_simw_null)
        n1 = plot(-10, -10, leg_symb_simw_null, ...
                  'markersize',marksize_null-2, 'LineWidth',lw_symbols_null, ...
                  'markerfacecolor',color_face_null, ...
                  'markeredgecolor',col_edge_null);
        handvec(end+1) = n1;
        strvec{end+1} = 'simw null';
    end

%--------------------------------------------------------------------------
    % plot legend

    h_leg_split = legend(handvec, strvec, 'Orientation','vertical');

    %......................................................................
    % position, [left bottom width height]

    if strcmp({staname_split},'BFO') && baz2plot_end==100
        % BAZ range to plot 0 deg to 100 deg
        % -> Ritter et al. 2022 Journal of Seismology
        h_leg_split.Position = [0.195 0.360 0.065 0.1];
    else
        % BAZ range to plot 0 deg to 360 deg
        % -> Fröhlich et al. 2022 Geophysical Journal International
        h_leg_split.Position = [0.695 0.610 0.065 0.1];
    end

    %......................................................................
    % distance symbole to text, default [30,30]
    h_leg_split.ItemTokenSize(1) = 20;

    %......................................................................
    % (semi-)trancparency of bakground color of legend
    % https://de.mathworks.com/matlabcentral/answers/320127-setting-transparancy-of-legend
    % last access 2021 Aug 09
    h_leg_split.BoxFace.ColorType = 'truecoloralpha';
    h_leg_split.BoxFace.ColorData = uint8(255*[1 1 1 0.70]');



%==========================================================================
    % adjust figure size and position
    pos = get(s1,'Position');
    set(s1,'Position',[pos(1) pos(2)+0 pos(3)*0.9 pos(4)]);

    pos = get(s2,'Position');
    set(s2,'Position',[pos(1) pos(2)+0.1 pos(3)*0.9 pos(4)]);


%==========================================================================
    % save plot
    filename = ['PLOT_SWSmod_' staname_split '_' ...
                num2str(plot_mod_max) 'modelsBaz' modtyp_str '_' ...
                cmap_rms_str '_' cmap_phi_str cbar_phi_str ...
                fit_str '_' data_used ...
                '_baz' num2str(modrange_low) 'to' num2str(modrange_upp) ...
                '_domper' num2str(domper) 's'];

    if strcmp(status_save,'yes')==1
        if vers_out==1 % R2020a or higher
            exportgraphics(fig_splitpara, [filename '.png'], ...
                            'Resolution',myresolution)
            exportgraphics(fig_splitpara, [filename '.eps'], ...
                            'ContentType','vector')
            exportgraphics(fig_splitpara, [filename '.pdf'], ...
                            'ContentType','vector')
        else % below R2020a
            saveas(fig_splitpara, [filename '.png'])
            print([filename '.eps'],'-depsc','-painters')
            print('-dpdf', '-painters', '-r600', [filename '.pdf'])
        end
        % if you work on a Linux machine, you can also try:
        % print('-depsc','-painters','-r600',[filename '.eps'])
        % dir_eps_file = dir([filename '.eps']);
        % [status,cmdout] = system(['epstopdf ' dir_eps_file.name]);
    end


end % mt



%==========================================================================
%% Figure II - statistical overview of model types
%==========================================================================

modtype1 = {models_sort.mod_type};
idx1 = find( strcmp(modtype1,'dipping') );
idx2 = find( strcmp(modtype1,'two_layers') );
idx3 = find( strcmp(modtype1,'single_layer') );

color_H1 = [0.9290 0.6940 0.1250];
color_H2 = [0.6350 0.0780 0.1840];
color_T1 = [0 0.4470 0.7410];
models_sort(1).color = 0;

for kk = 1:length(idx1)
   models_sort(idx1(kk)).color = color_T1;
end
for kk = 1:length(idx2)
   models_sort(idx2(kk)).color = color_H2;
end
for kk = 1:length(idx3)
   models_sort(idx3(kk)).color = color_H1;
end


fig_stat = figure('visible','off');


%==========================================================================
% subplot - panel 1 - models sorted by RMSE

%s11 =
subplot(2,1,1);
hold on
box on, grid on

for ii = 1:length(models_sort)
    plot(ii,models_sort(ii).RMSE, ...
        'color',models_sort(ii).color, ...
        'marker','o', ...
        'markerfacecolor',models_sort(ii).color, ...
        'markersize',10)
    hold on
end

%--------------------------------------------------------------------------
% axis
xlim([0,length(models_sort)])
ylim([0,max([models_sort.RMSE])])
xlabel('worst \leftarrow sorted models \rightarrow best', 'fontsize',fontsize)
ylabel(cb_label_rmse, 'fontsize',fontsize)
set(gca, 'XDir','reverse')
set(gca, 'TickLength',[0.01 0.01], 'XMinorTick','on', 'YMinorTick','on')
set(gca, 'TickDir','out')
set(gca, 'fontsize',fontsize)

%--------------------------------------------------------------------------
% title
if vers_out==1 % MATLAB R2020a or higher
    sgtitle([num2str(length(models_sort)) ' best models'], 'fontsize',fontsize)
else
    title([num2str(length(models_sort)) ' best models'], 'fontsize',fontsize)
end

%--------------------------------------------------------------------------
% plot station name
text(0.90, 0.93, staname_split, ...
     'Units','normalized', ...
     'HorizontalAlignment','left', 'VerticalAlignment','top', ...
     'fontsize',fontsize_RMSE, ...
     'backgroundcolor','w', 'edgecolor','k', 'color',colmod_bf_1);

%--------------------------------------------------------------------------
% legend for model types
ll1 = plot(-10, -10, 'ok', 'markersize',10, ...
           'markeredgecolor',color_T1, 'markerfacecolor',color_T1);
ll2 = plot(-10, -10, 'ok', 'markersize',10, ...
           'markeredgecolor',color_H2, 'markerfacecolor',color_H2);
ll3 = plot(-10, -10, 'ok', 'markersize',10, ...
           'markeredgecolor',color_H1, 'markerfacecolor',color_H1);

h_leg_distri = legend([ll1,ll2,ll3], ...
	{'1 dipping layer (T1)','2 horizontal layers (H2)', ...
	 '1 horizontal layer (H1)'}, ...
	 'Location','southwest');

% distance symbole to text, default [30,30]
h_leg_distri.ItemTokenSize(1) = 20;

% (semi-)trancparency of bakground color of legend
% https://de.mathworks.com/matlabcentral/answers/320127-setting-transparancy-of-legend
% last access 2021 Aug 09
h_leg_distri.BoxFace.ColorType = 'truecoloralpha';
h_leg_distri.BoxFace.ColorData = uint8(255*[1 1 1 0.70]');


%==========================================================================
% suplot - panel 2 - bar plot to show distribution of model types

%s22 =
subplot(2,1,2);
hold on
box on, grid on

b1 = bar(1,length(idx1));
b2 = bar(2,length(idx2));
b3 = bar(3,length(idx3));

set(b1, 'FaceColor',color_T1)
set(b2, 'FaceColor',color_H2)
set(b3, 'FaceColor',color_H1)

ylim([0 keep_mods])
ylabel('model count', 'fontsize',fontsize)
set(gca, 'xticklabel',[])
set(gca, 'xtick',[])
set(gca, 'TickDir','out')
set(gca, 'fontsize',fontsize)


%==========================================================================
% save plot
filename = ['PLOT_SWSmod_' staname_split '_' num2str(keep_mods) ...
            'modelsDist_' fit_str '_' data_used ...
            '_baz' num2str(modrange_low) 'to' num2str(modrange_upp) ...
            '_domper' num2str(domper) 's'];

if strcmp(status_save,'yes')==1
    if vers_out==1 % R2020a or higher
        exportgraphics(fig_stat, [filename '.png'], ...
                        'Resolution',myresolution)
        exportgraphics(fig_stat, [filename '.eps'], ...
                        'ContentType','vector')
        exportgraphics(fig_stat, [filename '.pdf'], ...
                        'ContentType','vector')
    else % below R2020a
        saveas(fig_stat, [filename '.png'])
        print([filename '.eps'], '-depsc', '-painters')
        print('-dpdf', '-painters', '-r600', [filename '.pdf'])
    end
    % if you work on a Linux machine, you can also try:
    % print('-depsc','-painters','-r600',[filename '.eps'])
    % dir_eps_file = dir([filename '.eps']);
    % [status,cmdout] = system(['epstopdf ' dir_eps_file.name]);
end




%==========================================================================
%% Figure III - model parameter distribution of model types seperatly
%==========================================================================

tick_label_step = 1;
tick_step = tick_label_step;

tick_pos_mid = transpose(0.5:tick_step:plot_mod_max);
tick_pos_low = transpose(0:1:plot_mod_max-1);

tick_H1_label = cell(plot_mod_max);
tick_H2_label = cell(plot_mod_max);
tick_T1_label = cell(plot_mod_max);


%==========================================================================
fig_H1 = figure('visible','off'); % one layer
    box on, grid on
    hold on

%--------------------------------------------------------------------------
    % RMSE color-coding
    if cmap_rms_ind~=0
        for ii = plot_mod_max:-1:1
            plot(dt_H1(ii),phi_H1(ii), 'o', ...
                 'MarkerSize',mymarkersize_symbols, ...
                 'LineWidth',mylinewidth_symbols, ...
                 'color',cmap_rms_sel(ii,:))
        end

        colormap(cmap_rms_sel)

        if strcmp(tick_label_pos_str,'low')
        % horizontal lines between color sections
            for tl = 1:tick_label_step:plot_mod_max-1
                tick_H1_label{tl} = num2str(round(RMSE_H1(tl),4), '%.4f');
                tick_H1_label{tl+1} = ' ';
            end
            cb = colorbar('Location','eastoutside', ...
                          'Ticks',tick_pos_low, ...
                          'TickLabels',tick_H1_label, ...
                          'TickDirection','in', ...
                          'TickLength',0.06);
        elseif strcmp(tick_label_pos_str,'mid')
        % tick & ticklabel in the middle of the corresponding color section
            tick_H1_label_temp = RMSE_H1(1:tick_label_step:plot_mod_max);
            tick_H1_label = num2str(transpose(round(tick_H1_label_temp,4)), ...
                                    '%.4f');
            cb = colorbar('Location','eastoutside', ...
                          'Ticks',tick_pos_mid, ...
                          'TickLabels',tick_H1_label, ...
                          'TickDirection','out', ...
                          'TickLength',0);
        end
        cb.FontSize = 8; % fontsize of TickLabels (not cb.Label)
        cb.Label.String = cb_label_rmse;
        cb.Label.Position = [0.70 -0.2 0]; % below vertical colorbar
        cb.Label.Rotation = 0; % horizontal
        colormap(cmap_rms_sel)
        caxis([0 plot_mod_max])

        % highligth best-fit model
        plot(dt_H1(1),phi_H1(1), 'o', ...
             'MarkerSize',mymarkersize_symbols/2, ...
             'LineWidth',mylinewidth_symbols, ...
             'color',colmod_bf_1)

    %......................................................................
    % no RMSE color-coding
    elseif cmap_rms_ind==0
        plot(dt_H1,phi_H1, 'o', ...
             'MarkerSize',mymarkersize_symbols, ...
             'LineWidth',mylinewidth_symbols, ...
             'color',colmod_bf_2max)
        % highligth best-fit model
        plot(dt_H1(1),phi_H1(1), 'o', ...
             'MarkerSize',mymarkersize_symbols, ...
             'LineWidth',mylinewidth_symbols, ...
             'color',colmod_bf_1)
    end

%--------------------------------------------------------------------------
    % axis
    xlim([0,4])
    ylim([-90,90])
    set(gca, 'ytick',-90:30:90)
    set(gca, 'TickLength',[0.01 0.01], 'XMinorTick','on', 'YMinorTick','on')
    set(gca, 'TickDir','out')
    set(gca, 'fontsize',fontsize)
    xlabel('delay time \delta{\itt} / s', 'fontsize',fontsize)
    ylabel('fast polarization direction \phi / N\circE', 'fontsize',fontsize)

%--------------------------------------------------------------------------
    % plot label (c)
    annotation('textbox', pos_label_c, ...
               'String','\bf(c)', 'FitBoxToText','on', ...
               'Margin',0, ...
               'HorizontalAlignment','center', ...
               'VerticalAlignment','middle', ...
               'Color','k', 'BackgroundColor','white', 'FaceAlpha',0.7)

%--------------------------------------------------------------------------
    % plot station name
    text(0.80, 0.97, staname_split, ...
         'Units','normalized', ...
         'HorizontalAlignment','left', ...
         'VerticalAlignment','top', ...
         'fontsize',fontsize_RMSE, ...
         'backgroundcolor','w', 'edgecolor','k', 'color',colmod_bf_1)

%--------------------------------------------------------------------------
    % plot model type
    text(0.93, 0.97, 'H1', ...
         'Units','normalized', ...
         'HorizontalAlignment','left', 'VerticalAlignment','top', ...
         'fontsize',fontsize_RMSE, ...
         'backgroundcolor','w', 'edgecolor','k', 'color',colmod_bf_1)

%--------------------------------------------------------------------------
    % save plot
    filename  = ['PLOT_SWSmod_' staname_split '_' ...
                 num2str(keep_mods_sep) 'modelsParaH1_' ...
                 cmap_rms_str '_' fit_str '_' data_used ...
                 '_baz' num2str(modrange_low) 'to' num2str(modrange_upp) ...
                 '_domper' num2str(domper) 's'];

    if strcmp(status_save,'yes')==1
        if vers_out==1 % R2020a or higher
            exportgraphics(fig_H1, [filename '.png'], ...
                            'Resolution',myresolution)
            exportgraphics(fig_H1, [filename '.eps'], ...
                            'ContentType','vector')
            exportgraphics(fig_H1, [filename '.pdf'], ...
                            'ContentType','vector')
        else % below R2020a
            saveas(fig_H1, [filename '.png'])
            print([filename '.eps'], '-depsc', '-painters')
            print('-dpdf', '-painters', '-r600', [filename '.pdf'])
        end
        % if you work on a Linux machine, you can also try:
        % print('-depsc','-painters','-r600',[filename '.eps'])
        % dir_eps_file = dir([filename '.eps']);
        % [status,cmdout] = system(['epstopdf ' dir_eps_file.name]);
    end


%==========================================================================
fig_H2 = figure('visible','off'); % two layers
    box on, grid on
    hold on

%--------------------------------------------------------------------------
%   % plot polyons around clusters
%   % Ritter et al. 2022 Journal of Seismology

%     % BFO
%     if strcmp(staname_split,"BFO")==1 && ...
%        modrange_upp==100 && ...
%        strcmp(data_used,"multiWS")==1 &&
%
%         x_pol_1 = [0.35 0.35 0.9 1.15 1.15 0.6];
%         y_pol_1 = [ -90  -75  -50 -50  -60 -90];
%         pol_1 = polyshape(x_pol_1,y_pol_1);
%
%         x_pol_2 = [0.9 0.35 0.35 1.2 1.9 1.9 1.4 0.9 0.9 1.9 1.9 0.9];
%         y_pol_2 = [-35  -15   25  45  45  35  35  15 -25 -25 -35 -35] ;
%         pol_2 = polyshape(x_pol_2,y_pol_2);
%
%         x_pol_3 = [ 1 0.40 0.40 0.9 1.7 2.6 2.6  1];
%         y_pol_3 = [50   70   90  90  60  60  50 50];
%         pol_3 = polyshape(x_pol_3,y_pol_3);
%
%         plot(pol_1, 'FaceColor','r', 'EdgeColor','k', 'Facealpha',0.25)
%         plot(pol_2, 'FaceColor','b', 'EdgeColor','k', 'Facealpha',0.25)
%         plot(pol_3, 'FaceColor','r', 'EdgeColor','k', 'Facealpha',0.25)
%
%         disp('BFOpolyons')
%     end

%--------------------------------------------------------------------------
    % RMSE color-coding
    if cmap_rms_ind~=0
        for ii = plot_mod_max:-1:1

%             % Fröhlich et al. 2022 Geophysical Journal International
%             if strcmp(staname_split,"ECH") && ...
%                strcmp(data_used,"multiWS") && phi_H2(2,ii)==90
%                 phi_H2(2,ii) = -90;
%                 disp('switch value')
%             end

            plot(dt_H2(1,ii),phi_H2(1,ii), 's', ...
                 'MarkerSize',mymarkersize_symbols, ...
                 'LineWidth',mylinewidth_symbols, ...
                 'color',cmap_rms_sel(ii,:))
            plot(dt_H2(2,ii),phi_H2(2,ii), 'o', ...
                 'MarkerSize',mymarkersize_symbols, ...
                 'LineWidth',mylinewidth_symbols, ...
                 'color',cmap_rms_sel(ii,:))
        end

        colormap(cmap_rms_sel)

        if strcmp(tick_label_pos_str,'low')
        % horizontal lines between color sections
            for tl = 1:tick_label_step:plot_mod_max-1
                tick_H2_label{tl} = num2str(round(RMSE_H2(tl),4), '%.4f');
                tick_H2_label{tl+1} = ' ';
            end
            cb = colorbar('Location','eastoutside', ...
                          'Ticks',tick_pos_low, ...
                          'TickLabels',tick_H2_label, ...
                          'TickDirection','in', ...
                          'TickLength',0.06);
        elseif strcmp(tick_label_pos_str,'mid')
        % tick & ticklabel in the middle of the corresponding color section
            tick_H2_label_temp = RMSE_H2(1:tick_label_step:plot_mod_max);
            tick_H2_label = num2str(transpose(round(tick_H2_label_temp,4)), ...
                                    '%.4f');
            cb = colorbar('Location','eastoutside', ...
                          'Ticks',tick_pos_mid, ...
                          'TickLabels',tick_H2_label, ...
                          'TickDirection','out', ...
                          'TickLength',0);
        end

        cb.FontSize = 8; % fontsize of TickLabels (not cb.Label)
        cb.Label.String = cb_label_rmse;
        cb.Label.Position = [0.70 -0.2 0];
        cb.Label.Rotation = 0; % horizontal
        colormap(cmap_rms_sel)
        caxis([0 plot_mod_max])

        % highligth best-fit model
        plot(dt_H2(1,1),phi_H2(1,1), 's', ...
             'MarkerSize',mymarkersize_symbols/2, ...
             'LineWidth',mylinewidth_symbols, ...
             'color',colmod_bf_1)
        plot(dt_H2(2,1),phi_H2(2,1), 'o', ...
             'MarkerSize',mymarkersize_symbols/2, ...
             'LineWidth',mylinewidth_symbols, ...
             'color',colmod_bf_1)

    %......................................................................
    % no RMSE color-coding
    elseif cmap_rms_ind==0
        plot(dt_H2(1,:),phi_H2(1,:), 's', ...
             'MarkerSize',mymarkersize_symbols, ...
             'LineWidth',mylinewidth_symbols, ...
             'color',colmod_bf_2max)
        plot(dt_H2(2,:),phi_H2(2,:), 'o', ...
             'MarkerSize',mymarkersize_symbols, ...
             'LineWidth',mylinewidth_symbols, ...
             'color',colmod_bf_2max)

        % highligth best-fit model
        plot(dt_H2(1,1),phi_H2(1,1), 's', ...
             'MarkerSize',mymarkersize_symbols, ...
             'LineWidth',mylinewidth_symbols, ...
             'color',colmod_bf_1)
        plot(dt_H2(2,1),phi_H2(2,1), 'o', ...
             'MarkerSize',mymarkersize_symbols, ...
             'LineWidth',mylinewidth_symbols, ...
             'color',colmod_bf_1)
    end

%--------------------------------------------------------------------------
    % axis
    xlim([0 4])
    ylim([-90 90])
    set(gca, 'ytick',-90:30:90)
    set(gca, 'TickLength',[0.01 0.01], 'XMinorTick','on', 'YMinorTick','on')
    set(gca, 'TickDir','out')
    set(gca, 'fontsize',fontsize)
    xlabel('delay time \delta{\itt} / s', 'fontsize',fontsize)
    ylabel('fast polarization direction \phi / N\circE', 'fontsize',fontsize)

%--------------------------------------------------------------------------
    % plot label (c)
    annotation('textbox', pos_label_c, ...
               'String','\bf(c)', 'FitBoxToText','on', ...
               'Margin',0, ...
               'HorizontalAlignment','center', ...
               'VerticalAlignment','middle', ...
               'Color','k', 'BackgroundColor','white', 'FaceAlpha',0.7)

%--------------------------------------------------------------------------
    % plot station name
    text(0.80, 0.97, staname_split, ...
         'Units','normalized', ...
         'HorizontalAlignment','left','VerticalAlignment','top', ...
         'fontsize',fontsize_RMSE, ...
         'backgroundcolor','w', 'edgecolor','k', 'color',colmod_bf_1)

%--------------------------------------------------------------------------
    % plot modeltype
    text(0.93, 0.97, 'H2', ...
         'Units','normalized', ...
         'HorizontalAlignment','left', 'VerticalAlignment','top', ...
         'fontsize',fontsize_RMSE, ...
         'backgroundcolor','w', 'edgecolor','k', 'color',colmod_bf_1)

%--------------------------------------------------------------------------
    % legend for lower layer and upper layer
    lupp = plot(-10, -10, 'ok', ...
                'markersize',mymarkersize_symbols, ...
                'LineWidth',lw_symbols, ...
                'markerfacecolor','w', 'markeredgecolor',[0 0 0]);
    llow = plot(-10, -10, 'sk', ...
                'markersize',mymarkersize_symbols, ...
                'LineWidth',lw_symbols, ...
                'markerfacecolor','w', 'markeredgecolor',[0 0 0]);
    h_leg_H2 = legend([lupp,llow], {'upper layer','lower layer'}, ...
                        'Orientation','vertical', ...
                        'Location','southeast', ...
                        'fontsize',fontsize_RMSE);

    % distance symbole to text, default [30,30]
    h_leg_H2.ItemTokenSize(1) = 20;

%--------------------------------------------------------------------------
    % save plot
    filename = ['PLOT_SWSmod_' staname_split '_' ...
                num2str(keep_mods_sep) 'modelsParaH2_' ...
                cmap_rms_str '_' fit_str '_' data_used ...
                '_baz' num2str(modrange_low) 'to' num2str(modrange_upp) ...
                '_domper' num2str(domper) 's'];

    if strcmp(status_save,'yes')==1
        if vers_out==1 % R2020a or higher
            exportgraphics(fig_H2, [filename '.png'], ...
                            'Resolution',myresolution)
            exportgraphics(fig_H2, [filename '.eps'], ...
                            'ContentType','vector')
            exportgraphics(fig_H2, [filename '.pdf'], ...
                            'ContentType','vector')
        else % below R2020a
            saveas(fig_H2, [filename '.png'])
            print([filename '.eps'], '-depsc', '-painters')
            print('-dpdf', '-painters', '-r600', [filename '.pdf'])
        end
        % if you work on a Linux machine, you can also try:
        % print('-depsc','-painters','-r600',[filename '.eps'])
        % dir_eps_file = dir([filename '.eps']);
        % [status,cmdout] = system(['epstopdf ' dir_eps_file.name]);
    end


%==========================================================================
% >>> under development <<<
fig_T1 = figure('visible','off'); % one dipping layer
    box on, grid on
    hold on

%--------------------------------------------------------------------------
    % RMSE color-coding
    if cmap_rms_ind~=0
        for ii = plot_mod_max:-1:1
            plot(downdipdir_T1(ii),dip_T1(ii), 'o', ...
                 'MarkerSize',mymarkersize_symbols, ...
                 'LineWidth',mylinewidth_symbols, ...
                 'color',cmap_rms_sel(ii,:))
        end

        colormap(cmap_rms_sel)

        if strcmp(tick_label_pos_str,'low')
            for tl = 1:tick_label_step:plot_mod_max-1
                tick_T1_label{tl} = num2str(round(RMSE_T1(tl),4),'%.4f');
                tick_T1_label{tl+1} = ' ';
            end
            cb = colorbar('Location','eastoutside', ...
                          'Ticks',tick_pos_low, ...
                          'TickLabels',tick_T1_label, ...
                          'TickDirection','in', ...
                          'TickLength',0.06);
        elseif strcmp(tick_label_pos_str,'mid')
        % tick & ticklabel in the middle of the corresponding color section
            tick_T1_label_temp = RMSE_T1(1:tick_label_step:plot_mod_max);
            tick_T1_label =  num2str(transpose(round(tick_T1_label_temp,4)), ...
                                     '%.4f');
            cb = colorbar('Location','eastoutside', ...
                          'Ticks',tick_pos_mid, ...
                          'TickLabels',tick_T1_label, ...
                          'TickDirection','out', ...
                          'TickLength',0);
        end

        cb.FontSize = 8; % fontsize of TickLabels (not cb.Label)
        cb.Label.String = cb_label_rmse;
        cb.Label.Position = [0.70 -0.2 0]; % below vertical colorbar
        cb.Label.Rotation = 0; % horizontal
        colormap(cmap_rms_sel)
        caxis([0 plot_mod_max])

        % highligth best-fit model
        plot(downdipdir_T1(1),dip_T1(1), 'o', ...
             'MarkerSize',mymarkersize_symbols/2, ...
             'LineWidth',mylinewidth_symbols, ...
             'color',colmod_bf_1)

    %......................................................................
    % no RMSE color-coding
    elseif cmap_rms_ind==0
        plot(downdipdir_T1,dip_T1, 'o', ...
             'MarkerSize',mymarkersize_symbols, ...
             'LineWidth',mylinewidth_symbols, ...
             'color',colmod_bf_2max)

        % highligth best-fit model
        plot(downdipdir_T1(1),dip_T1(1), 'o', ...
             'MarkerSize',mymarkersize_symbols, ...
             'LineWidth',mylinewidth_symbols, ...
             'color',colmod_bf_1)
        % scatter(downdipdir_T1,dip_T1, ...
                % thick_T1,colmod_bf_2max,'o', ...
                % 'LineWidth',mylinewidth_symbols)
    end

%--------------------------------------------------------------------------
    % axis
    xlim([0 360])
    ylim([0 90])
    set(gca, 'TickLength',[0.01 0.01], 'XMinorTick','on', 'YMinorTick','on')
    set(gca, 'TickDir','out')
    set(gca, 'fontsize',fontsize)
    xlabel('downdip direction (= strike direction \sigma+90\circ) / N\circE', ...
            'fontsize',fontsize)
    ylabel('dip angle \Psi / \circ', 'fontsize',fontsize)

%--------------------------------------------------------------------------
    % plot label (c)
    annotation('textbox', pos_label_c, ...
               'String','\bf(c)', 'FitBoxToText','on', ...
               'Margin',0, ...
               'HorizontalAlignment','center', ...
               'VerticalAlignment','middle', ...
               'Color','k', 'BackgroundColor','white', 'FaceAlpha',0.7)

%--------------------------------------------------------------------------
    % plot station name
    text(0.80,0.97, staname_split, ...
         'Units','normalized', ...
         'HorizontalAlignment','left', 'VerticalAlignment','top', ...
         'fontsize',fontsize_RMSE, ...
         'backgroundcolor','w', 'edgecolor','k', 'color',colmod_bf_1);

%--------------------------------------------------------------------------
   % plot model type
    text(0.93,0.97, 'T1', ...
         'Units','normalized', ...
         'HorizontalAlignment','left', 'VerticalAlignment','top', ...
         'fontsize',fontsize_RMSE, ...
         'backgroundcolor','w', 'edgecolor','k', 'color',colmod_bf_1);

%--------------------------------------------------------------------------
    % save plot
    filename = ['PLOT_SWSmod_' staname_split '_' ...
                num2str(keep_mods_sep) 'modelsParaT1_' ...
                cmap_rms_str '_' fit_str '_' data_used ...
                '_baz' num2str(modrange_low) 'to' num2str(modrange_upp) ...
                '_domper' num2str(domper) 's'];

    if strcmp(status_save,'yes')==1
        if vers_out==1 % R2020a or higher
            exportgraphics(fig_T1, [filename '.png'], ...
                            'Resolution',myresolution)
            exportgraphics(fig_T1, [filename '.eps'], ...
                            'ContentType','vector')
             exportgraphics(fig_T1, [filename '.pdf'], ...
                            'ContentType','vector')
        else % below R2020a
            saveas(fig_T1, [filename '.png'])
            print([filename '.eps'], '-depsc', '-painters')
            print('-dpdf', '-painters', '-r600', [filename '.pdf'])
        end
        % if you work on a Linux machine, you can also try:
        % print('-depsc','-painters','-r600',[filename '.eps'])
        % dir_eps_file = dir([filename '.eps']);
        % [status,cmdout] = system(['epstopdf ' dir_eps_file.name]);
    end


%==========================================================================
end % EOF
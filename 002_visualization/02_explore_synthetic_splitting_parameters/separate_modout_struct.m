% #############################################################################
% Explore forwardt calculated splitting parameters for structural anisotropy
% -----------------------------------------------------------------------------
% Supported model types
% - One horizontal layer (with HTI): H1
% - Two horizontal layer (with HTI): H2
% - One tilted layer (with TTI): T1
% Previously calculated synthetic splitting parameters
% - https://github.com/yvonnefroehlich/sws-visualization-and-modeling/tree/main/003_modeling
% - The output MATLAB struct is split into separate structs for the different model types
% -----------------------------------------------------------------------------
% History
% - Created: 2025/04/08
% -----------------------------------------------------------------------------
% Versions
%   MATLAB R2024b
% -----------------------------------------------------------------------------
% Contact
% - Author: Yvonne Fröhlich
% - ORCID: https://orcid.org/0000-0002-8566-0619
% - GitHub: https://github.com/yvonnefroehlich/sws-visualization-and-modeling
% -----------------------------------------------------------------------------
% Related to
% - Fröhlich Y. & Ritter J. R. R. (2024)
%   Vertical and Small-scale Lateral Varying Seismic Anisotropy in the Upper Mantle
%   Underneath the Upper Rhine Graben, Central Europe.
%   Annual Meeting of the American Geophysical Union,
%   Division Session Exploring Innovations and New Directions in Seismic Anisotropy and
%   Attenuation: Observations, Models, and Experiments I Oral, DI21A-02.
%   Abstract ID 1578275, https://agu.confex.com/agu/agu24/meetingapp.cgi/Paper/1578275.
%   http://dx.doi.org/10.5281/zenodo.14510993.
% - Fröhlich Y., Grund M. & Ritter J. R. R. (2024)
%   Lateral and vertical variations of seismic anisotropy in the lithosphere-asthenosphere
%   system underneath Central Europe from long-term splitting measurements.
%   Geophysical Journal International, 239(1), 112-135.
%   https://doi.org/10.1093/gji/ggae245.
% #############################################################################


% -----------------------------------------------------------------------------
%%
clear all
close all

folder_data = '01_in_data';
name_struct = 'sws_modout_domper8s';
models_all = load([folder_data '/' name_struct '.mat']);


% -----------------------------------------------------------------------------
%%
models_all_split = models_all.splitmods;

models_H1 = models_all_split(strcmp({models_all_split.type}, 'single_layer'));
models_H2 = models_all_split(strcmp({models_all_split.type}, 'two_layers'));
models_T1 = models_all_split(strcmp({models_all_split.type}, 'dipping'));

models_H1_in = [models_H1.mod_paras];
models_H2_in = [models_H2.mod_paras];
models_T1_in = [models_T1.mod_paras];

models_H1_in_phi = [models_H1_in.phi_in];
models_H1_in_dt = [models_H1_in.dt_in];

models_H2_in_phis = [models_H2_in.phi_in];
models_H2_in_dts = [models_H2_in.dt_in];
models_H2_in_phi1 = models_H2_in_phis(1,:);
models_H2_in_dt1 = models_H2_in_dts(1,:);
models_H2_in_phi2 = models_H2_in_phis(2,:);
models_H2_in_dt2 = models_H2_in_dts(2,:);

models_T1_in_downdipdir = [models_T1_in.downdipdir];
models_T1_in_dip = [models_T1_in.dip];
models_T1_in_thick = [models_T1_in.thick];


% -----------------------------------------------------------------------------
%%
clear model_out
for i_mod = 1:1:length(models_H1)
    model_out(i_mod).phi_eff = models_H1(i_mod).phi_eff;
    model_out(i_mod).dt_eff = models_H1(i_mod).dt_eff;
    model_out(i_mod).phi_in = models_H1_in_phi(i_mod);
    model_out(i_mod).dt_in = models_H1_in_dt(i_mod);
end
disp(model_out)
save([folder_data '/' name_struct '_H1.mat'],'model_out')
% test_H1 = load([folder_data '/' name_struct '_H1.mat']);

%%
clear model_out
for i_mod = 1:1:length(models_H2)
    model_out(i_mod).phi_eff = models_H2(i_mod).phi_eff;
    model_out(i_mod).dt_eff = models_H2(i_mod).dt_eff;
    model_out(i_mod).phis_in = models_H2_in_phis(:,i_mod);
    model_out(i_mod).dts_in = models_H2_in_dts(:,i_mod);
end
disp(model_out)
save([folder_data '/' name_struct '_H2.mat'],'model_out')
% test_H2 = load([folder_data '/' name_struct '_H2.mat']);

%%
clear model_out
for i_mod = 1:1:length(models_T1)
    model_out(i_mod).phi_eff = [models_T1(i_mod).phi_eff];
    model_out(i_mod).dt_eff = [models_T1(i_mod).dt_eff];
    model_out(i_mod).downdipdir_in = models_T1_in_downdipdir(i_mod);
    model_out(i_mod).dip_in = models_T1_in_dip(i_mod);
    model_out(i_mod).thick_in = models_T1_in_thick(i_mod);
end
disp(model_out)
save([folder_data '/' name_struct '_T1.mat'],'model_out')
% test_T1 = load([folder_data '/' name_struct '_T1.mat']);

# ==========================================================================
# Load forward calculated splitting parameters for structural anisotropy
# -----------------------------------------------------------------------------
# Supported model types
# - One horizontal layer (with HTI): H1
# - Two horizontal layers (with HTI): H2
# - One tilted layer (with TTI): T1
# Previously calculated synthetic splitting parameters
# - https://github.com/yvonnefroehlich/sws-visualization-and-modeling/tree/main/003_modeling
# - The output MATLAB struct is split into separate structs for the three model types
# -----------------------------------------------------------------------------
# History
# - Created: 2026/01/28
# -----------------------------------------------------------------------------
# Versions
#   PyGMT v0.18.0 -> https://www.pygmt.org/v0.18.0 | https://www.pygmt.org
#   GMT 6.6.0 -> https://www.generic-mapping-tools.org
# -----------------------------------------------------------------------------
# Contact
# - Author: Yvonne Fröhlich
# - ORCID: https://orcid.org/0000-0002-8566-0619
# - GitHub: https://github.com/yvonnefroehlich/sws-visualization-and-modeling
# -----------------------------------------------------------------------------
# Related to
# - Fröhlich (2025) Dissertation
#   https://doi.org/10.5445/IR/1000183786
# - Fröhlich, Ritter (2024) Annual Meeting of the American Geophysical Union
#   http://dx.doi.org/10.5281/zenodo.14510993
# - Fröhlich, Grund, Ritter (2024) Geophysical Journal International
#   https://doi.org/10.1093/gji/ggae245
# -----------------------------------------------------------------------------
# LICENSE
#
# Copyright (C) 2026  Yvonne Fröhlich (v2.0)
# https://github.com/yvonnefroehlich/sws-visualization-and-modeling
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.
# --------------------------------------------------------------------------
# TERMS OF USE
#
# The loading routines are provided "as is" and without any warranty.
# The author cannot be held responsible for anything that happens to you
# or your equipment. Use it at your own risk.
# --------------------------------------------------------------------------
# CONTRIBUTING
#
# Feel free to modify/adjust the code for your needs. Submit improvements
# and report bugs by opening a "New issue" in the GitHub repository (:
# ==========================================================================



import numpy as np
import pandas as pd
from scipy import io

# %%
def load_models(
    model_type,  ## H1 | H2 | T1
    dom_per=8,  ## 6 | 8 | 10  # in seconds  (TEST data provided for 8 s)
    path_models="000_test_data",
    # "test": provided test data
    # "default": naming structur from forwardt calculation
    # <model_name>: user defined name
    file_models="test",
# -----------------------------------------------------------------------------
    # Limits for model parameters
    # H1
    phi_min=-90,
    phi_max=90,
    dt_min=0,
    dt_max=4,
    # H2 (index 1 lower layer, index 2 upper layer)
    phi1_min=-90,
    phi1_max=90,
    phi2_min=-90,
    phi2_max=90,
    dt1_min=0,
    dt1_max=4,
    dt2_min=0,
    dt2_max=4,
    # T1
    dip_min=0,
    dip_max=90,
    thick_min=0,
    thick_max=700,
    downdipdir_min=0,
    downdipdir_max=360,
):


# %%
# -----------------------------------------------------------------------------
# General stuff
# -----------------------------------------------------------------------------
    if file_models == "test":
        file_models = f"sws_modout_domper{dom_per}s_{model_type}_TEST.mat"
    elif file_models == "default":
        file_models = f"sws_modout_domper{dom_per}s_{model_type}.mat"


# %%
# -----------------------------------------------------------------------------
# Prepare data
# -----------------------------------------------------------------------------
    print(f"Dominant period {dom_per} s - Model type {model_type}")
    models_mat = io.loadmat(f"{path_models}/{file_models}")
    models_dict = models_mat["model_out"][0]
    models_df_raw = pd.DataFrame(models_dict)
    N_total = len(models_df_raw)
    models_df_raw["i_total"] = np.arange(N_total).tolist()

# -----------------------------------------------------------------------------
    models_df = models_df_raw
    phi_in = []
    dt_in = []
    phi1_in = []
    phi2_in = []
    dt1_in = []
    dt2_in = []
    dip_in = []
    thick_in = []
    downdipdir_in = []
    phi_gmt = []
    phi1_gmt = []
    phi2_gmt = []
    phi_T1 = []
    baz_nulls = [0] * N_total
    baz_null1 = []
    baz_null2 = []
    baz_null3 = []
    baz_null4 = []

    for i_model in range(N_total):

        match model_type:
# .............................................................................
            case "H1":
                phi_in_temp = int(str(models_df_raw["phi_in"][i_model][0][0]))
                dt_in_temp = float(str(models_df_raw["dt_in"][i_model][0][0]))
                phi_in.append(phi_in_temp)
                dt_in.append(dt_in_temp)

                if phi_in_temp > 0:
                    phi_gmt_temp = 90 - phi_in_temp
                else:
                    phi_gmt_temp = 90 + (-phi_in_temp)
                phi_gmt.append(phi_gmt_temp)

                # nulls (occur in steps of 90 deg)
                baz_nulls_neg = [
                    phi_in_temp, phi_in_temp + 90, phi_in_temp + 180, phi_in_temp + 270
                ]
                baz_nulls_temp = []
                for baz_null in baz_nulls_neg:
                    baz_null_pos = baz_null
                    if baz_null < 0:
                        baz_null_pos = 360 + baz_null  # -90 to 90 deg
                    baz_nulls_temp.append(baz_null_pos)
# .............................................................................
            case "H2":
                phis_in = np.squeeze(models_df_raw["phis_in"][i_model])
                dts_in = np.squeeze(models_df_raw["dts_in"][i_model])
                phi1_in_temp = phis_in[0]
                phi2_in_temp = phis_in[1]
                dt1_in_temp = dts_in[0]
                dt2_in_temp = dts_in[1]
                phi1_in.append(phi1_in_temp)
                phi2_in.append(phi2_in_temp)
                dt1_in.append(dt1_in_temp)
                dt2_in.append(dt2_in_temp)

                if phi1_in_temp > 0:
                    phi1_gmt_temp = 90 - phi1_in_temp
                else:
                    phi1_gmt_temp = 90 + (-phi1_in_temp)
                if phi2_in_temp > 0:
                    phi2_gmt_temp = 90 - phi2_in_temp
                else:
                    phi2_gmt_temp = 90 + (-phi2_in_temp)
                phi1_gmt.append(phi1_gmt_temp)
                phi2_gmt.append(phi2_gmt_temp)

                # nulls (occur in steps of 90 deg)
                dt_a = np.squeeze(np.squeeze(models_df_raw["dt_eff"][i_model]))
                ind_dt_max_null = np.argmax(dt_a)
                baz_nulls_theo = [
                    ind_dt_max_null - 270,
                    ind_dt_max_null - 180,
                    ind_dt_max_null - 90,
                    ind_dt_max_null,
                    ind_dt_max_null + 90,
                    ind_dt_max_null + 180,
                    ind_dt_max_null + 270,
                ]
                baz_nulls_temp = []
                for i_null in range(len(baz_nulls_theo)):
                    if baz_nulls_theo[i_null] >= 0 and baz_nulls_theo[i_null] <= 360:
                        baz_nulls_temp.append(baz_nulls_theo[i_null])
# .............................................................................
            case "T1":
                dip_in_temp = int(str(models_df_raw["dip_in"][i_model][0][0]))
                thick_in_temp = int(str(models_df_raw["thick_in"][i_model][0][0]))
                downdipdir_in_temp = int(str(models_df_raw["downdipdir_in"][i_model][0][0]))
                dip_in.append(dip_in_temp)
                thick_in.append(thick_in_temp)
                downdipdir_in.append(downdipdir_in_temp)

                if downdipdir_in_temp <= 90:
                    phi_T1_temp = downdipdir_in_temp
                elif downdipdir_in_temp > 90 and downdipdir_in_temp <= 270:
                    phi_T1_temp = downdipdir_in_temp - 180
                elif downdipdir_in_temp > 270:
                    phi_T1_temp = -(360 - downdipdir_in_temp)
                phi_T1.append(phi_T1_temp)

                # nulls (occur NOT in steps of 90 deg)
                phi_a = np.squeeze(np.squeeze(models_df_raw["phi_eff"][i_model]))
                # in the down-dip direction
                baz_nulls_neg = [downdipdir_in_temp, downdipdir_in_temp + 180]
                baz_nulls_1 = []
                for baz_null in baz_nulls_neg:
                    baz_null_pos = baz_null
                    if baz_null > 360:
                        baz_null_pos = baz_null - 360  # 0 to 360 deg
                    baz_nulls_1.append(baz_null_pos)
                # for down-dip direction orthogonal to backazimuth
                if downdipdir_in_temp < 90:
                    null_1_ortho = phi_a + 90
                    null_1 = min(null_1_ortho)
                    null_2_ortho = phi_a - 90
                    null_2 = max(null_2_ortho)
                elif downdipdir_in_temp > 90 and downdipdir_in_temp < 180:
                    null_1_ortho = phi_a + 90
                    null_1 = max(null_1_ortho)
                    null_2_ortho = phi_a - 90
                    null_2 = min(null_2_ortho)
                elif downdipdir_in_temp > 270:
                    null_1_ortho = phi_a - 90
                    null_1 = max(null_1_ortho)
                    null_2_ortho = phi_a + 90
                    null_2 = min(null_2_ortho)
                elif downdipdir_in_temp > 180 and downdipdir_in_temp < 270:
                    null_1_ortho = phi_a + 90
                    null_1 = max(null_1_ortho)
                    null_2_ortho = phi_a - 90
                    null_2 = min(null_2_ortho)
                else:
                    null_1_ortho = phi_a + 90
                    null_1 = max(null_1_ortho)
                    null_2_ortho = phi_a - 90
                    null_2 = min(null_2_ortho)

                if null_1 < 0:
                    null_1 = 360 + null_1
                if null_2 < 0:
                    null_2 = 360 + null_2
                baz_nulls_2 = [null_1, null_2]

                # Strange special case nulls at 90, 270, 180, 180+delta deg
                if abs(baz_nulls_2[1] - baz_nulls_2[0]) < 1:
                    baz_nulls_2[1] = 0

                baz_nulls_temp = [
                    baz_nulls_1[0], baz_nulls_1[1], baz_nulls_2[0], baz_nulls_2[1]
                ]

# .............................................................................
        baz_nulls_temp = sorted(baz_nulls_temp)
        baz_nulls[i_model] = baz_nulls_temp
        baz_null1.append(baz_nulls_temp[0])
        baz_null2.append(baz_nulls_temp[1])
        baz_null3.append(baz_nulls_temp[2])
        baz_null4.append(baz_nulls_temp[3])

# -----------------------------------------------------------------------------
    match model_type:
        case "H1":
            models_df["phi_in"] = phi_in
            models_df["dt_in"] = dt_in
            models_df["phi_gmt"] = phi_gmt_temp
        case "H2":
            models_df["phi1_in"] = phi1_in
            models_df["phi2_in"] = phi2_in
            models_df["dt1_in"] = dt1_in
            models_df["dt2_in"] = dt2_in
            models_df["phi1_gmt"] = phi1_gmt
            models_df["phi2_gmt"] = phi2_gmt
        case "T1":
            models_df["dip_in"] = dip_in
            models_df["thick_in"] = thick_in
            models_df["downdipdir_in"] = downdipdir_in
            models_df["phi_T1"] = phi_T1

    models_df["baz_nulls"] = baz_nulls
    models_df["baz_null1"] = baz_null1
    models_df["baz_null2"] = baz_null2
    models_df["baz_null3"] = baz_null3
    models_df["baz_null4"] = baz_null4


# %%
# -----------------------------------------------------------------------------
# Select models with model parameters in the selected ranges
# -----------------------------------------------------------------------------
    match model_type:
        case "H1":
            models_df_select = models_df.loc[
                (models_df["phi_in"] >= phi_min)
                & (models_df["phi_in"] <= phi_max)
                & (models_df["dt_in"] >= dt_min)
                & (models_df["dt_in"] <= dt_max)
            ]
        case "H2":
            models_df_select = models_df.loc[
                (models_df["phi1_in"] >= phi1_min)
                & (models_df["phi1_in"] <= phi1_max)
                & (models_df["dt1_in"] >= dt1_min)
                & (models_df["dt1_in"] <= dt1_max)
                & (models_df["phi2_in"] >= phi2_min)
                & (models_df["phi2_in"] <= phi2_max)
                & (models_df["dt2_in"] >= dt2_min)
                & (models_df["dt2_in"] <= dt2_max)
            ]
        case "T1":
            models_df_select = models_df.loc[
                (models_df["dip_in"] >= dip_min)
                & (models_df["dip_in"] <= dip_max)
                & (models_df["thick_in"] >= thick_min)
                & (models_df["thick_in"] <= thick_max)
                & (models_df["downdipdir_in"] >= downdipdir_min)
                & (models_df["downdipdir_in"] <= downdipdir_max)
            ]

# -----------------------------------------------------------------------------
    N_select = len(models_df_select)
    models_df_select["i_select"] = np.arange(N_select).tolist()

    print(f"Data loaded: in total {N_total} models, selected {N_select} models.")
    if N_select == 0:
        print("No models select!")

    return(models_df, models_df_select, model_type, dom_per, N_total, N_select)

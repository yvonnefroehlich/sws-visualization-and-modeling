# ==========================================================================
# Explore forward calculated splitting parameters for structural anisotropy
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
# - Created: 2024/07/09
# - Continued: 2025/01/07
# - Continued: 2025/04/06-08
# - Continued: 2026/01/23 - Use PyGMT v0.18.0 with GMT 6.6.0
# - Continued: 2026/01/24 - Allow setting ranges for model parameters
# - Updated: 2026/01/25 - Improve data preparation, shorten code for plotting nulls
# - Updated: 2026/01/27 - Improve determination of null directions
# - Updated: 2026/01/18 - Move loading anisotropy model to separate function
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
import pygmt
from pygmt.params import Position

from load_models import load_models


# %%
# -----------------------------------------------------------------------------
# Adjust for your needs
# -----------------------------------------------------------------------------
model_type = "H2"  # "H1" | "H2" | "T1"
dom_per = 8  ## 6 | 8 | 10  # in seconds  (TEST data provided for 8 s)
root_path = "C:/Users/Admin/C2/008_github_repos/github_sws"  # Adjust for your file sturcture
path_models = f"{root_path}/sws-visualization-and-modeling/000_test_data"
file_models = "default"
path_out = "02_out_figs"

# -----------------------------------------------------------------------------
# Limits for model parameters
# H1
phi_min = -90
phi_max = 90
dt_min = 0
dt_max = 4
# H2 (index 1 lower layer, index 2 upper layer)
phi1_min = -90
phi1_max = 90
phi2_min = -90
phi2_max = 90
dt1_min = 0
dt1_max = 4
dt2_min = 0
dt2_max = 4
# T1
dip_min = 0
dip_max = 70
thick_min = 0
thick_max = 400
downdipdir_min = 0
downdipdir_max = 360

# -----------------------------------------------------------------------------
status_cb = True  ## True | False
status_per = False  ## True | False
font_map = 9  # in points
font_cb = 14.5  # in points


# %%
# -----------------------------------------------------------------------------
# Load anisotropy models
# -----------------------------------------------------------------------------
models_df, models_df_select, model_type, dom_per, N_total, N_select = load_models(
    model_type=model_type,
    dom_per=dom_per,
    path_models=path_models,
    phi_min=phi_min,
    phi_max=phi_max,
    dt_min=dt_min,
    dt_max=dt_max,
    phi1_min=phi1_min,
    phi1_max=phi1_max,
    phi2_min=phi2_min,
    phi2_max=phi2_max,
    dt1_min=dt1_min,
    dt1_max=dt1_max,
    dt2_min=dt2_min,
    dt2_max=dt2_max,
    dip_min=dip_min,
    dip_max=dip_max,
    thick_min=thick_min,
    thick_max=thick_max,
    downdipdir_min=downdipdir_min,
    downdipdir_max=downdipdir_max,
)

# -----------------------------------------------------------------------------
model_start = 0
model_end = "NaN"  # total number of models
model_step = 1

if model_end == "NaN":
    model_end = N_select

# backazimuth in degrees North to East
baz_step = 1
baz = np.arange(0, 360 + baz_step, baz_step)


# %%
# -----------------------------------------------------------------------------
# General stuff
# -----------------------------------------------------------------------------
box_standard = "+glightgray@30+p0.1p,gray30+r1p"

# Colors based on Fröhlich et al. (2024) GJI
color_highlight = "255/90/0"  # -> orange
color_H1 = "127/140/95"  # 1 horizontal layer (H1) -> green
color_H2l = "178/34/34"  # 2 horizontal layers (H2) lower (first) layer -> red
color_H2u = "24/116/205"  # 2 horizontal layers (H2) upper (second) layer -> blue
color_T1 = "218/163/109"  # 1 tilted layer (T1) -> brown

args_nulls = {"style": "c0.15c", "fill": "white", "pen": "0.7p"}
args_nulls_cath = {"style": "c0.07c", "fill": "white", "pen": "0.5p"}
phi_ys = np.arange(-90, 90 + 10, 10)
dt_ys = np.arange(0, 4 + 0.1, 0.2)
baz_null_add = 5


# %%
# -----------------------------------------------------------------------------
# Make plots of anisotropy models
# -----------------------------------------------------------------------------
for i_model in range(model_start, model_end, model_step):

    model_out = models_df_select[models_df_select["i_select"] == i_model]
    i_total = int(model_out["i_total"].iloc[0])

    # apparent splitting parameters
    phi_a = np.squeeze(np.squeeze(model_out["phi_eff"]))
    dt_a = np.squeeze(np.squeeze(model_out["dt_eff"]))

    # nulls
    baz_nulls = model_out["baz_nulls"][i_total]

    # model parameters
    match model_type:
        case "H1":
            phi = model_out["phi_in"][i_total]
            dt = model_out["dt_in"][i_total]
            phi_gmt = model_out["phi_gmt"][i_total]
        case "H2":
            phi_1 = model_out["phi1_in"][i_total]
            phi_2 = model_out["phi2_in"][i_total]
            dt_1 = model_out["dt1_in"][i_total]
            dt_2 = model_out["dt2_in"][i_total]
            phi_1_gmt = model_out["phi1_gmt"][i_total]
            phi_2_gmt = model_out["phi2_gmt"][i_total]
        case "T1":
            downdipdir = model_out["downdipdir_in"][i_total]
            dip = model_out["dip_in"][i_total]
            thick = model_out["thick_in"][i_total]
            downdipdir_gmt = 90 - downdipdir
            strike_gmt = downdipdir_gmt + 90
            phi_T1 = model_out["phi_T1"][i_total]

# -----------------------------------------------------------------------------
    fig = pygmt.Figure()
    pygmt.config(MAP_GRID_PEN_PRIMARY="0.01p,gray80", FONT=f"{font_map}p")

    pygmt.makecpt(cmap="phase", series=[-90, 90], cyclic=True)

# .............................................................................
# Left: Cartesian plots of splitting parameters
# .............................................................................
    x_hline = [-10, 360]
    proj_stereo = "X10c/4c"
    args_mp_line = {"x": x_hline, "no_clip": True}
    args_nulls_fill = {"y": [-90, -90, 90, 90, -90], "fill": "gray80@50"}
    args_nulls_line = {"y": [-90, 90], "pen": "1p,gray30,2_4"}

    # Top Left: fast polarization direction
    label_phi = "app. fast pol. dir. @~f@~@-a@- / N°E"
    if model_type == "H1":
        label_phi = "fast pol. dir. @~f@~@-a@- / N°E"
    fig.basemap(
        region=[0, 360, -90, 90],
        projection=proj_stereo,
        frame=["WSne", "xf10g30", f"ya30f10g30+l{label_phi}"],
    )

    match model_type:
        case "H1":
            fig.plot(y=[phi] * 2, pen=f"1p,{color_H1},dashed", **args_mp_line)
        case "H2":
            fig.plot(y=[phi_1] * 2, pen=f"1p,{color_H2l},dashed", **args_mp_line)
            fig.plot(y=[phi_2] * 2, pen=f"1p,{color_H2u},dashed", **args_mp_line)
        case "T1":
            fig.plot(y=[phi_T1] * 2, pen=f"1p,{color_T1},dashed", **args_mp_line)

    fig.plot(x=baz, y=phi_a, pen="0.1p")
    fig.plot(x=baz, y=phi_a, style="c0.07c", fill=phi_a, cmap=True)

    for i_null in range(len(baz_nulls)):
        fig.plot(
            x=[
                baz_nulls[i_null] - baz_null_add,
                baz_nulls[i_null] + baz_null_add,
                baz_nulls[i_null] + baz_null_add,
                baz_nulls[i_null] - baz_null_add,
                baz_nulls[i_null] - baz_null_add,
            ],
            **args_nulls_fill,
        )
        fig.plot(x=[baz_nulls[i_null], baz_nulls[i_null]], **args_nulls_line)

    fig.shift_origin(yshift="-h-0.5c")

    # Bottom Left: delay time
    label_dt = "app. delay time @~d@~t@-a@- / s"
    if model_type == "H1":
        label_dt = "delay time @~d@~t@-a@- / s"
    fig.basemap(
        region=[0, 360, 0, 4],
        projection=proj_stereo,
        frame=["WSne", "xa30f10g30+lbackazimuth / °", f"ya1f0.25g0.5+l{label_dt}"],
    )

    match model_type:
        case "H1":
            fig.plot(y=[dt] * 2, pen=f"1p,{color_H1},dashed", **args_mp_line)
        case "H2":
            fig.plot(y=[dt_1] * 2, pen=f"1p,{color_H2l},dashed", **args_mp_line)
            fig.plot(y=[dt_2] * 2, pen=f"1p,{color_H2u},dashed", **args_mp_line)

    fig.plot(x=baz, y=dt_a, pen="0.1p")
    fig.plot(x=baz, y=dt_a, style="c0.07c", fill=phi_a, cmap=True)

    for i_null in range(len(baz_nulls)):
        fig.plot(
            x=[
                baz_nulls[i_null] - baz_null_add,
                baz_nulls[i_null] + baz_null_add,
                baz_nulls[i_null] + baz_null_add,
                baz_nulls[i_null] - baz_null_add,
                baz_nulls[i_null] - baz_null_add,
            ],
            **args_nulls_fill,
        )
        fig.plot(x=[baz_nulls[i_null], baz_nulls[i_null]], **args_nulls_line)

    fig.shift_origin(xshift="+w+1.5c", yshift="4.5c")

# .............................................................................
# Top Right: model parameter of anisotropy model
# .............................................................................
    size = 2
    region_mp = [-size, size] * 2
    fig.basemap(region=region_mp, projection="X4c/4c", frame="g0.5")

    if model_type in ["H1", "H2"]:
        # Plot top view sign
        fig.plot(x=-1.5, y=-1.5, style="x0.3c", pen="1.2p")
        fig.plot(x=-1.5, y=-1.5, style="c0.3c", pen="1.2p")
        # Mark central lines
        fig.plot(x=[-size, size], y=[0] * 2, pen="0.5p")
        fig.plot(x=[0] * 2, y=[-size, size], pen="0.5p")

    args_leg_bar = {"x": 4, "y": 4, "style": "j0/0.5/0.1"}
    match model_type:
        case "H1":
            bar_H1 = f"j{phi_gmt}/{dt}/0.1"
            label_H1 = f"{phi} N°E | {dt} s"

            fig.plot(x=0, y=0, style=bar_H1, fill=color_H1)

            fig.plot(fill=color_H1, label=label_H1, **args_leg_bar)
        case "H2":
            bar_H2l = f"j{phi_1_gmt}/{dt_1}/0.1"
            bar_H2u = f"j{phi_2_gmt}/{dt_2}/0.1"
            label_H2l = f"lower: {phi_1} N°E | {dt_1} s"
            label_H2u = f"upper: {phi_2} N°E | {dt_2} s"

            # Note order in main plot: stacking approach -> lower then upper layer
            fig.plot(x=0, y=0, style=bar_H2l, fill=color_H2l, pen="0.5p")
            fig.plot(x=0, y=0, style=bar_H2u, fill=color_H2u, pen="0.5p,white")

            # Note order in legend: from top to bottom -> upper then lower layer
            fig.plot(
                fill=color_H2u, label=label_H2u, pen="0.1p,white", **args_leg_bar
            )
            fig.plot(fill=color_H2l, label=label_H2l, pen="0.1p", **args_leg_bar)
        case "T1":
            bar_T1 = f"j{strike_gmt}/1/0.1"
            vec_T1_ddd = ([downdipdir_gmt], [1])  # Input for vector must be a list
            vec_T1_dip = ([-dip], [1.5])
            label_T1 = f"{dip} ° | {downdipdir} N°E | {thick} km"

            fig.plot(
                x=-0.25,
                y=0,
                style=f"w1.5/{-float(dip)}/0",
                fill="lightgray",
                pen="1p,gray30",
            )
            fig.plot(
                x=-0.25,
                y=0,
                style="v0c",
                direction=vec_T1_dip,
                fill=color_T1,
                pen=f"4p,{color_T1}",
            )
            fig.plot(x=[-0.4, -0.1], y=[0.025] * 2, pen="2p,white")
            fig.plot(x=0, y=0.23, style="i0.4c", fill="gold", pen="0.2p,black")
            fig.plot(x=[-1, 1], y=[0] * 2, pen="2.25p,black")

            fig.plot(fill=color_T1, label=label_T1, pen="0.1p", **args_leg_bar)

    with pygmt.config(FONT="8p"):
        fig.legend(position="jTC+w3.8c+o0c/0.1c", box=box_standard)

    # fig.shift_origin(yshift="-h-1c")
    fig.shift_origin(xshift="-0.5c", yshift="-h-0.7c")

# .............................................................................
# Bottom Right: Stereoplot for splitting parameter
# .............................................................................
    rho = 0.6
    step = 3
    baz_stereo = baz[0::step]
    phi_a_stereo = phi_a[0::step]
    dt_a_stereo = dt_a[0::step]

    with pygmt.config(FORMAT_GEO_MAP="+D"):  # 0 to 360 deg
        fig.basemap(region=[0, 360, 0, 1], projection="P4c+a", frame="xa30f10g30")

    for i_bar in range(len(baz_stereo)):
        if phi_a_stereo[i_bar] > 0:
            phi_a_stereo_gmt = 90 - phi_a_stereo[i_bar]
        else:
            phi_a_stereo_gmt = 90 + (-phi_a_stereo[i_bar])
        data_bar = [
            [
                baz_stereo[i_bar],  # angle
                rho,  # rho
                phi_a_stereo[i_bar],  # color
                phi_a_stereo_gmt,  # orientation
                dt_a_stereo[i_bar] / 2,  # length
                0.02,  # thickness
            ]
        ]
        fig.plot(data=data_bar, style="j", cmap=True, no_clip=True)

    # Mark theoretical null directions
    fig.plot(x=baz_nulls, y=[rho] * 4, no_clip=True, **args_nulls)

    # For T1 add arrow
    if model_type == "T1":
        fig.plot(x=baz_nulls, y=[rho] * 4, no_clip=True, **args_nulls)
        fig.plot(
            x=0,
            y=0,
            style="v0.3c+e+h0+a45",
            direction=vec_T1_ddd,
            fill=color_T1,
            pen=f"3p,{color_T1}",
        )
        fig.plot(x=0, y=0, style=bar_T1, fill="black")

# -----------------------------------------------------------------------------
    # Add colorbar for fast polarization direction
    if status_cb:
        with pygmt.config(
            FONT=font_cb, MAP_TICK_LENGTH_PRIMARY="4p", MAP_FRAME_PEN="0.5p"
        ):
            fig.colorbar(
                position=Position("CT", cstype="inside", offset=("-2.33c", "-4.69c")),
                length="4c",
                width="0.2c",
                orientation="vertical",
                move_text="label",
                # white space as y-label to move cyclic arrow symbol down
                frame=["xa30f10+lapp. fast pol. dir. @~f@~@-a@- / N@.E", "y+l "],
            )

    # Add label for dominant period
    if status_per:
        fig.text(
            text=f"{dom_per} s",
            position="TL",
            justify="MC",
            offset="-0.35c/0c",
            font=f"7.5p,{color_highlight}",
            fill="white@30",
            pen="0.01p,black",
            clearance="0.08c/0.08c+tO",
            no_clip=True,
        )

# -----------------------------------------------------------------------------
    fig.show()
    fig_name_basic = f"forwardt_syn_sp_period{dom_per}s_{model_type}"

    str_cb = ""
    if not status_cb:
        str_cb = "NO"
    str_per = ""
    if not status_per:
        str_per = "NO"

    fig_name_mt = ""
    match model_type:
        case "H1":
            fig_name_mt = f"phi{phi}deg_dt{dt}s"
        case "H2":
            fig_name_mt = f"phil{phi_1}deg_phiu{phi_2}deg_dtl{dt_1}s_dtu{dt_2}s"
        case "T1":
            fig_name_mt = f"thick{thick}km_dip{dip}deg_ddd{downdipdir}deg"

    for ext in ["png", "pdf", "eps"]:
        fig_name = f"{fig_name_basic}_{fig_name_mt}_cb{str_cb}_per{str_per}"
        if ext == "png":
            fig_name = f"{i_total}_{fig_name}"
        # fig.savefig(fname=f"{path_out}/{model_type}/{fig_name}.{ext}", dpi=720)

    print(f"{i_total}_{fig_name}")

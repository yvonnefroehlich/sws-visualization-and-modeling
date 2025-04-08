# #############################################################################
# Explore forwardt calculated splitting parameters for structural anisotropy
# -----------------------------------------------------------------------------
# Supported model types
# - One horizontal layer (with HTI): H1
# - Two horizontal layer (with HTI): H2
# - One tilted layer (with TTI): T1
# Previously calculated synthetic splitting parameters
# - https://github.com/yvonnefroehlich/sws-visualization-and-modeling/tree/main/003_modeling
# - The output MATLAB struct is split into separate structs for the different model types
# -----------------------------------------------------------------------------
# History
# - Created: 2024/07/09
# - Continued: 2025/01/07
# - Continued: 2025/04/06-08
# -----------------------------------------------------------------------------
# Versions
#   PyGMT v0.14.0 -> https://www.pygmt.org/v0.14.0/ | https://www.pygmt.org/
#   GMT 6.5.0 -> https://www.generic-mapping-tools.org/
# -----------------------------------------------------------------------------
# Contact
# - Author: Yvonne Fröhlich
# - ORCID: https://orcid.org/0000-0002-8566-0619
# - GitHub: https://github.com/yvonnefroehlich/sws-visualization-and-modeling
# -----------------------------------------------------------------------------
# Related to
# - Fröhlich Y. & Ritter J. R. R. (2024)
#   Vertical and Small-scale Lateral Varying Seismic Anisotropy in the Upper Mantle
#   Underneath the Upper Rhine Graben, Central Europe.
#   Annual Meeting of the American Geophysical Union,
#   Division Session Exploring Innovations and New Directions in Seismic Anisotropy and
#   Attenuation: Observations, Models, and Experiments I Oral, DI21A-02.
#   Abstract ID 1578275, https://agu.confex.com/agu/agu24/meetingapp.cgi/Paper/1578275.
#   http://dx.doi.org/10.5281/zenodo.14510993.
# - Fröhlich Y., Grund M. & Ritter J. R. R. (2024)
#   Lateral and vertical variations of seismic anisotropy in the lithosphere-asthenosphere
#   system underneath Central Europe from long-term splitting measurements.
#   Geophysical Journal International, 239(1), 112-135.
#   https://doi.org/10.1093/gji/ggae245.
# #############################################################################


import numpy as np
import pygmt
from scipy import io


# %%
# -----------------------------------------------------------------------------
# General stuff
# -----------------------------------------------------------------------------
path_in = "01_in_data"
path_out = "02_out_figs"

status_cb = True  ## True, False
status_per = True  ## True, False
font_size = 8  # in points

dom_per = 8  ## 6, 8, 10  # in seconds
model_type = "H1"  ## H1, H2, T1
print(f"Dominant period {dom_per} s - Model type {model_type}")

models = f"sws_modout_domper{dom_per}s_{model_type}.mat"
models_mat = io.loadmat(f"{path_in}/{models}")
N_models = len(models_mat["model_out"][0])
print(f"Data loaded - {N_models} models!\nStarting with making plots!")

# -----------------------------------------------------------------------------
model_start = 0
model_end = N_models
model_step = 1

baz_step = 1
baz = np.arange(0, 360 + baz_step, baz_step)  # backazimuth in degrees North to East

# -----------------------------------------------------------------------------
box_standard = "+glightgray@30+p0.1p,gray30+r1p"

# Colors based on Fröhlich et al. (2024) GJI
color_highlight = "255/90/0"  # -> orange
color_H1 = "127/140/95"  # 1 horizontal layer (H1) -> green
color_H2lower = "178/34/34"  # 2 horizontal layers (H2) lower (first) layer -> red
color_H2upper = "24/116/205"  # 2 horizontal layers (H2) upper (second) layer -> blue
color_T1 = "218/163/109"  # 1 tilted layer (T1) -> brown


# %%
# -----------------------------------------------------------------------------
# Make plots of anisotropy models
# -----------------------------------------------------------------------------

for i_model in range(model_start, model_end + model_step, model_step):

    model_out = models_mat["model_out"][0][i_model]

    phi_a = np.squeeze(model_out[0])
    dt_a = np.squeeze(model_out[1])

    match model_type:
        case "H1":
            phi = str(model_out[2][0])[1:-1]
            dt = str(model_out[3][0])[1:-1]
            if float(phi) > 0:
                phi_gmt = 90 - float(phi)
            else:
                phi_gmt = 90 + (-float(phi))
        case "H2":
            phi_1 = str(model_out[2][0])[1:-1]
            phi_2 = str(model_out[2][1])[1:-1]
            dt_1 = str(model_out[3][0])[1:-1]
            dt_2 = str(model_out[3][1])[1:-1]
            if float(phi_1) > 0:
                phi_1_gmt = 90 - float(phi_1)
            else:
                phi_1_gmt = 90 + (-float(phi_1))
            if float(phi_2) > 0:
                phi_2_gmt = 90 - float(phi_2)
            else:
                phi_2_gmt = 90 + (-float(phi_2))
        case "T1":
            downdipdir = str(model_out[2][0])[1:-1]
            dip = str(model_out[3][0])[1:-1]
            thick = str(model_out[4][0])[1:-1]
            downdipdir_gmt = 90 - float(downdipdir)
            strike_gmt = downdipdir_gmt + 90
            if float(downdipdir) <= 90:
                phi = float(downdipdir)
            elif float(downdipdir) > 90 and float(downdipdir) <= 270:
                phi = float(downdipdir) - 180
            elif float(downdipdir) > 270:
                phi = -(360 - float(downdipdir))

# -----------------------------------------------------------------------------
    fig = pygmt.Figure()
    pygmt.config(MAP_GRID_PEN_PRIMARY="0.01p,gray80", FONT=f"{font_size}p")

    pygmt.makecpt(cmap="phase", series=[-90, 90], cyclic=True)

# .............................................................................
    # Left: Cartesian plots of splitting parameters
# .............................................................................
    x_hline = [-10, 360]
    proj_stereo = "X10c/4c"

    # Top Left: fast polarization direction
    fig.basemap(
        region=[0, 360, -90, 90],
        projection=proj_stereo,
        frame=["WSne", "xf10g30", "ya30f10g30+lapp. fast pol. dir. @~f@~@-a@- / N°E"],
    )
    match model_type:
        case "H1":
            fig.plot(x=x_hline, y=[phi] * 2, pen=f"1p,{color_H1},dashed", no_clip=True)
        case "H2":
            fig.plot(x=x_hline, y=[phi_1] * 2, pen=f"1p,{color_H2lower},dashed", no_clip=True)
            fig.plot(x=x_hline, y=[phi_2] * 2, pen=f"1p,{color_H2upper},dashed", no_clip=True)
        case "T1":
            fig.plot(x=x_hline, y=[phi] * 2, pen=f"1p,{color_T1},dashed", no_clip=True)
    fig.plot(x=baz, y=phi_a, pen="0.1p")
    fig.plot(x=baz, y=phi_a, style="c0.07c", fill=phi_a, cmap=True)

    fig.shift_origin(yshift="-h-0.5c")

    # Bottom Left: delay time
    fig.basemap(
        region=[0, 360, 0, 4],
        projection=proj_stereo,
        frame=[
            "WSne", "xa30f10g30+lbackazimuth / °", "ya1f0.25g0.5+lapp. delay time @~d@~t@-a@- / s",
        ],
    )
    match model_type:
        case "H1":
            fig.plot(x=x_hline, y=[dt] * 2, pen=f"1p,{color_H1},dashed", no_clip=True)
        case "H2":
            fig.plot(x=x_hline, y=[dt_1] * 2, pen=f"1p,{color_H2lower},dashed", no_clip=True)
            fig.plot(x=x_hline, y=[dt_2] * 2, pen=f"1p,{color_H2upper},dashed", no_clip=True)
    fig.plot(x=baz, y=dt_a, pen="0.1p")
    fig.plot(x=baz, y=dt_a, style="c0.07c", fill=phi_a, cmap=True)

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
            bar_H2lower = f"j{phi_1_gmt}/{dt_1}/0.1"
            bar_H2upper = f"j{phi_2_gmt}/{dt_2}/0.1"
            label_H2lower = f"lower: {phi_1} N°E | {dt_1} s"
            label_H2upper = f"upper: {phi_2} N°E | {dt_2} s"

            # Note order in main plot: stacking approach -> lower then upper layer
            fig.plot(x=0, y=0, style=bar_H2lower, fill=color_H2lower, pen="0.5p")
            fig.plot(x=0, y=0, style=bar_H2upper, fill=color_H2upper, pen="0.5p,white")

            # Not order in legend: from top to bottom -> upper then lower layer
            fig.plot(fill=color_H2upper, label=label_H2upper, pen="0.1p,white", **args_leg_bar)
            fig.plot(fill=color_H2lower, label=label_H2lower, pen="0.1p", **args_leg_bar)
        case "T1":
            bar_T1 = f"j{strike_gmt}/1/0.1"
            vec_T1_ddd = ([downdipdir_gmt], [1])  # Input for vector must be a list
            vec_T1_dip = ([-float(dip)], [1.5])
            label_T1 = f"{downdipdir} N°E | {dip} ° | {thick} km"

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

    with pygmt.config(FORMAT_GEO_MAP="+D"):  # 0°-360°
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

    match model_type:
        # Mark theoretical null directions
        case "H1":
            phi = int(phi)
            baz_nulls = np.array([phi, phi + 90, phi + 180, phi + 270])
            fig.plot(x=baz_nulls, y=[rho] * 4, style="c0.15c", fill="white", pen="0.7p")
        # TODO case "H2":
        # Arrow showing strike / down dip direction
        case "T1":
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
    if status_cb == True:
        with pygmt.config(FONT="13p", MAP_TICK_LENGTH_PRIMARY="4p", MAP_FRAME_PEN="0.5p"):
            fig.colorbar(
                position="jCT+w4c/0.2c+o-2.33c/-4.69c+v+ml",
                # white space as y label to move cyclic arrow symbole down
                frame=["xa30f10+lapp. fast pol. dir. @~f@~@-a@- / N@.E", "y+l "],
        )

    # Add label for dominant period
    if status_per == True:
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
    if status_cb == False:
        str_cb = "NO"
    str_per = ""
    if status_per == False:
        str_per = "NO"

    fig_name_mt =  ""
    match model_type:
        case "H1":
            fig_name_mt = f"phi{phi}deg_dt{dt}s"
        case "H2":
            fig_name_mt = f"phil{phi_1}deg_phiu{phi_2}deg_dtl{dt_1}s_dtu{dt_2}s"
            fig_name_mt = f"phi{phi_1}deg_phi{phi_2}deg_dt{dt_1}s_dt{dt_2}s"
        case "T1":
            fig_name_mt = f"thick{thick}km_dip{dip}deg_ddd{downdipdir}deg"

    fig_name = f"{fig_name_basic}_{fig_name_mt}_cb{str_cb}_per{str_per}"
    for ext in ["png"]: #, "pdf", "eps"]:
        fig.savefig(fname=f"{path_out}/{model_type}/{fig_name}.{ext}", dpi=720)
    print(f"{i_model} | {fig_name}")

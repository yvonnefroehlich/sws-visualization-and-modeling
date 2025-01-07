# #############################################################################
# Forwardt calculated splitting parameters
# -----------------------------------------------------------------------------
# Author: Yvonne Fröhlich
# ORCID: https://orcid.org/0000-0002-8566-0619
# GitHub: https://github.com/yvonnefroehlich/gmt-pygmt-plotting
# -----------------------------------------------------------------------------
# - Created: 2024/07/09
# - Continued: 2025/01/07
#   PyGMT v0.14.0 -> https://www.pygmt.org/v0.14.0/ | https://www.pygmt.org/
#   GMT 6.5.0 -> https://www.generic-mapping-tools.org/
# #############################################################################


import numpy as np
import pygmt
from scipy import io

# %%
# -----------------------------------------------------------------------------
# General stuff
# ----------------------------------------------------------------------------
path_in = "01_in_data"
path_out = "02_out_figs"

dom_per = 8  ## 6, 8, 10  # in seconds
model_type = "H2"  ## H1, H2, T1
print(f"Period {dom_per} s - Model Type {model_type}")

models = f"sws_modout_domper{dom_per}s_{model_type}.mat"
models_mat = io.loadmat(f"{path_in}/{models}")
print("Data loaded! Starting with making plots!")

box_standard = "+gwhite@30+p0.1p,gray30+r2p"
color_highlight = "255/90/0"  # ->  orange | URG paper

color_H1 = "127/140/95"  # -> 1 horizontal layer; green
color_lower = "178/34/34"  # -> 2 horizontal layers lower / first layer; red
color_upper = "24/116/205"  # -> 2 horizontal layers upper / second layer; blue
color_T1 = "218/163/109"  # -> 1 tiled layer brown


# %%
# -----------------------------------------------------------------------------
# Make plot of splitting parameters
# -----------------------------------------------------------------------------
baz = np.arange(0, 361, 1)

# -----------------------------------------------------------------------------
model_start = 1 # 6760
model_step = 1
# - H1
#   * phi: [-80:10:90] deg | 17 |
#   * dt: [0.25:0.25:4] s | 15 |
# - H2
#   * phi1: [-80:10:90] deg | 17 |
#   * phi2: [-80:10:90] deg | 17 |
#   * dt1: [0.25:0.25:4] s | 15 |
#   * dt2: [0.25:0.25:4] s | 15 |
# - T1
#   * thickness: [10:10:250] km | 25 | 252
#   * dip angle: [10:10:70] deg | 7 |
#   * down-dip direction: [0:10:350] deg | 35 |

for i_model in range(
    model_start, len(models_mat["model_out"][0]) + model_step, model_step
):
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
    pygmt.config(MAP_GRID_PEN_PRIMARY="0.05p,gray30")

    pygmt.makecpt(cmap="phase", series=[-90, 90], cyclic=True)

    # -----------------------------------------------------------------------------
    # Top Left: fast polarization direction
    fig.basemap(
        region=[0, 360, -90, 90],
        projection="X10c/4c",
        frame=["WSne", "xf10g30", "ya30f10g30+lapp. fast pol. dir. / N°E"],
    )
    match model_type:
        case "H1":
            fig.plot(x=[-10, 370], y=[phi, phi], pen=f"1p,{color_H1},dashed", no_clip=True)
        case "H2":
            fig.plot(x=[-10, 370], y=[phi_1, phi_1], pen=f"1p,{color_lower},dashed", no_clip=True)
            fig.plot(x=[-10, 370], y=[phi_2, phi_2], pen=f"1p,{color_upper},dashed", no_clip=True)
        case "T1":
            fig.plot(x=[-10, 370], y=[phi, phi], pen=f"1p,{color_T1},dashed", no_clip=True)
    fig.plot(x=baz, y=phi_a, pen="0.1p")
    fig.plot(x=baz, y=phi_a, style="c0.07c", fill=phi_a, cmap=True)

    fig.shift_origin(yshift="-h-0.5c")

    # -----------------------------------------------------------------------------
    # Bottom Left: delay time
    fig.basemap(
        region=[0, 360, 0, 4],
        projection="X10c/4c",
        frame=[
            "WSne",
            "xa30f10g30+lbackazimuth / °",
            "ya1f0.2g0.5+lapp. delay time / s",
        ],
    )
    if model_type == "H1":
        fig.plot(x=[-10, 370], y=[dt, dt], pen=f"1p,{color_H1},dashed", no_clip=True)
    if model_type == "H2":
        fig.plot(x=[-10, 370], y=[dt_1, dt_1], pen=f"1p,{color_lower},dashed", no_clip=True)
        fig.plot(x=[-10, 370], y=[dt_2, dt_2], pen=f"1p,{color_upper},dashed", no_clip=True)
    fig.plot(x=baz, y=dt_a, pen="0.1p")
    fig.plot(x=baz, y=dt_a, style="c0.07c", fill=phi_a, cmap=True)

    fig.shift_origin(xshift="+w+1c", yshift="4.5c")

    # -----------------------------------------------------------------------------
    # Top Right: model parameter
    size = 2
    region_mp = [-size, size, -size, size]
    fig.basemap(region=region_mp, projection="X4c/4c", frame="g0.5")
    # frame=["WSne", "a1f0.1g0.5+u s"],

    if model_type in ["H1", "H2"]:
        fig.plot(x=-1.5, y=-1.5, style="x0.3c", pen="1.2p")
        fig.plot(x=-1.5, y=-1.5, style="c0.3c", pen="1.2p")
        fig.plot(x=[-size, size], y=[0, 0], pen="0.5p")
        fig.plot(x=[0, 0], y=[-size, size], pen="0.5p")
    match model_type:
        case "H1":
            bar_H1 = f"j{phi_gmt}/{dt}/0.1"
            label_H1 = f"{phi} N°E | {dt} s"

            fig.plot(x=0, y=0, style=bar_H1, fill=color_H1, pen="0.5p")

            fig.plot(x=4, y=4, style="j0/0.5/0.05", fill=color_H1, label=label_H1)
            fig.legend(position="jTC+w3.9c+o0c/0.05c", box=box_standard)
        case "H2":
            bar_lower = f"j{phi_1_gmt}/{dt_1}/0.1"
            bar_upper = f"j{phi_2_gmt}/{dt_2}/0.1"
            label_lower = f"lower: {phi_1} N°E | {dt_1} s"
            label_upper = f"upper: {phi_2} N°E | {dt_2} s"

            fig.plot(x=0, y=0, style=bar_lower, fill=color_lower, pen="0.5p")
            fig.plot(x=0, y=0, style=bar_upper, fill=color_upper, pen="0.5p,white")

            fig.plot(x=4, y=4, style="j0/0.5/0.05", fill=color_upper, label=label_upper)
            fig.plot(x=4, y=4, style="j0/0.5/0.05", fill=color_lower, label=label_lower)
            fig.legend(position="jTC+w3.9c+o0c/0.05c", box=box_standard)
        case "T1":
            bar_T1 = f"j{strike_gmt}/1/0.1"
            vec_T1_ddd = ([downdipdir_gmt], [1])
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
            fig.plot(x=[-0.4, -0.1], y=[0.025, 0.025], pen="2p,white")
            fig.plot(x=0, y=0.23, style="i0.4c", fill="gold", pen="0.2p,black")
            fig.plot(x=[-1, 1], y=[0, 0], pen="2.25p,black")

            fig.plot(x=4, y=4, style="j0/0.5/0.05", fill=color_T1, label=label_T1)
            fig.legend(position="jTC+w3.9c+o0c/0.05c", box=box_standard)

    fig.shift_origin(yshift="-h-1c")

    # -----------------------------------------------------------------------------
    # Bottom Right: stereoplot
    rho = 0.6
    step = 3
    baz_stereo = baz[0::step]
    phi_a_stereo = phi_a[0::step]
    dt_a_stereo = dt_a[0::step]

    pygmt.config(FORMAT_GEO_MAP="+D")  # 0°-360°
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
        case "H1":
            phi = int(phi)
            baz_nulls = np.array([phi, phi + 90, phi + 180, phi + 270])
            fig.plot(
                x=baz_nulls, y=np.ones(4) * rho, style="c0.1c", fill="white", pen="0.7p"
            )
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

    fig.text(
        text=f"{dom_per} s",
        position="TL",
        justify="MC",
        offset="-0.45c/0.73c",
        font=f"9p,{color_highlight}",
        fill="white@30",
        pen="0.01p,black",
        clearance="0.08c/0.08c+tO",
        no_clip=True,
    )

    # -----------------------------------------------------------------------------
    fig.show()
    fig_name = f"forwardt_syn_sp_period{dom_per}s_{model_type}"
    fig_name_add = ""
    match model_type:
        case "H1":
            fig_name_add = f"_phi{phi}deg_dt{dt}s"
        case "H2":
            fig_name_add = f"_{phi_1}deg_{phi_2}deg_{dt_1}s_{dt_2}s"
        case "T1":
            fig_name_add = f"_d{thick}km_dip{dip}deg_ddd{downdipdir}deg"
    for ext in ["png"]: #, "pdf", "eps"]:
        fig.savefig(fname=f"{path_out}/{model_type}/{fig_name}{fig_name_add}.{ext}")
    print(f"{fig_name}{fig_name_add}")

# #############################################################################
# Load CSV files with shear wave splitting measurements (SWSMs)
#
# Related to the publications
# (1) Ritter J. R. R., Fröhlich Y., Sanz Alonso Y. & Grund M. (2022).
#     Short-scale laterally varying SK(K)S shear wave splitting at BFO,
#     Germany – implications for the determination of anisotropic structures.
#     Journal of Seismology, 26, 1137-1156.
#     https://doi.org/10.1007/s10950-022-10112-w.
# (2) Fröhlich Y., Grund M., Ritter J. R. R. (2024).
#     Lateral and vertical variations of seismic anisotropy in the
#     lithosphere-asthenosphere system underneath Central Europe from
#     long-term splitting measurements. Geophysical Journal International
#     accepted June 27 2024.
#
# Data is available from RADAR4KIT
# (1) https://dx.doi.org/10.35097/684
# (2) https://dx.doi.org/10.35097/685
# -----------------------------------------------------------------------------
# created: 2024/07/08
# author: Yvonne Fröhlich
# contact: yvonne.froehlich@kit.edu
# #############################################################################


# %%
# -----------------------------------------------------------------------------
# Funktion to load SWSMs into pandas.dataframe
# -----------------------------------------------------------------------------

def read_radar4kit_data(station, network, obstyp):
    """
    Parameters
    ----------
    station : str
        Station code of recording station.
    network : str
        Network code of seismological network.
    obstyp : str
        Observation typ, either "NULLS" or "SPLITS".

    Returns
    -------
    df_swsm : pandas.dataframe
        Dataframe with shear wave splitting measurements at one single
        recording station; either nulls or splits.
    """

    import pandas as pd

    df_swsm = pd.read_csv(
        f"splitresults_{obstyp}_goodfair_{network}_{station}.csv",
        sep=";",
        header=15,
        dtype={'year_jday': 'str'},  # Do not load as floating point number
    )

    df_swsm[["year", "jday"]] = df_swsm["year_jday"].str.split(".", expand=True)

    return df_swsm


# %%
# -----------------------------------------------------------------------------
# Example
# -----------------------------------------------------------------------------

station = "BFO"
network = "GR"
obstyp = "NULLS"

df_swsm = read_radar4kit_data(station=station, network=network, obstyp=obstyp)
# df_swsm.head()
# df_swsm.columns

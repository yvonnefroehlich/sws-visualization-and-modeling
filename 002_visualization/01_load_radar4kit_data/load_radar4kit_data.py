# #############################################################################
# Load CSV files with shear wave splitting measurements (SWSMs)
# -----------------------------------------------------------------------------
# Datasets are available from RADAR4KIT
# - Upper Rhine Graben area: https://dx.doi.org/10.35097/685
#   related to Fröhlich et al. (2024)
# - Blackforest Observatory: https://dx.doi.org/10.35097/684
#   related to Ritter et al. (2022)
# -----------------------------------------------------------------------------
# Related to the publications
# - Fröhlich Y., Grund M. & Ritter J. R. R. (2024).
#   Lateral and vertical variations of seismic anisotropy in the
#   lithosphere-asthenosphere system underneath Central Europe from
#   long-term splitting measurements.
#   Geophysical Journal International, 239(1), 112-135.
#   https://doi.org/10.1093/gji/ggae245.
# - Ritter J. R. R., Fröhlich Y., Sanz Alonso Y. & Grund M. (2022).
#   Short-scale laterally varying SK(K)S shear wave splitting at BFO,
#   Germany – implications for the determination of anisotropic structures.
#   Journal of Seismology, 26, 1137-1156.
#   https://doi.org/10.1007/s10950-022-10112-w.
# -----------------------------------------------------------------------------
# Created: 2024/07/08
# Author: Yvonne Fröhlich
# Contact: https://orcid.org/0000-0002-8566-0619
# #############################################################################


import pandas as pd

# Example
root_path = ""
header = 15  # Number of header lines, given in README
station = "BFO"
network = "GR"
obstyp = "NULLS"

# Read CSV file into pandas.dataframe
df_swsm = pd.read_csv(
    f"{root_path}splitresults_{obstyp}_goodfair_{network}_{station}.csv",
    sep=";",
    header=header,
    # Do not load column "year_jday" as floating point number
    dtype={'year_jday': 'str'},
)

# Split column "year_jday" into two new columns "year" and "jday"
df_swsm[["year", "jday"]] = df_swsm["year_jday"].str.split(".", expand=True)

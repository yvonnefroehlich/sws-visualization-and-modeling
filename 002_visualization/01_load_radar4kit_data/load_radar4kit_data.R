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


# Example 
station = "BFO"
network = "GR"
obstyp = "NULLS"

# Load CSV file into dataframe
swsm_data <- read.csv(
  file=paste("splitresults_",obstyp,"_goodfair_",network,"_",station,".csv", sep=""),
  sep=";",
  skip=15,
)

# Convert from Colum "year_jday" from floating point number to string
swsm_data$year_jday <- as.character(swsm_data$year_jday)

# Split column "year_jday" into "year" and "jday"
data_split <- unlist(strsplit("2012.125", "\\."))
year <- data_split[1]
jday <- data_split[2]

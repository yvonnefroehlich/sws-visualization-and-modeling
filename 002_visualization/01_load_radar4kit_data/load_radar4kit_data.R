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
#   Geophysical Journal International.
#   https://doi.org/10.1093/gji/ggae245.
# - Ritter J. R. R., Fröhlich Y., Sanz Alonso Y. & Grund M. (2022).
#   Short-scale laterally varying SK(K)S shear wave splitting at BFO,
#   Germany – implications for the determination of anisotropic structures.
#   Journal of Seismology, 26, 1137-1156.
#   https://doi.org/10.1007/s10950-022-10112-w.
# -----------------------------------------------------------------------------
# created: 2024/07/08
# author: Yvonne Fröhlich
# contact: yvonne.froehlich@kit.edu
# #############################################################################


library(dplyr)
library(tidyr)

# Example
root_path = ""
header = 15  # Number of header lines, given in README
station = "BFO"
network = "GR"
obstyp = "NULLS"

# Load CSV file into dataframe
df_swsm <- read.csv(
  file=paste(root_path,"splitresults_",obstyp,"_goodfair_",network,"_",station,".csv", sep=""),
  sep=";",
  skip=header,
)

# Convert from column "year_jday" from floating point number to string
df_swsm$year_jday <- sprintf("%0.3f", df_swsm$year_jday)
df_swsm$year_jday <- as.character(df_swsm$year_jday)

# Split column "year_jday" into "year" and "jday"
df_swsm <- df_swsm %>% separate(year_jday, c("year", "jday"), remove=FALSE)

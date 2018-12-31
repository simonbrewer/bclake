###############################################################################
## reGrid.R
## 
## Script to regrid the output of getClimDailyGrid.R onto the dem to be used
## for modeling
##
###############################################################################

require(raster)

dem.r = raster("dem/bclake_dem.nc")

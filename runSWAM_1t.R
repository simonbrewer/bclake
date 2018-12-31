###############################################################################
## Version of the SWAM model using hydroCA to move water around
##
## Ref: Guidolin et al. (2016). A weighted cellular automata 2D inundation model 
## for rapid flood analysis Env. Mod. Soft., 84, 378-394
##
###############################################################################

###############################################################################
## Libraries
library(rgdal)
library(raster)
library(RColorBrewer)
source("helpers.R")

## Files
dem.r = raster("./dem/bclake_dem.nc")
bas.r = raster("./dem/bclake_bas.nc")
## Clip lake basin
mask.r = bas.r == 19238

###############################################################################
## Estimate cell areas
area.r = area(dem.r) * 1e6
## Calculate slope
slope.r = terrain(dem.r, "slope", unit="degrees")
## Distance between cells (set as constant)
## Needs to be calculated in fortran from coordinates
cellem = 857 ## Approximately 857m cell centers

###############################################################################
## Grid sizes for outout
gridx = dim(dem.r)[1]
gridy = dim(dem.r)[2]

###############################################################################
## Assign outlet cell
out.x = 20.364
out.y = -30.88
out.sp = SpatialPoints(cbind(out.x,out.y))
out.cell <- cellFromXY(dem.r, out.sp)
outlet.r = setValues(dem.r, 0)
outlet.r[out.cell] <- 1
plot(outlet.r)

###############################################################################
## Time steps
delt = 60*60 ## Hourly integration

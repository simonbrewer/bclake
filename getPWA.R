###############################################################################
## Potential water area estimates
##
###############################################################################

require(raster)
require(rgdal)
source("helpers.R")

dem.r = raster("~/Documents/grassdata/hydrosheds/bclake/af_dem_30g.nc")
ldd.r = raster("~/Documents/grassdata/hydrosheds/bclake/af_ldd_30g.nc")

myext = extent(c(19.9,20.3,-30.95,-30.6))
plot(crop(dem.r, myext))

dem.r = crop(dem.r, myext)
ldd.r = crop(ldd.r, myext)
###############################################################################
## Lake border
bclake = readOGR("~/Documents/grassdata/hydrosheds/bclake/929m lake level.kml")

plot(dem.r)
plot(bclake, add=TRUE)

#dem.r = crop(dem.r, bclake)
#ldd.r = crop(ldd.r, bclake)
## r.watershed directions
## "Provides the "aspect" for each cell measured CCW from East. 
## Multiplying positive values by 45 will give the direction in degrees 
## that the surface runoff will travel from that cell."
# 1: 45:  NE
# 2: 90:  N
# 3: 135: NW
# 4: 180: W
# 5: 225: SW
# 6: 270: S
# 7: 315: SE
# 8: 360: E

offx = c(1, 0, -1, -1, -1, 0, 1, 1) ## cols
offy = c(-1, -1, -1, 0, 1, 1, 1, 0) ## rows

mask.r = setValues(dem.r, 1)
gridx = nrow(dem.r)
gridy = ncol(dem.r)
dem = as.matrix(dem.r)
ldd = as.matrix(ldd.r)
mask = as.matrix(mask.r)

pwa.out = fpwa(gridx, gridy, dem, ldd, mask)

drain.r = setValues(dem.r, matrix(pwa.out$drain, 
                                nrow=dim(dem.r)[1], ncol=dim(dem.r)[2]))
plot(drain.r)
plot(bclake, add=TRUE)

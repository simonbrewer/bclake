###############################################################################
## reGrid.R
## 
## Script to regrid the output of getClimDailyGrid.R onto the dem to be used
## for modeling
##
###############################################################################

require(raster)

dem.r = raster("dem/bclake_dem.nc")

## PRE
dpre.cru = brick("./inputs/dpre.nc")
dpre.bc = resample(dpre.cru, dem.r)
writeRaster(dpre.bc, "./inputs/dpre_cru.nc", format="CDF", overwrite=TRUE)

## PET
dpet.cru = brick("./inputs/dpet.nc")
dpet.bc = resample(dpet.cru, dem.r)
writeRaster(dpet.bc, "./inputs/dpet_cru.nc", format="CDF", overwrite=TRUE)

## EVP
devp.cru = brick("./inputs/devp.nc")
devp.bc = resample(devp.cru, dem.r)
writeRaster(devp.bc, "./inputs/devp_cru.nc", format="CDF", overwrite=TRUE)

## CN
dcn.cru = brick("./inputs/dcn.nc")
dcn.bc = resample(dcn.cru, dem.r)
writeRaster(dcn.bc, "./inputs/dcn_cru.nc", format="CDF", overwrite=TRUE)
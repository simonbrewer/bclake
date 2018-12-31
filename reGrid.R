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
dpre.cru = brick("./inputs/dpre_cru.nc")
dpre.bc = resample(dpre.cru, dem.r)
writeRaster(dpre.bc, "./inputs/dpre_bc.nc", format="CDF", overwrite=TRUE)

## PET
dpet.cru = brick("./inputs/dpet_cru.nc")
dpet.bc = resample(dpet.cru, dem.r)
writeRaster(dpet.bc, "./inputs/dpet_bc.nc", format="CDF", overwrite=TRUE)

## EVP
devp.cru = brick("./inputs/devp_cru.nc")
devp.bc = resample(devp.cru, dem.r)
writeRaster(devp.bc, "./inputs/devp_bc.nc", format="CDF", overwrite=TRUE)

## CN
dcn.cru = brick("./inputs/dcn_cru.nc")
dcn.bc = resample(dcn.cru, dem.r)
writeRaster(dcn.bc, "./inputs/dcn_bc.nc", format="CDF", overwrite=TRUE)
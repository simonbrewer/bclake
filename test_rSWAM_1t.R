###############################################################################
## Version of the SWAM model based on original code
##
## This is prototype code in R for getting expected valuesx
##
## Ref: M.T. Coe (1998). A linked global model of terrestrial hydrologic processes:
## Simulation of modern rivers, lakes, and wetlands. JGR, 103, D8, 8885-8899
##
###############################################################################

###############################################################################
## Libraries
library(rgdal)
library(raster)
library(RColorBrewer)
source("helpers.R")

###############################################################################
## Lake border
bclake = readOGR("~/Documents/grassdata/hydrosheds/bclake/929m lake level.kml")

## Lake center 
lon = 20.0833
lat = -30.75
bccent = SpatialPoints(cbind(lon,lat))

###############################################################################
## MODEL SETUP
## Parameters
delt = 60*60 ## Time step (s)
deltu = 24 ## Number of time steps to run model for
bpf = 0 ## Proportion of runoff to put in baseflow [0-1]
effvol = 0.0003 

## Base files
dem.r = raster("./dem/bclake_dem.nc")
bas.r = raster("./dem/bclake_bas.nc")
ldd.r = raster("~/Documents/grassdata/hydrosheds/bclake/af_ldd_30g.nc")

myext = extent(c(20.01,20.03,-30.725,-30.72))
plot(crop(dem.r, myext))
text(crop(dem.r, myext))

dem.r = crop(dem.r, extent(myext))
bas.r = crop(bas.r, extent(myext))
ldd.r = crop(ldd.r, extent(myext))

## Clip lake basin
mask.r = bas.r == 19238

###############################################################################
## Estimate cell areas
area.r = area(dem.r) * 1e6

## Distance between cells (set as constant)
## Needs to be calculated in fortran from coordinates
dist.r =  setValues(dem.r, 857) ## Approximately 857m cell centers

###############################################################################
## Grid sizes for outout
gridx = dim(dem.r)[2]
gridy = dim(dem.r)[1]

###############################################################################
## Assign outlet cell
out.x = 20.364
out.y = -30.88
out.sp = SpatialPoints(cbind(out.x,out.y))
out.cell <- cellFromXY(dem.r, out.sp)
outlet.r = setValues(dem.r, 0)
outlet.r[out.cell] <- 1

###############################################################################
## Forcing data
# dpre.stk = brick("./inputs/dpre_bc.nc") * 5
# dpet.stk = brick("./inputs/dpet_bc.nc")
# devp.stk = brick("./inputs/devp_bc.nc")
# dcn.stk = brick("./inputs/dcn_bc.nc")
dpre.stk = crop(brick("./inputs/dpre_bc.nc"),myext) * 10
dpet.stk = crop(brick("./inputs/dpet_bc.nc"),myext)
devp.stk = crop(brick("./inputs/devp_bc.nc"),myext)
dcn.stk = crop(brick("./inputs/dcn_bc.nc"),myext)

###############################################################################
## Calculate runoff
dro.stk = (dpre.stk + dcn.stk) - dpet.stk

###############################################################################
## Water storage rasters
wse.r = wvl.r = war.r = setValues(dem.r, 0)

###############################################################################
## Convert to matrices
dem = as.matrix(dem.r)
ldd = as.matrix(ldd.r)
ldd[,] <- 8 ## Constant flow direction
outelev = as.matrix(dem.r)  # Needs changing
mask = as.matrix(mask.r)
cella = as.matrix(area.r)
celld = as.matrix(dist.r)
wvl = as.matrix(wvl.r)
wse = as.matrix(wse.r)
war = as.matrix(war.r)

offx = c(1, 0, -1, -1, -1, 0, 1, 1) ## cols
offy = c(-1, -1, -1, 0, 1, 1, 1, 0) ## rows

## Test with 5m everywhere
# wse = wse + 5

###############################################################################
cols <- colorRampPalette(brewer.pal(9,"Blues"))(100)

###############################################################################
nyrs = 5
bclevel = matrix(NA, nrow=365, ncol=nyrs)
## Convert forcing to matrices
for (yr in 1:nyrs) {
  
  for (dy in 1:179) {
    print(paste(yr,"Lake level:", bclevel[dy,yr]))
    print(paste("Doing",yr,dy))
    
    pre = as.matrix(raster(dpre.stk, dy))
    evp = as.matrix(raster(devp.stk, dy))
    ro = clamp(raster(dro.stk, dy), lower=0, useValues=TRUE)
    ro = as.matrix(ro)
    print(paste("ro",ro))
    sro = ro * (1-bpf)
    bro = ro * bpf
    
    fin = fout = matrix(0 ,gridy, gridx)
    drain = 0
    
    for (k in 1:24) {
      for (i in 1:gridx) {
        for (j in 1:gridy) {
          
          if (mask[j,i] == 1) {
            edge = 0 # Set edge
            # Find target downstream cell
            move = ldd[j,i]
            ii = i + offx[move]
            jj = j + offy[move]
            
            # Test to see if we've left the basin
            if ((ii < 1) | (ii > gridx) | (jj < 1) | (jj > gridy)) {
              edge = 1
            }

            # Test to see if this is a coastal cell
            if (ldd[j,i] < 0) {
              edge = 1
            }

            # Calculate outflow
            # Get wei and wed
            if ((wse[j,i]) > 0.) {
              # If not an edge, then get wel and wed (water elevations)
              if (edge == 0)  {
                wel = wse[j,i] + dem[j,i]
                wed = wse[jj,ii] + dem[jj,ii]
                fout[j,i] = max((wel-wed)*cella[j,i],0.) * effvol/celld[j,i]
                fin[jj,ii] = fin[jj,ii] + fout[j,i]
              } else { # If at edge account for water loss in drain reservoir
                wel = wse[j,i] + dem[j,i]
                wed = wel - 1 # constant drop for edge cells (needs fixing)
                fout[j,i] = max((wel-wed)*cella[j,i],0.) * effvol/celld[j,i]
                drain = drain + fout[j,i]
              }
            } 
            

          } # Mask loop
        } # y loop
      } # x loop
      # print(fin)
      # print(fout)
      # if (min(fout) > 0) stop()
      
###############################################################################
## Update linear reservoir
      for (i in 1:gridx) {
        for (j in 1:gridy) {
          
          if (mask[j,i] == 1) {
            ro_tmp = ((sro[j,i] * 1e-3) / (60*60*24)) * delt * cella[j,i]
            bf_tmp = ((bro[j,i] * 1e-3) / (60*60*24)) * delt * cella[j,i]
            ppt_tmp = ((pre[j,i] * 1e-3) / (60*60*24)) * delt * cella[j,i]
            evap_tmp = ((evp[j,i] * 1e-3) / (60*60*24)) * delt * cella[j,i]
            
            dwv = ( (ro_tmp + bf_tmp ) * (1 - war[j,i])  ) +
              ( (ppt_tmp - evap_tmp) * war[j,i]  ) +
              ( fin[j,i] - fout[j,i] )
            
            wvl[j,i] = max(c(wvl[j,i] + dwv, 0))
            wse[j,i] = wvl[j,i] / cella[j,i]
            
            # if (wse[j,i] > 0) {
            #   war[j,i] = min(c(wse[j,i], 1))
            # }
            
          } # Mask loop
        } # y loop
      } # x loop
      
    } # Hourly loop
  } # Daily loop
    stop()
  
} # Yearly loop

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

## Base files
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
## Forcing data
dpre.stk = brick("./inputs/dpre_bc.nc") #* 5
dpet.stk = brick("./inputs/dpet_bc.nc")
devp.stk = brick("./inputs/devp_bc.nc")
dcn.stk = brick("./inputs/dcn_bc.nc")

###############################################################################
## Calculate runoff
dro.stk = (dpre.stk + dcn.stk) - dpet.stk

###############################################################################
## Grid to record total inflow and outflow from each timestep
itot.r = setValues(dem.r, 0)
otot.r = setValues(dem.r, 0)

###############################################################################
## Convert to matrices
dem = as.matrix(dem.r)
cella = as.matrix(area.r)
mask = as.matrix(mask.r)
outlet = as.matrix(outlet.r)
itot = as.matrix(itot.r)
otot = as.matrix(otot.r)
# wse = as.matrix(dem.r)

## Test with 5m everywhere
# wse = wse + 5

###############################################################################
cols <- colorRampPalette(brewer.pal(9,"Blues"))(100)

###############################################################################
nyrs = 20
bclevel = matrix(NA, nrow=365, ncol=nyrs)
## Convert forcing to matrices
for (j in 1:nyrs) {

  print(paste(j,"Lake level:", bclevel[i,j]))
  for (i in 1:365) {
    # print(paste("Doing",i))
    
    pre = as.matrix(raster(dpre.stk, i))
    evp = as.matrix(raster(devp.stk, i))
    ro = clamp(raster(dro.stk, i), lower=0, useValues=TRUE)
    ro = as.matrix(ro)
    sro = ro * (1-bpf)
    bro = ro * bpf
    sim.out = swamCA_1t(gridx, gridy, dem, mask, cella, outlet,
                        pre, evp, sro, bro,
                        wse, otot, itot, delt, deltu,
                        mannN=0.05, cellem=cellem,
                        tolwd=0.0001, tolslope=0.001)
    wse = sim.out$wse
    # print(sum(sim.out$wse-sim.out$dem))
    wse.r = setValues(dem.r, matrix(sim.out$wse - sim.out$dem, 
                                    nrow=dim(dem.r)[1], ncol=dim(dem.r)[2]))
    # print(paste(i,"Max depth:", cellStats(wse.r, max)))
    # plot(wse.r, col=cols)
    
    bclevel[i,j] = extract(wse.r, bccent)
  }
  
}
plot(wse.r, col=cols)
plot(bclake, add=TRUE)

plot(log10(wse.r), col=cols)
plot(bclake, add=TRUE)

plot(c(bclevel), type='l')

pit.r = focal(dem.r, w=matrix(1,3,3), findPit)
pit.sp = data.frame(coordinates(pit.r), pit = getValues(pit.r))
pit.sp = subset(pit.sp, pit==1)
coordinates(pit.sp) <- ~x+y 
plot(wse.r)
plot(log10(wse.r), col=cols)
plot(bclake, add=TRUE)
plot(pit.sp, add=TRUE)

## Output

bclevel.ts = ts(c(bclevel), start = 1, freq=365)
pdf("bclevel_test.pdf")
plot(bclevel.ts, xlab = "Time (yrs)", ylab="Lake level (m)")
dev.off()


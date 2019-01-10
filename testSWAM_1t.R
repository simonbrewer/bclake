###############################################################################
## Version of the SWAM model based on original code
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
effvol = 0.003 

## Base files
dem.r = raster("./dem/bclake_dem.nc")
bas.r = raster("./dem/bclake_bas.nc")
ldd.r = raster("~/Documents/grassdata/hydrosheds/bclake/af_ldd_30g.nc")
ldd.r = crop(ldd.r, extent(dem.r))

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

###############################################################################
## Forcing data
dpre.stk = brick("./inputs/dpre_bc.nc") * 5
dpet.stk = brick("./inputs/dpet_bc.nc")
devp.stk = brick("./inputs/devp_bc.nc")
dcn.stk = brick("./inputs/dcn_bc.nc")

###############################################################################
## Calculate runoff
dro.stk = (dpre.stk + dcn.stk) - dpet.stk

###############################################################################
## Water storage rasters
wse.r = wvl.r = war.r = setValues(dem.r, 0)

###############################################################################
## Convert to matrices
dem = as.matrix(dem.r)
ldd = as.matrix(ldd.r)  # Needs changing
outelev = as.matrix(dem.r)  # Needs changing
mask = as.matrix(mask.r)
cella = as.matrix(area.r)
celld = as.matrix(dist.r)
wvl = as.matrix(wvl.r)
wse = as.matrix(wse.r)
war = as.matrix(war.r)

## Test with 5m everywhere
# wse = wse + 5

###############################################################################
cols <- colorRampPalette(brewer.pal(9,"Blues"))(100)

###############################################################################
nyrs = 5
bclevel = matrix(NA, nrow=365, ncol=nyrs)
## Convert forcing to matrices
for (j in 1:nyrs) {
  
  print(paste(j,"Lake level:", bclevel[i,j]))
  for (i in 1:365) {
    print(paste("Doing",j,i))
    
    pre = as.matrix(raster(dpre.stk, i))
    evp = as.matrix(raster(devp.stk, i))
    ro = clamp(raster(dro.stk, i), lower=0, useValues=TRUE)
    ro = as.matrix(ro)
    sro = ro * (1-bpf)
    bro = ro * bpf
    sim.out = swam_1t(gridx, gridy, dem, ldd, outelev, 
                      mask, cella, celld,
                      pre, evp, sro, bro,
                      delt, deltu, effvol,
                      wvl, wse, war)
    
    wvl = sim.out$wvl
    wse = sim.out$wse
    # print(sum(sim.out$wse-sim.out$dem))
    print(max(sim.out$fout))
    wse.r = setValues(dem.r, matrix(sim.out$wse, 
                                    nrow=dim(dem.r)[1], ncol=dim(dem.r)[2]))
    # print(paste(i,"Max depth:", cellStats(wse.r, max)))
    # plot(wse.r, col=cols)
    
    bclevel[i,j] = extract(wse.r, bccent)
  }
  
}
plot(wse.r, col=cols)
plot(bclake, add=TRUE)

war.r = setValues(dem.r, matrix(sim.out$war, 
                                nrow=dim(dem.r)[1], ncol=dim(dem.r)[2]))
plot(war.r, col=cols)
plot(bclake, add=TRUE)
stop()
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


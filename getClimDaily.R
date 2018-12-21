## Extract climate data and run splash
require(raster)
require(rgdal)
require(Evapotranspiration)
source("./helpers.R")

## Location
lon = 20.0833
lat = -30.75
yr = 1975
elv = 929
whc = 150.0
mdays = c(31,28,31,30,31,30,31,31,30,31,30,31)

## Lake center 
bccent = SpatialPoints(cbind(lon,lat))

## Lake polygon
bclake = readOGR("~/Documents/grassdata/hydrosheds/bclake/929m lake level.kml")
plot(bclake)
plot(bccent, add=TRUE)

## Get monthly series
tmp.r = stack("./cru/cl/bc_cru_10min_tmp.nc")
dtr.r = stack("./cru/cl/bc_cru_10min_dtr.nc")
pre.r = stack("./cru/cl/bc_cru_10min_pre.nc")
rhm.r = stack("./cru/cl/bc_cru_10min_reh.nc")
sun.r = stack("./cru/cl/bc_cru_10min_sun.nc")

## Region
tmp.bc = extract(tmp.r, bccent) #- 273.15
dtr.bc = extract(dtr.r, bccent) 
tmn.bc = tmp.bc - (dtr.bc/2)
tmx.bc = tmp.bc + (dtr.bc/2)
pre.bc = extract(pre.r, bccent) #* 60*60*24
rhm.bc = extract(rhm.r, bccent) 
sun.bc = extract(sun.r, bccent) / 100
# sf.bc = 1 - (cld.bc / 100)

## Get daily values
dtmp = daily(c(tmp.bc))$dly
dtmn = daily(c(tmn.bc))$dly
dtmx = daily(c(tmx.bc))$dly
dpre = daily(c(pre.bc/mdays))$dly
drhm = daily(c(rhm.bc))$dly
dsun = daily(c(sun.bc))$dly

## sunp is listed as 
## sunp	sunshine percent of maximum possible (percent of daylength)

ndays = length(dtmp)
dpet = rep(NA, ndays)

## Run SPLASH
aetpet.df = splashf(dtmp,dpre,dsun,lat=lat,
                      yr=yr, elv=elv)
  
## Estimate Tdew
## From humidity.to.dewpoint from the weathermetrics package
## dewpoint <- (rh/100)^(1/8) * (112 + (0.9 * t)) - 112 + (0.1 * t)
## This is originally from the source code for the 
## US National Weather Service's online heat index calculator.

ddew = (drhm/100)^(1/8) * (112 + (0.9 * dtmp)) - 112 + (0.1 * dtmp)

## Run Morton CRWE
bcdates = seq.Date(as.Date("1975-01-01"), length.out = 365, by = 1)

## Using n for sunshine hours
clim.df = data.frame(Year = as.numeric(format(bcdates, "%Y")),
                     Month = as.numeric(format(bcdates, "%m")),
                     Day = as.numeric(format(bcdates, "%d")),
                     Tmax = dtmx, Tmin = dtmn, 
                     n = aetpet.df$dsl, Tdew = ddew, Precip = dpre)

clim.in = ReadInputs(varnames = c("Tmax","Tmin","Tdew","n"),
                     clim.df, 
                     stopmissing=c(10,10,3))
data(constants)
constants$lat = lat
constants$lat_rad = lat * pi / 180
constants$Elev = elv
lake.out <- ET.MortonCRWE(clim.in, constants, ts="monthly",
                         est="shallow lake ET", solar="sunshine hours", Tdew= TRUE, 
                         alpha = NULL, message="yes", save.csv="no") 
plot(lake.out$ET.MonthlyAve)
# write.csv(out.df, "bclake.csv", row.names = FALSE)

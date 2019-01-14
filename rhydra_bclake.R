## Test case (runs Chad example)

require(raster)
source("helpers.R")

load("./data/chadExample_hydra.RData")

## Parameters
nyrs = 1 # of years to run the model after January 1937
startyear = 1 # year to start run (from January 1937) 
converg = 1 # 0 if set convergence helper, 1 if not
laket = 0 # 0 if predict lakes, 1 if parameterized with Cogley obs, 2 if with previously simulated
spin = 30  # number of spin-up years to run before reading transient data
normal = 1 # 0 if use normalization, 1 if not
leap = 2 # of years to the first leap year
irrig = 1 # 0 if use irrigation, function 1 if not

## Grid dimensions
gridxf = 216
gridyf = 216

hydro.out = rhydra(nyrs, startyear, converg, laket, spin,
                   normal, leap, irrig,
                   outnewi, outnewj, basin, dem, rivdir, mflac,
                   prcpi, evapi, runin, drainin,
                   gridxf, gridyf)

save(hydro.out, file=paste0("rhydra_out_",spin,".Rdata"))

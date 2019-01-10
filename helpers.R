###############################################################################
## Helper functions for "runSplash.R"
dyn.load("./fortran/splash.so")
dyn.load("./fortran/swamCA.so")
dyn.load("./fortran/getpwa.so")
dyn.load("./fortran/swam.so")
###############################################################################

###############################################################################
## aetpet: Function to calculate aet and pet
splashf <- function(dtemp,dprec,dsun,lat,yr,elv) {
  dyn.load('fortran/splash.so')
  retdata <- .Fortran("spin_up_sm",
                      yr = as.integer(yr),
                      lat = as.double(lat),
                      elv = as.double(elv),
                      pr = as.double(dprec),
                      tc = as.double(dtemp),
                      sf = as.double(dsun),
                      maet = double(12),
                      mpet = double(12),
                      mcn = double(12),
                      mro = double(12),
                      msm = double(12),
                      daet = double(365),
                      dpet = double(365),
                      dcn = double(365),
                      dro = double(365),
                      dsm = double(365),
                      sm = double(1),
                      ddl = double(365),
                      dsl = double(365))
  return(retdata)
}
###############################################################################

###############################################################################
# ## DAILY: Function to interpolate from monthly to daily
# ## Replace with 'approx'?
# 
daily <- function(mly) {
  retdata <- .Fortran("daily",
                      mly = as.double(mly),
                      dly = double(365))
  
  return(retdata)
}
###############################################################################

###############################################################################
## SWAM CA model code
swamCA_1t <- function(gridx, gridy, dem, mask, cella, outlet,
                      ppt, evap, runoff, baseflow,
                      wse, otot, itot, delt, deltu,
                      mannN=0.05, cellem=50, 
                      tolwd=0.0001, tolslope=0.001) {
  
  simcf = .Fortran("swamca_1t",
                   m = as.integer(gridx), n = as.integer(gridy),
                   dem = as.double(dem), 
                   mask = as.integer(mask), 
                   cella = as.double(cella), 
                   outlet = as.integer(outlet), 
                   ppt = as.double(ppt),
                   evap = as.double(evap),
                   runoff = as.double(runoff),
                   baseflow = as.double(baseflow),
                   wse = as.double(wse),
                   otot = as.double(dem), itot = as.double(dem),
                   dt = as.double(delt), dtu = as.integer(deltu), 
                   mannn = as.double(mannN), cellx = as.double(cellem),
                   cellem = as.double(cellem), 
                   tolwd = as.double(tolwd), tolslope = as.double(tolslope))
  return(simcf)
  
}
###############################################################################

###############################################################################
## SWAM original model code
swam_1t <- function(gridx, gridy, dem, ldd, outelev, 
                    mask, cella, celld,
                    ppt, evap, runoff, baseflow,
                    delt, deltu, effvol,
                    wvl, wse, war) {
  
  simcf = .Fortran("swam_1t",
                   m = as.integer(gridx), n = as.integer(gridy),
                   dem = as.double(dem), 
                   ldd = as.integer(ldd), 
                   outelv = as.double(outelev), 
                   mask = as.integer(mask), 
                   cella = as.double(cella), 
                   celld = as.double(celld), 
                   ppt = as.double(ppt),
                   evap = as.double(evap),
                   runoff = as.double(runoff),
                   baseflow = as.double(baseflow),
                   dt = as.double(delt), dtu = as.integer(deltu), 
                   u = as.double(effvol), 
                   wvl = as.double(wvl), 
                   wse = as.double(wse), 
                   war = as.double(war), 
                   fin = double(gridx*gridy), 
                   fout = double(gridx*gridy))
  return(simcf)
  
}
###############################################################################

###############################################################################
## Potential water area code (fortran)
fpwa <- function(gridx, gridy, dem, ldd, mask) {
  
  pwaout = .Fortran("getpwa",
                   m = as.integer(gridx), n = as.integer(gridy),
                   dem = as.double(dem), 
                   ldd = as.integer(ldd),
                   mask = as.integer(mask),
                   pwa = integer(gridx*gridy),
                   drain = integer(gridx*gridy),
                   outelev = double(gridx*gridy),
                   iout = integer(gridx*gridy),
                   jout = integer(gridx*gridy)
                   )
  return(pwaout)
  
}
###############################################################################

###############################################################################
## binOutlet
## Function to find outlet cells
## Defined as cells next to the border mask, with no lower elevation neighbors
binOutlet <- function (x) {
  outlet = FALSE
  nna = length(is.na(x))
  if (length(nna) >= 1) { ## Are we next to an edge?
    if (!is.na(x[5])) { ## Is the cell an edge cell?
      if (which.min(x) == 5) {
        outlet = TRUE
      }
    } 
  } 
  return(outlet)
}

## Older version
# binOutlet <- function (x) {
#   outlet = FALSE
#   if (max(x) == 1e6) { ## Are we next to an edge?
#     if (x[5] != 1e6) { ## Is the cell an edge cell?
#       if (which.min(x) == 5) {
#         outlet = TRUE
#       }
#     } 
#   } 
# }
###############################################################################

###############################################################################
## Function to convert radians to degrees
rad2deg <- function(x) {
  x*180/pi
}
###############################################################################

###############################################################################
## Function to convert radians to degrees
deg2rad <- function(x) {
  x/180*pi
}
###############################################################################

###############################################################################
## Function to calculate great circle distances
gcDist <- function(lon1, lat1, lon2, lat2, r=6378) {
  ## Degree conversion
  lon1 <- lon1 * pi/180
  lon2 <- lon2 * pi/180
  lat1 <- lat1 * pi/180
  lat2 <- lat2 * pi/180
  ## Central angle
  ca <- acos((sin(lat1)*sin(lat2)) + 
               (cos(lat1)*cos(lat2) * cos(abs(lon1-lon2))))
  d <- r * ca
  return(d)
}
###############################################################################

###############################################################################
## Function to lower border next to outlets
modBorder <- function(outlet, dem, mask) {
  require(geosphere)
  oID <- Which(outlet==1, cells=TRUE)
  if (length(oID) > 0) {
    for (i in 1:length(oID)) {
      cen.crds = xyFromCell(mask, oID[i])
      ngb.rc = adjacent(mask, oID[i], directions = 8)
      ngb.crds = xyFromCell(mask, ngb.rc[,2])
      ngb.vals = extract(mask, ngb.rc[,2])
      ngb.crds <- ngb.crds[which(ngb.vals==1),]
      ngb.rc <- ngb.rc[which(ngb.vals==1),]
      ngb.dist = distCosine(cen.crds, ngb.crds)
      
      bID <- which.min(ngb.dist)
      dem[ngb.rc[bID,2]] <- dem[ngb.rc[bID,2]]*-1
      mask[ngb.rc[bID,2]] <- mask[ngb.rc[bID,2]]*-1
      
    }
    
  }
  return(list(dem=dem,mask=mask))
}
###############################################################################

###############################################################################
## Function to find pits in DEM
findPit <- function(x) {
  pit = FALSE
  if (which.min(x) == 5) {
    pit = TRUE
  }
  return(pit)
}
###############################################################################

############################################################################
# Matrix manipulation methods
#
# For simplicity we have avoided to create generic functions for 'flip' etc.
# and therefore we have to call the corresponding methods coupled to the
# 'matrix' class explicitly, i.e. flip.matrix().
############################################################################
# Flip matrix (upside-down)
flip.matrix <- function(x) {
  mirror.matrix(rotate180.matrix(x))
}

# Mirror matrix (left-right)
mirror.matrix <- function(x) {
  xx <- as.data.frame(x);
  xx <- rev(xx);
  xx <- as.matrix(xx);
  xx;
}

# Rotate matrix 90 clockworks
rotate90.matrix <- function(x) {
  t(mirror.matrix(x))
}

# Rotate matrix 180 clockworks
rotate180.matrix <- function(x) { 
  xx <- rev(x);
  dim(xx) <- dim(x);
  xx;
}

# Rotate matrix 270 clockworks
rotate270.matrix <- function(x) {
  mirror.matrix(t(x))
}


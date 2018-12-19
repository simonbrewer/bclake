library(raster)

crudir = "~/Dropbox/Data/climate/cru_cl_2.00/"
filelist = list.files(crudir, pattern=".nc")

r = raster(paste0(crudir,"cru_10min_elv.nc"))

## The. Worst. Biome. Ever….
## 19E – 27.5E
## 27.5S – 33S
# myext = extent(c(10,40,-35,-10))
myext = extent(c(19,27.5,-33,-27.5))

for (i in 1:length(filelist)) {
  varname = unlist(strsplit(unlist(strsplit(filelist[i], "_"))[3],"\\."))[1]
  outfile = paste0("bc_cru_10min_",varname,".nc")
  tmp.stk = stack(paste0(crudir,filelist[i]))
  tmp.stk.2 = crop(tmp.stk, myext)
  
  writeRaster(tmp.stk.2, filename = outfile, varname=varname)
}

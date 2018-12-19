library(raster)

filelist = list.files(".", pattern=".nc")

r = raster("cru_ts4.02.1901.2017.pre.dat.nc")

myext = extent(c(10,40,-35,-10))

for (i in 1:length(filelist)) {
  varname = unlist(strsplit(filelist[i], "\\."))[5]
  outfile = paste0("safr_1901_2017_",varname,".nc")
  tmp.stk = stack(filelist[i])
  tmp.stk.2 = crop(tmp.stk, myext)
  writeRaster(tmp.stk.2, filename = outfile, varname=varname)
}

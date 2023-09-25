# Simulated example with leaf-level data
library(terra)
library(ncdf4)
library(dplyr)
source(file.path(Sys.getenv('specdiv'),'specdiv_functions.R'))


# LPJ data ----------------------------------------------------------------
lpj.path <- Sys.getenv('lpjpath')
dp <- Sys.getenv('datapath')

lpj.nc <- Sys.getenv('lpjnc')
shp.path <- Sys.getenv('shpfile')

# LPJ array ---------------------------------------------------------------


shp <- vect(shp.path)
crs(shp) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
print(shp)
lpj.nc <- nc_open(file.path(lpj.path,lpj.nc)) %>%
    ncvar_get(varid = 'DR') %>%
    aperm(c(2,1,3,4))


for (m in seq(12)) {
    
    # Read in the data
    lpj.r <- lpj.nc[,,,m] %>% rast(extent = c(-180,180,-90,90)) 
    crs(lpj.r) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
    r <- lpj.r %>%
        crop(shp) %>%
        mask(shp)
   print(r) 
    lpj.r[lpj.r==0] <- NA
    SpectralBands <- seq(400,2500,10)
    BandNames <- paste('Band', SpectralBands, sep = '_')
    names(r) <- BandNames
    
    # Apply specdiv algorithm
    sdiv <- specdiv(r, fact = 1) # No alpha diversity
    
    
    if (m == 1) {
        lpj.beta.lcsd <- sdiv$rasters$beta_lcsd
        lpj.beta.lcss <- sdiv$rasters$beta_lcss
        
    } else {
        lpj.beta.lcsd <- c(lpj.beta.lcsd, sdiv$rasters$beta_lcsd)
        lpj.beta.lcss <- c(lpj.beta.lcss, sdiv$rasters$beta_lcss)
        
    }
    
} # ..end month

mons <- c('Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec')
names(lpj.beta.lcsd) <- mons
names(lpj.beta.lcss) <- mons

writeRaster(lpj.beta.lcsd, file.path(dp, 'lpj_prosail_Version_021_global_lcsd_2022.tif')
writeRaster(lpj.beta.lcsd, file.path(dp, 'lpj_prosail_Version_021_global_lcss_2022.tif')

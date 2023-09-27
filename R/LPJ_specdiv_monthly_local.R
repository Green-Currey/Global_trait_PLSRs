# Simulated example with leaf-level data
library(terra)
library(ncdf4)
library(dplyr)
source(file.path('~/Other Projects/specdiv/functions/specdiv_functions.R'))



# LPJ data ----------------------------------------------------------------
shp <- vect('~/Current Projects/SBG/LPJ/Misc data/Shapefile/Would_boundaries_noGreenland.shp')
lpj.path <- '~/Current Projects/SBG/LPJ/Reflectance_Data/version2'
lpj.nc <- nc_open(list.files(lpj.path, pattern = '_DR_', full.names = T)[1]) %>%
    ncvar_get(varid = 'DR')


# LPJ array ---------------------------------------------------------------

crs(shp) <- 'EPSG:4326'

lpj.nc <- nc_open(list.files(lpj.path, pattern = '_DR_', full.names = T)[1]) %>%
    ncvar_get(varid = 'DR') %>%
    aperm(c(2,1,3,4))


for (m in 7) {
    
    # Read in the data
    lpj.r <- lpj.nc[,,,m] %>% rast(extent = c(-180,180,-90,90), crs = 'EPSG:4326') 
    r <- lpj.r %>%
        crop(shp) %>%
        mask(shp)
    
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

plot(lpj.beta.lcsd)
plot(lpj.beta.lcss)

mons <- c('Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec')
names(lpj.beta.lcsd) <- mons
names(lpj.beta.lcss) <- mons

writeRaster(lpj.beta.lcsd, file.path(dp, 'lpj_prosail_Version_021_global_lcsd_2022.tif'))
writeRaster(lpj.beta.lcsd, file.path(dp, 'lpj_prosail_Version_021_global_lcss_2022.tif'))
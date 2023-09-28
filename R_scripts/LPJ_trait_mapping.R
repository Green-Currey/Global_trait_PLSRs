library(terra)
library(ncdf4)
library(tidyverse)
source('/discover/nobackup/bcurrey/Global_trait_PLSRs/R_scripts/lpj_plsr_functions.R')


print('Starting trait mapping')
t <- Sys.time()
# paths -------------------------------------------------------------------
lpj.path <- Sys.getenv('lpjpath')
dp <- Sys.getenv('datapath')
lpj.nc <- Sys.getenv('lpjnc')

print(lpj.path)
print(dp)
print(lpj.nc)

# var names
lma.name <- 'LDMC_3km_v1'
n.name <- 'LNC_3km_v1'
p.name <- 'LPC_3km_v1'
sla.name <- 'SLA_3km_v1'


# create PLSR data.frame --------------------------------------------------
print('Reading in LPJ array')
lpj.array <- nc_open(file.path(lpj.path, lpj.nc)) %>%
    ncvar_get('DR', start = c(1,1,1,7), count = c(-1,-1,-1,1)) %>%
    aperm(c(2,1,3))

lpj.r <- rast(lpj.array, ext = c(-180,180,-90,90))
crs(lpj.r) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"

lpj.r[lpj.r>0.7] <- NA
lpj.r[lpj.r==0] <- NA

cells <- vect(crds(lpj.r, na.rm = F))
crs(cells) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"


print('Extract points from TRY trait maps')
ldmc <- rast(file.path(dp, paste0(lma.name,'.tif')))
ldmc[ldmc<0.1 | ldmc>0.5] <- NA
ldmc <- ldmc %>% terra::extract(cells, ID = F)

lnc <- rast(file.path(dp, paste0(n.name,'.tif')))
lnc[lnc<13 | lnc>25] <- NA
lnc <- lnc %>% terra::extract(cells, ID = F)

lpc <- rast(file.path(dp, paste0(p.name,'.tif'))) 
lpc[lpc<0.8 | lpc>2.1] <- NA
lpc <- lpc %>% terra::extract(cells, ID = F)

sla <- rast(file.path(dp, paste0(sla.name,'.tif')))
sla[sla<7 | sla>21] <- NA
sla <- sla %>% terra::extract(cells, ID = F)

lpj.df <- as.data.frame(lpj.r, xy = T, na.rm = F)
plsr.data <- cbind.data.frame(lpj.df, ldmc, lnc, lpc, sla) %>% na.exclude()
names(plsr.data) <- c('x', 'y', paste0('wave',seq(400,2500,10)), lma.name, n.name, p.name, sla.name)

# run PLSR ----------------------------------------------------------------

# This is for July
print('LMA PLSR')
lma.coefs <- runPLSR(plsr.data, data.var = lma.name, band.prefix = 'wave', train.size = 5000, plots = F,
                     jk.test = F, jk.prop = 0.15, jk.iterations = 20, jk.comps = 5)
print('N PLSR')
n.coefs <- runPLSR(plsr.data, data.var = n.name, band.prefix = 'wave', train.size = 5000, plots = F,
                   jk.test = F, jk.prop = 0.15, jk.iterations = 20, jk.comps = 5)
print('P PLSR')
p.coefs <- runPLSR(plsr.data, data.var = p.name, band.prefix = 'wave', train.size = 5000, plots = F,
                   jk.test = F, jk.prop = 0.15, jk.iterations = 20, jk.comps = 5)
print('SLA PLSR')
sla.coefs <- runPLSR(plsr.data, data.var = sla.name, band.prefix = 'wave', train.size = 5000, plots = F,
                     jk.test = F, jk.prop = 0.15, jk.iterations = 20, jk.comps = 5)

# writing coeffs
print('Writing coefficients')
coeff.df <- data.frame(coeff = c('Intercept', seq(400,2500,10)), lma = lma.coefs, n = n.coefs, p = p.coefs, sla = sla.coefs)
write_csv(coeff.df, paste0(dp, 'LPJ-PROSAIL_TRY_coeffs_July.csv'))



print('Applying coefficients')
lma.map <- trait.map(lpj.r, coeffs = coeff.df$lma[-1], intercept = coeff.df$lma[1], coeffs_wl = seq(400,2500,10))
n.map <- trait.map(lpj.r,  coeffs = coeff.df$n[-1], intercept = coeff.df$n[1], coeffs_wl = seq(400,2500,10))
p.map <- trait.map(lpj.r, coeff.df$p[-1], coeff.df$p[1], coeffs_wl = seq(400,2500,10))
sla.map <- trait.map(lpj.r, coeff.df$sla[-1], coeff.df$sla[1], coeffs_wl = seq(400,2500,10))


print('Writing raster stack')
trait.stack <- c(lma.map, n.map, p.map, sla.map)
names(trait.stack) <- c(lma.name, n.name, p.name, sla.name)
writeRaster(trait.stack, file.path(dp, 'LPJ_v21_TRY_trait_maps_July_2020.tif'))

print('Finished')
print(Sys.time() - t)

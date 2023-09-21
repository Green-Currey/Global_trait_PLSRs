library(terra)
library(ncdf4)
library(tidyverse)
source('/discover/nobackup/bcurrey/Global_trait_PLSRs/R/lpj_plsr_functions.R')


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
lpj.array <- nc_open(file.path(lpj.path, lpj.nc)) %>% ncvar_get('DR') %>% aperm(c(2,1,3,4))
r <- rast(lpj.array[,,1,1], ext = c(-180,180,-90,90))
crs(r) <-'EPSG:4326'
cells <- crds(r) %>% terra::vect()
crs(cells) <- 'EPSG:4326'

print('Extract points from TRY trait maps')
ldmc <- rast(file.path(dp, paste0(lma.name,'.tif'))) %>% terra::extract(cells, ID = F)
lnc <- rast(file.path(dp, paste0(n.name,'.tif'))) %>% terra::extract(cells, ID = F)
lpc <- rast(file.path(dp, paste0(p.name,'.tif'))) %>% terra::extract(cells, ID = F)
sla <- rast(file.path(dp, paste0(sla.name,'.tif'))) %>% terra::extract(cells, ID = F)

# Create PLSRs for July
lpj.df <- as.data.frame(rast(lpj.array[,,,7]), xy = T)

plsr.df <- cbind.data.frame(lpj.df, ldmc, lnc, lpc, sla) %>% na.exclude()
names(plsr.df) <- c('x', 'y', paste0('wave',seq(400,2500,10)), lma.name, n.name, p.name, sla.name)
# run PLSR ----------------------------------------------------------------

# This is for July
print('LMA PLSR')
lma.coefs <- runPLSR(plsr.df, data.var = lma.name, train.size = 5000, plots = F,
                 jk.prop = 0.15, jk.iterations = 30)
print('N PLSR')
n.coefs <- runPLSR(plsr.df, data.var = n.name, train.size = 5000, plots = F,
                 jk.prop = 0.15, jk.iterations = 30)
print('P PLSR')
p.coefs <- runPLSR(plsr.df, data.var = p.name, train.size = 5000, plots = F,
                 jk.prop = 0.15, jk.iterations = 30)
print('SLA PLSR')
sla.coefs <- runPLSR(plsr.df, data.var = sla.name, train.size = 5000, plots = F,
                 jk.prop = 0.15, jk.iterations = 30)

# writing coeffs
print('Writing coefficients')
coeff.df <- data.frame(coeff = c('Intercept', seq(400,2500,10)), lma = lma.coefs, n = n.coefs, p = p.coefs, sla = sla.coefs)
write_csv(coeff.df, file.path(fp, 'LPJ_TRY_coeffs_July.csv'))


print('Applying coefficients')
lpj.r <- rast(lpj.array[,,,7])
lma.map <- trait.map(lpj.r, coeffs = coeff.df$lma, coeffs_wl = seq(400,2500,10))
n.map <- trait.map(lpj.r, coeffs = coeff.df$n, coeffs_wl = seq(400,2500,10))
p.map <- trait.map(lpj.r, coeffs = coeff.df$p, coeffs_wl = seq(400,2500,10))
sla.map <- trait.map(lpj.r, coeffs = coeff.df$sla, coeffs_wl = seq(400,2500,10))

print('Writing raster stack')
trait.stack <- c(lma.map, n.map, p.map, sla.coefs)
names(trait.stack) <- c(lma.name, n.name, c.name, sla.name)
writeRaster(trait.stack, file.path(dp, 'LPJ_v21_TRY_trait_maps_July_2020.tif'))

print('Finished')
print(Sys.time() - t)

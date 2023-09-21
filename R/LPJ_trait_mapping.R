library(terra)
library(ncdf4)
library(tidyverse)
source('~/Current Projects/SBG/LPJ/Global_trait_PLSRs/R/lpj_plsr_functions.R')



# paths -------------------------------------------------------------------


trait.path <- '~/Current Projects/SBG/Trait mapping/Global_trait_maps_Moreno_Martinez_2018_Version2_3km_resolution'
# lpj.path <- '~/Current Projects/SBG/LPJ/Reflectance_Data/version2/'
lpj.path <- Sys.getenv('lpjpath')
fp <- Sys.getenv('traitpath')


lma.name <- 'LDMC_3km_v1'
n.name <- 'LNC_3km_v1'
p.name <- 'LPC_3km_v1'
sla.name <- 'SLA_3km_v1'

lpj.nc <- 'lpj-prosail_levelC_DR_Version021_m_2020.nc'

# create PLSR data.frame --------------------------------------------------

lpj.array <- nc_open(file.path(lpj.path, lpj.nc)) %>% ncvar_get('DR') %>% aperm(c(2,1,3,4))

cells <- crds( rast(lpj.array[,,1,1], ext = c(-180,180,-90,90), crs = crs('EPSG:4326')) ) %>% vect(crs = crs('EPSG:4326'))

ldmc <- rast(file.path(trait.path,lma.name, paste0(lma.name,'.tif'))) %>% terra::extract(cells, ID = F)
lnc <- rast(file.path(trait.path,n.name, paste0(n.name,'.tif'))) %>% terra::extract(cells, ID = F)
lpc <- rast(file.path(trait.path,p.name, paste0(p.name,'.tif'))) %>% terra::extract(cells, ID = F)
sla <- rast(file.path(trait.path,sla.name, paste0(sla.name,'.tif'))) %>% terra::extract(cells, ID = F)

lpj.df <- as.data.frame(rast(lpj.array[,,,7]), xy = T)

plsr.df <- cbind.data.frame(lpj.df, ldmc, lnc, lpc, sla) %>% na.exclude()
names(plsr.df) <- c('x', 'y', paste0('wave',seq(400,2500,10)), lma.name, n.name, p.name, sla.name)
# run PLSR ----------------------------------------------------------------

# This is for July
print('LMA PLSR')
lma.coefs <- runPLSR(plsr.df, data.var = lma.name, train.size = 5000, plots = F,
                 jk.prop = 0.20)
print('N PLSR')
n.coefs <- runPLSR(plsr.df, data.var = n.name, train.size = 5000, plots = F,
                 jk.prop = 0.20)
print('P PLSR')
p.coefs <- runPLSR(plsr.df, data.var = p.name, train.size = 5000, plots = F,
                 jk.prop = 0.20)
print('SLA PLSR')
sla.coefs <- runPLSR(plsr.df, data.var = sla.name, train.size = 5000, plots = F,
                 jk.prop = 0.20)

# writing coeffs
print('Writing coefficients')
coeff.df <- data.frame(coeff = c('Intercept', seq(400,2500,10)), lma = lma.coefs, n = n.coefs, p = p.coefs, sla = sla.coefs)
write_csv(coeff.df, file.path(fp, 'LPJ_TRY_coeffs_July.csv'))




source('~/R/clean.R')
library(terra)
library(ncdf4)
library(tidyverse)
source('~/Current Projects/SBG/LPJ/Global_trait_PLSRs/R_scripts/lpj_plsr_functions.R')


# paths -------------------------------------------------------------------


trait.path <- '~/Current Projects/SBG/Trait mapping/Global_trait_maps_Moreno_Martinez_2018_Version2_3km_resolution'
lpj.path <- '~/Current Projects/SBG/LPJ/Reflectance_Data/version2/'



lma.name <- 'LDMC_3km_v1'
n.name <- 'LNC_3km_v1'
p.name <- 'LPC_3km_v1'
sla.name <- 'SLA_3km_v1'

lpj.nc <- 'lpj-prosail_levelC_DR_Version021_m_2020.nc'

# create PLSR data.frame --------------------------------------------------

lpj.array <- nc_open(file.path(lpj.path, lpj.nc)) %>% ncvar_get('DR') %>% aperm(c(2,1,3,4))

lpj.r <- rast(lpj.array[,,,7], crs = crs('EPSG:4326'), ext = c(-180,180,-90,90))
lpj.r[lpj.r>0.7] <- NA
lpj.r[lpj.r==0] <- NA

cells <- crds(lpj.r, na.rm = F) %>% vect(crs = crs('EPSG:4326'))

ldmc <- rast(file.path(trait.path, lma.name, paste0(lma.name,'.tif')))
ldmc[ldmc<0.1 | ldmc>0.5] <- NA
ldmc <- ldmc %>% terra::extract(cells, ID = F)

lnc <- rast(file.path(trait.path,n.name, paste0(n.name,'.tif')))
lnc[lnc<13 | lnc>25] <- NA
lnc <- lnc %>% terra::extract(cells, ID = F)

lpc <- rast(file.path(trait.path,p.name, paste0(p.name,'.tif'))) 
lpc[lpc<0.8 | lpc>2.1] <- NA
lpc <- lpc %>% terra::extract(cells, ID = F)

sla <- rast(file.path(trait.path,sla.name, paste0(sla.name,'.tif')))
sla[sla<7 | sla>21] <- NA
sla <- sla %>% terra::extract(cells, ID = F)

lpj.df <- as.data.frame(lpj.r, xy = T, na.rm = F)
plsr.data <- cbind.data.frame(lpj.df, ldmc, lnc, lpc, sla) %>% na.exclude()
names(plsr.data) <- c('x', 'y', paste0('wave',seq(400,2500,10)), lma.name, n.name, p.name, sla.name)


# run PLSR ----------------------------------------------------------------


# This is for July

# LMA
lma.coefs <- runPLSR(plsr.data, data.var = lma.name, band.prefix = 'wave', train.size = 2000, plots = F,
                 jk.test = F, jk.prop = 0.15, jk.iterations = 20, jk.comps = 5)
n.coefs <- runPLSR(plsr.data, data.var = n.name, band.prefix = 'wave', train.size = 2000, plots = F,
                     jk.test = F, jk.prop = 0.15, jk.iterations = 20, jk.comps = 5)
p.coefs <- runPLSR(plsr.data, data.var = p.name, band.prefix = 'wave', train.size = 2000, plots = F,
                   jk.test = F, jk.prop = 0.15, jk.iterations = 20, jk.comps = 5)
sla.coefs <- runPLSR(plsr.data, data.var = sla.name, band.prefix = 'wave', train.size = 2000, plots = F,
                   jk.test = F, jk.prop = 0.15, jk.iterations = 20, jk.comps = 5)
coeff.df <- data.frame(coeff = c('Intercept', seq(400,2500,10)), lma = lma.coefs, n = n.coefs, p = p.coefs, sla = sla.coefs)
write_csv(coeff.df, '~/Current Projects/SBG/LPJ/Global_trait_PLSRs/Example_LPJ-PROSAIL_PLSR_coefficients.csv')

lma.map <- trait.map(lpj.r, coeffs = coeff.df$lma[-1], intercept = coeff.df$lma[1], coeffs_wl = seq(400,2500,10))
n.map <- trait.map(lpj.r,  coeffs = coeff.df$n[-1], intercept = coeff.df$n[1], coeffs_wl = seq(400,2500,10))
p.map <- trait.map(lpj.r, coeff.df$p[-1], coeff.df$p[1], coeffs_wl = seq(400,2500,10))
sla.map <- trait.map(lpj.r, coeff.df$sla[-1], coeff.df$sla[1], coeffs_wl = seq(400,2500,10))


# plotting ----------------------------------------------------------------

library(ggplot2)
library(tidyterra)

ggplot() +
    geom_spatraster(data = lma.map) +
    scale_fill_gradientn(colors = c("wheat2", "darkgreen"), limits = c(0.25, 0.4), na.value = 'transparent') +
    theme_void() +
    labs(title = 'LPJ-PROSAIL estimated LMA (g/g)')

ggplot() +
    geom_spatraster(data = n.map) +
    scale_fill_gradientn(colors = c("wheat2", "darkgreen"), limits = c(15, 25), na.value = 'transparent') +
    theme_void() +
    labs(title = 'LPJ-PROSAIL estimated Leaf N (mg/g)')

ggplot() +
    geom_spatraster(data = p.map) +
    scale_fill_gradientn(colors = c("wheat2", "darkgreen"), limits = c(0.9, 1.6), na.value = 'transparent') +
    theme_void() +
    labs(title = 'LPJ-PROSAIL estimated Leaf P (mg/g)')

ggplot() +
    geom_spatraster(data = sla.map) +
    scale_fill_gradientn(colors = c("wheat2", "darkgreen"), limits = c(7, 20), na.value = 'transparent') +
    theme_void() +
    labs(title = 'LPJ-PROSAIL estimated SLA (mm2/mg)')


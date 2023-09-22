#!/bin/bash




export specdiv='/discover/nobackup/bcurrey/Git/specdiv/functions/'
export path='/discover/nobackup/bcurrey/Global_trait_PLSRs/'
export lpjpath='/discover/nobackup/projects/SBG-DO/bcurrey/global_run_simulations/lpj_version_021/ncdf_outputs/'
export lpjnc='lpj_prosail_Version021_DR_m_2022.nc'
export shpfile='discover/nobackup/bcurrey/Animation/Shapefile/Would_boundaries_noGreenland.shp'
export datapath="$path/data"

Rstudio "$path/R/LPJ_specdiv_monthly.R"

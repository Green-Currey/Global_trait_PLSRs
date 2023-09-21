#!/bin/bash

#......................................................................
#  E D I T:   S L U R M   A R R A Y   J O B    D E T A I L S
#......................................................................
#SBATCH --job-name=Traits
#SBATCH --time=01:59:00
#SBATCH --account=s3673

export lpjpath='/discover/nobackup/projects/SBG-DO/bcurrey/global_run_simulations/lpj_prosail_v21/ncdf_outputs/'
export datapath='/discover/nobackup/bcurrey/Global_trait_PLSRs/data/'
export lpjnc='lpj-prosail_levelC_DR_Version021_m_2020.nc'

Rscript '/discover/nobackup/bcurrey/Global_trait_PLSRs/R/LPJ_trait_mapping.R'




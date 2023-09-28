#!/bin/bash

#......................................................................
#  E D I T:   S L U R M   A R R A Y   J O B    D E T A I L S
#......................................................................
#SBATCH --job-name=Traits
#SBATCH --time=00:59:00
#SBATCH --account=s3673
#SBATCH --mail-user brycecurrey93@gmail.com
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

export lpjpath='/discover/nobackup/projects/SBG-DO/bcurrey/global_run_simulations/lpj_prosail_v21/ncdf_outputs/'
export datapath='/discover/nobackup/bcurrey/Global_trait_PLSRs/data/'
export lpjnc='lpj-prosail_levelC_DR_Version021_m_2022.nc'

Rscript '/discover/nobackup/bcurrey/Global_trait_PLSRs/R_scripts/LPJ_trait_mapping.R'




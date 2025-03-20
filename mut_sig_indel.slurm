#!/bin/bash
#SBATCH --account=miller_mbm31
#SBATCH -J NMF
#SBATCH -p priority
#SBATCH -t 1-00:00:00
#SBATCH --mem=16g
#SBATCH -c 20
#SBATCH -o NMF_%j.out


export PATH=$HOME/Tools/shapeit.v2.904.3.10.0-693.11.6.el7.x86_64/bin:$PATH
export PATH=/home/boj924/miniconda3/envs/scan2/bin:$PATH
export PATH=/n/app/python/3.10.11.conda/bin:$PATH


Rscript /home/boj924/AD_Tau_PTA/analysis_scripts/mutsig_indel.R



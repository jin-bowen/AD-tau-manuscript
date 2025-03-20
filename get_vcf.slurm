#!/bin/bash
#SBATCH --account=miller_mbm31
#SBATCH -J test
#SBATCH -p short
#SBATCH -t 05:00
#SBATCH --mem=16g
#SBATCH -c 1

export PATH=$HOME/Tools/shapeit.v2.904.3.10.0-693.11.6.el7.x86_64/bin:$PATH
export PATH=/home/boj924/miniconda3/envs/scan2/bin:$PATH
export PATH=/n/app/python/3.10.11.conda/bin:$PATH

outdir=/home/boj924/PrimTau/results
file=$1
#ls /n/scratch/groups/bwh/miller_mbm31/PrimTau/00*/call_mutations/*/scan2_object_update.rda | while read file
#do
id=`echo $file | awk -F"/" '{print $(NF-1)}'`
Rscript /home/boj924/PrimTau/analysis_scripts/tab2vcf.R  ${file} \
		/n/scratch/groups/bwh/miller_mbm31/PrimTau/0025/scan.yaml \
		${id} ${outdir}/${id}
#done

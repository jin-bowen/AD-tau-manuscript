#!/bin/bash
#SBATCH --account=miller_mbm31
#SBATCH -J suminfor
#SBATCH -p short
#SBATCH -t 1:00:00
#SBATCH --mem=16g
#SBATCH -c 1
#SBATCH -o scan2.%j

module load gcc/6.2.0
module load gatk

export PATH=$HOME/shapeit.v2.904.3.10.0-693.11.6.el7.x86_64/bin:$PATH
export PATH=/home/boj924/miniconda3/envs/scan2/bin:$PATH
export PATH=/n/app/python/3.10.11.conda/bin:$PATH
export PATH=/home/boj924/Tools/Eagle_v2.4.1/:$PATH

/home/boj924/Tools/annovar-20150322/annotate_variation.pl -buildver hg19 -downdb -webfrom annovar dbnsfp30a /n/data1/bwh/pathology/miller/lab/ref/humandb/


workdir=/n/data1/bwh/pathology/miller/lab/ctrl_hg19
Rscript /home/boj924/AD_Tau_PTA/analysis_scripts/extract_scan2_results.R \
	${workdir} \
	/home/boj924/AD_Tau_PTA/Control_config/scan2_summary.csv

workdir=/n/data1/bwh/pathology/miller/lab/AD_hg19
Rscript /home/boj924/AD_Tau_PTA/analysis_scripts/extract_scan2_results.R \
	${workdir} \
	/home/boj924/AD_Tau_PTA/AD_tau_config/scan2_meta_summary.csv

workdir=/n/data1/bwh/pathology/miller/lab/AD_Tau/results_hg19
Rscript /home/boj924/AD_Tau_PTA/analysis_scripts/extract_scan2_results.R \
	${workdir} \
	/home/boj924/AD_Tau_PTA/AD_tau_config/scan2_meta_summary.csv


mkdir -p ${workdir}/annovar
sed 1d ${workdir}/snv_list.combined.txt | awk '{print $3, $4, $4, $5, $6, $1, $2}' > ${workdir}/annovar/snv_list.combined.annovar 
/home/boj924/Tools/annovar-20150322/annotate_variation.pl -out ${workdir}/annovar/snv_list.annovar  \
			-build hg19 ${workdir}/annovar/snv_list.combined.annovar \
			/n/data1/bch/genetics/lee/shulin/softwares/annovar/humandb/


workdir=/n/data1/bwh/pathology/miller/lab/ctrl_hg19
workdir=/n/data1/bwh/pathology/miller/lab/AD_hg19
workdir=/n/data1/bwh/pathology/miller/lab/AD_Tau/results_hg19
/home/boj924/Tools/annovar-20150322/table_annovar.pl ${workdir}/annovar/snv_list.combined.annovar \
			/n/data1/bwh/pathology/miller/lab/ref/humandb/ \
			-protocol dbnsfp47a -operation f --otherinfo -build hg19 -nastring .			

cat $miller/AD_Tau/results_hg19/annovar/snv_list.combined.annovar.hg19_multianno.txt \
	$miller/ctrl_hg19/annovar/snv_list.combined.annovar.hg19_multianno.txt \
	$miller/AD_hg19/annovar/snv_list.combined.annovar.hg19_multianno.txt > ~/all_snv_list.annovar.hg19_multianno


workdir=/n/data1/bwh/pathology/miller/lab/ctrl_hg19
workdir=/n/data1/bwh/pathology/miller/lab/AD_hg19
workdir=/n/data1/bwh/pathology/miller/lab/AD_Tau/results_hg19
awk 'OFS=FS{print $3, $4, ".", $5,$6, ".","PASS", $1","$2}' ${workdir}/indel_list.combined.txt  | \
	sed 1d | sed 's/\s\+/\t/g'  > ${workdir}/annovar/indel_list.combined.vcf
sed -i '1 i#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO ' ${workdir}/annovar/indel_list.combined.vcf 
sed -i '1 i\\##fileformat=VCF' ${workdir}/annovar/indel_list.combined.vcf
sed -i 's/\s\+/\t/g' ${workdir}/annovar/indel_list.combined.vcf

/n/data1/bch/genetics/lee/shulin/softwares/annovar/convert2annovar.pl  -format vcf4 ${workdir}/annovar/indel_list.combined.vcf \
		-includeinfo > ${workdir}/annovar/indel_list.annovar

/home/boj924/Tools/annovar-20150322/annotate_variation.pl -out ${workdir}/annovar/indel_list.annovar  \
			-build hg19 ${workdir}/annovar/indel_list.annovar \
			/n/data1/bch/genetics/lee/shulin/softwares/annovar/humandb/

/home/boj924/Tools/annovar-20150322/table_annovar.pl ${workdir}/annovar/indel_list.annovar \
			/n/data1/bwh/pathology/miller/lab/ref/humandb  \
			-protocol dbnsfp47a --operation f --otherinfo  -buildver hg19 -nastring . 

cat $miller/AD_Tau/results_hg19/annovar/indel_list.annovar.hg19_multianno.txt \
	$miller/ctrl_hg19/annovar/indel_list.annovar.hg19_multianno.txt \
	$miller/AD_hg19/annovar/indel_list.annovar.hg19_multianno.txt > ~/all_indel_list.annovar.hg19_multianno


cat $miller/AD_Tau/results_hg19/annovar/indel_list.annovar.exonic_variant_function \
	$miller/ctrl_hg19/annovar/indel_list.annovar.exonic_variant_function \
	$miller/AD_hg19/annovar/indel_list.annovar.exonic_variant_function > ~/all_indel_list.annovar.exonic_variant_function


 
###t2t
#dir=/n/data1/bwh/pathology/miller/lab/AD_Tau_PTA/results_T2T
#outdir=/home/boj924/AD_Tau_PTA/results/t2t_AD_neurons
#Rscript /home/boj924/AD_Tau_PTA/get_indel_infor.R \
#	${dir}/2207/call_mutations/2207_noTau_A10/scan2_object.rda ${outdir}/2207_noTau_A10_t2t
#Rscript /home/boj924/AD_Tau_PTA/get_indel_infor.R \
#	${dir}/2207/call_mutations/2207_Tau_A5/scan2_object.rda ${outdir}/2207_Tau_A5_t2t
#Rscript /home/boj924/AD_Tau_PTA/get_indel_infor.R \
#	${dir}/2207/call_mutations/2207_noTau_A12/scan2_object.rda ${outdir}/2207_noTau_A12_t2t
#Rscript /home/boj924/AD_Tau_PTA/get_indel_infor.R \
#	${dir}/2207/call_mutations/2207_Tau_A2/scan2_object.rda ${outdir}/2207_Tau_A2_t2t
#
#
#Rscript /home/boj924/AD_Tau_PTA/get_indel_infor.R \
#	${dir}/2207_IndelRealign/call_mutations/2207_noTau_A10/scan2_object.rda ${dir}/2207_noTau_A10_t2t_IndelRealign
#Rscript /home/boj924/AD_Tau_PTA/get_indel_infor.R \
#	${dir}/2207_IndelRealign/call_mutations/2207_Tau_A5/scan2_object.rda ${dir}/2207_Tau_A5_t2t_IndelRealign
#Rscript /home/boj924/AD_Tau_PTA/get_indel_infor.R \
#	${dir}/2207_IndelRealign/call_mutations/2207_noTau_A12/scan2_object.rda ${dir}/2207_noTau_A12_t2t_IndelRealign
#Rscript /home/boj924/AD_Tau_PTA/get_indel_infor.R \
#	${dir}/2207_IndelRealign/call_mutations/2207_Tau_A2/scan2_object.rda ${dir}/2207_Tau_A2_t2t_IndelRealign
#
#
#gatk3 -T HaplotypeCaller --dontUseSoftClippedBases -l INFO --dbsnp /n/data1/bwh/pathology/miller/lab/ref/common_all_20180423.vcf.gz -mmq 60 -ERC BP_RESOLUTION --min_base_quality_score 20 -L  2:182816550-182816600  -R /n/data1/bwh/pathology/miller/lab/ref/hs37d5/hs37d5.fa -I /n/data1/bch/genetics/lee/lan/CTE/data/bam/1278_heart_bulk.bam -o 1278_bulk_short.gvcf.gz   -G Standard    -G AS_Standard
#
#gatk3 -T HaplotypeCaller --dontUseSoftClippedBases -l INFO --dbsnp /n/data1/bwh/pathology/miller/lab/ref/common_all_20180423.vcf.gz -mmq 60 -ERC BP_RESOLUTION --min_base_quality_score 20 -L 2:182816550-182816600 -R /n/data1/bwh/pathology/miller/lab/ref/hs37d5/hs37d5.fa -I /n/data1/bch/genetics/lee/lan/CTE/data/bam/1278BA9-A.b37.bam -o 1278BA9-A_short.gvcf.gz   -G Standard    -G AS_Standard
#
#gatk3 -T CombineGVCFs  -R $miller/ref/hs37d5/hs37d5.fa -V 1278BA9-A_short.gvcf.gz -V 1278_bulk_short.gvcf.gz -o combine2.vcf.gz -G StandardAnnotation -G AS_StandardAnnotation
#
#gatk3 -T GenotypeGVCFs  --dbsnp /n/data1/bwh/pathology/miller/lab/ref/common_all_20180423.vcf.gz  -R $miller/ref/hs37d5/hs37d5.fa  -V combine2.vcf.gz -o combine2_genotype.vcf.gz -G StandardAnnotation -G AS_StandardAnnotation
#
#####AD
yaml=/n/scratch/users/b/boj924/AD_hg19/1995/scan.yaml
meta=/home/boj924/AD_Tau_PTA/AD_case_config/scan2_panel_meta_AD.csv 
res_dir=/home/boj924/AD_Tau_PTA/results/AD_neurons
yaml_dir=/n/scratch/users/b/boj924/AD_hg19
out_dir=/home/boj924/AD_Tau_PTA/results/AD_neurons

######control
meta=/home/boj924/AD_Tau_PTA/AD_control_config/scan2.panel_meta.control.csv 
res_dir=/home/boj924/AD_Tau_PTA/results/control_neurons
yaml_dir=/n/scratch/users/b/boj924/ctrl_hg19
out_dir=/home/boj924/AD_Tau_PTA/results/control_neurons

######tau
meta=/home/boj924/AD_Tau_PTA/AD_tau_config/scan2_panel_meta_tau.csv
res_dir=/home/boj924/AD_Tau_PTA/results/tau_neurons
yaml_dir=/n/data1/bwh/pathology/miller/lab/AD_Tau_PTA/results_hg19
out_dir=/home/boj924/AD_Tau_PTA/results/tau_neurons
cat ${meta} | sed 1d | cut -d, -f1-2 | while read line
do 
  id=`echo $line | awk -F"," '{print $1}'`; cell=`echo $line | awk -F"," '{print $2}'`
  Rscript /home/boj924/AD_Tau_PTA/analysis_scripts/tab2vcf.R ${res_dir}/${cell}_snv.tab \
		${yaml} ${cell} ${out_dir}/${cell} snv
  bgzip ${out_dir}/${cell}_snv.vcf
done

cat ${meta} | sed 1d | cut -d, -f1-2 | while read line
do 
  id=`echo $line | awk -F"," '{print $1}'`; cell=`echo $line | awk -F"," '{print $2}'`
  Rscript /home/boj924/AD_Tau_PTA/analysis_scripts/tab2vcf.R ${res_dir}/${cell}_indel.tab \
		${yaml} ${cell} ${out_dir}/${cell} indel
  bgzip ${out_dir}/${cell}_indel.vcf
done


workdir=/n/data1/bwh/pathology/miller/lab/AD_Tau_PTA_results
cat /n/data1/bwh/pathology/miller/lab/AD_Tau/results_hg19/indel_list.combined.txt \
	/n/data1/bwh/pathology/miller/lab/AD_hg19/indel_list.combined.txt > /n/data1/bwh/pathology/miller/lab/AD_Tau_PTA_results/all_indel_list.combined.txt

awk 'OFS=FS{print $3, $4, ".", $5,$6, ".","PASS", $1","$2}' ${workdir}/all_indel_list.combined.txt  | \
	sed 1d | sed 's/\s\+/\t/g'  > ${workdir}/all_indel_list.combined.vcf
sed -i '1 i#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO ' ${workdir}/all_indel_list.combined.vcf 
sed -i '1 i\\##fileformat=VCF' ${workdir}/all_indel_list.combined.vcf
sed -i 's/\s\+/\t/g' ${workdir}/all_indel_list.combined.vcf

/n/data1/bch/genetics/lee/shulin/softwares/annovar/convert2annovar.pl  -format vcf4 ${workdir}/all_indel_list.combined.vcf \
		-includeinfo > ${workdir}/indel_list.annovar

/home/boj924/Tools/annovar-20150322/annotate_variation.pl -out ${workdir}/indel_list.annovar  \
			-build hg19 ${workdir}/indel_list.annovar \
			/n/data1/bwh/pathology/miller/lab/ref/humandb/


cat  /n/data1/bwh/pathology/miller/lab/AD_Tau/results_hg19/Permute_Tau/perms_by_cell_indel.combined.txt \
	/n/data1/bwh/pathology/miller/lab/AD_Tau/results_hg19/Permute_noTau/perms_by_cell_indel.combined.txt \
	/n/data1/bwh/pathology/miller/lab/AD_hg19/permute/perms_by_cell_indel.combined.txt \
		> /n/data1/bwh/pathology/miller/lab/AD_Tau_PTA_results/all_perms_by_cell_indel.combined.txt


awk 'OFS=FS{print "chr"$1, $2, ".", ".", ".", ".", "PASS", $7";"$5}' ${workdir}/all_perms_by_cell_indel.combined.txt | \
        sed 1d | sed 's/\s\+/\t/g'  > ${workdir}/all_perms_by_cell_indel.combined.vcf
sed -i '1 i#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO ' ${workdir}/all_perms_by_cell_indel.combined.vcf
sed -i '1 i\\##fileformat=VCF' ${workdir}/all_perms_by_cell_indel.combined.vcf
sed -i 's/\s\+/\t/g' ${workdir}/all_perms_by_cell_indel.combined.vcf

/n/data1/bch/genetics/lee/shulin/softwares/annovar/convert2annovar.pl -format vcf4 ${workdir}/all_perms_by_cell_indel.combined.vcf \
                -includeinfo > ${workdir}/all_perms_by_cell_indel.combined.annovar

/home/boj924/Tools/annovar-20150322/annotate_variation.pl -out ${workdir}/all_perms_by_cell_indel.combined  \
                        -build hg19 ${workdir}/all_perms_by_cell_indel.combined.annovar \
			/n/data1/bwh/pathology/miller/lab/ref/humandb/







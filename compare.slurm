#!/bin/bash
#SBATCH --account=miller_mbm31
#SBATCH -J panel
#SBATCH -p short
#SBATCH -t 1:00:00
#SBATCH --mem=8g
#SBATCH -c 1
#SBATCH -o summary_%j.out

dir=$1
output=all_res.tab
echo ",muttype,ncalls,callable.sens,callable.bp,callable.burden,rate.per.gb,burden,somatic.sens,pre.genotyping.burden,unsupported.filters" > ${output}
ls ${dir}/*.output | while read file
do 
  id=`basename $file | sed 's/\.output//g'`
  cat $file | sed 1d | sed "1s/^/${id},snv,/" | sed "2s/^/${id},indel,/" >> ${output}
done



cat ../AD_control_config/scan2.panel_meta.control.csv | sed 1d | grep -v "1278" | grep -v "bulk" | while read row
do
  id=`echo $row | awk -F"," '{print $2}'`
  sample=`echo $row | awk -F"," '{print $1}'`

  echo $id $sample

  id=$1
  sample=$2
  mut=$3

  comm -13 <(cut -d, -f1-2,4-5 control_neurons/${id}_${mut}.tab | sort )  \
	<( cut -d, -f1-2,4-5 lan_control_neurons/${id}_${mut}.tab | sort ) | \
	 sed 's/"//g' > lan_control_neurons/${id}_${mut}_unique.tab
  

  cat lan_control_neurons/${id}_${mut}_unique.tab | \
	cut -d, -f2 | while read pos
	do 
		zcat /n/scratch/users/b/boj924/ctrl_hg19/${sample}_gatk3_gvcf/call_mutations/integrated_table_update.tab.gz | \
		awk -F"\t" -v var="$pos" '$2==var{print}' >> control_neurons/${id}_${mut}_LanUnique.tab
	done

  line=`zcat /n/scratch/users/b/boj924/ctrl_hg19/${sample}_gatk3_gvcf/call_mutations/integrated_table_update.tab.gz | head -1`
  sed -i "1s/^/$line\n/" control_neurons/${id}_${mut}_LanUnique.tab
  python /home/boj924/AD_Tau_PTA/panel_construction/panel_test.py  control_neurons/${id}_${mut}_LanUnique.tab control_neurons/${id}_${mut}_LanUnique.test
  cat lan_control_neurons/${id}_${mut}_unique.tab | wc -l
  cat control_neurons/${id}_${mut}_LanUnique.tab | wc -l

  comm -23 <(cut -d, -f1-2,4-5 control_neurons/${id}_${mut}.tab | sort )        \
	 <( cut -d, -f1-2,4-5 lan_control_neurons/${id}_${mut}.tab | sort ) |    \
	sed 's/"//g' >  control_neurons/${id}_${mut}_unique.tab

  cat control_neurons/${id}_${mut}_unique.tab | \
	cut -d, -f2 | while read pos
	do 
		zcat /n/data1/bch/genetics/lee/lan/CTE/scan2_control/${sample}/call_mutations/integrated_table.tab.gz | \
		awk -F"\t" -v var="$pos" '$2==var{print}' >> lan_control_neurons/${id}_${mut}_CurrUnique.tab
	done
  line=`zcat /n/data1/bch/genetics/lee/lan/CTE/scan2_control/${sample}/call_mutations/integrated_table.tab.gz | head -1`
  sed -i "1s/^/$line\n/" lan_control_neurons/${id}_${mut}_CurrUnique.tab
  python /home/boj924/AD_Tau_PTA/panel_construction/panel_test.py  lan_control_neurons/${id}_${mut}_CurrUnique.tab lan_control_neurons/${id}_${mut}_CurrUnique.test
  cat control_neurons/${id}_${mut}_unique.tab | wc -l
  cat lan_control_neurons/${id}_${mut}_CurrUnique.tab | wc -l

done




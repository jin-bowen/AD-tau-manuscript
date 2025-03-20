

analysis_scripts=/home/boj924/AD_Tau_PTA/analysis_scripts
results=/home/boj924/AD_Tau_PTA/results/
meta=/home/boj924/AD_Tau_PTA/metafiles/

python ${analysis_scripts}/compare_lan.py \
	${meta}/control_metafile_lan \
	${results}/control_neurons/all_res.tab \
	${meta}/all_clinic \
	${results}/ctrlNeuN	

python ${analysis_scripts}/compare_joe.py \
	${meta}/control_metafile_joe \
	${results}/control_neurons/all_res.tab \
	${meta}/all_clinic \
	${results}/ctrlNeuN_joe	

python ${analysis_scripts}/compare_lan.py \
	${meta}/AD_metafile_lan \
	${results}/AD_neurons/all_res.tab \
	${meta}/all_clinic \
	${results}/ADNeuN	

python ${analysis_scripts}/compare_LiRA.py \
	${results}/control_neurons/all_res.tab \
	${meta}/control_metafile_2022nature \
	${results}/ctrlNeuN	

python ${analysis_scripts}/compare_LiRA.py \
	${results}/AD_neurons/all_res.tab \
	${meta}/ADnCtrl_metafile_august \
	${results}/ADNeuN

python ${analysis_scripts}/utility.py sensitivity \
        ${results}/AD_neurons/all_res.tab \
        ${results}/control_neurons/all_res.tab \
        ${results}/outlier

python ${analysis_scripts}/age_burden.py \
	${results}/tau_AD_ctrl_res.tab \
	${meta}/all_meta.tab \
	${meta}/all_clinic \
	${results}/tau_AD_ctrl_burden

python ${analysis_scripts}/sig_contribute.py \
	${results}/tauADCtrl_knowsig_contribute.tab \
	${results}/tau_AD_ctrl_res.tab \
	${meta}/all_meta.tab \
	${meta}/all_clinic \
	${results}/tau_AD_ctrl_burden

python ${analysis_scripts}/scatter.py \
	${results}/tau_AD_ctrl_res.tab \
	${meta}/all_meta.tab \
	${meta}/all_clinic \
	${results}/tau_AD_burden

python ${analysis_scripts}/scatter.py \
	${results}/indel_sig.tab \
	${results}/tau_AD_ctrl_res.tab \
	${meta}/all_meta.tab \
	${results}/indel_tau


python SigExtra.py \
	~/AD_Tau_PTA/results/SBS96/Suggested_Solution/COSMIC_SBS96_Decomposed_Solution/Activities/COSMIC_SBS96_Activities.txt \
	~/AD_Tau_PTA/results snv cosmic

python SigExtra.py \
	 ~/AD_Tau_PTA/results/SBS96/Suggested_Solution/SBS96_De-Novo_Solution/Activities/SBS96_De-Novo_Activities_refit.txt \
	~/AD_Tau_PTA/results snv denovo

python SigExtra.py \
	~/AD_Tau_PTA/results/ID83/Suggested_Solution/De_Novo_Solution/Activities/De_Novo_Activities.txt \
	~/AD_Tau_PTA/results indel denovo

python SigExtra.py \
	~/AD_Tau_PTA/results/ID83/Suggested_Solution/Decompose_Solution/Activities/Decompose_Solution_Activities.txt \
	~/AD_Tau_PTA/results indel cosmic






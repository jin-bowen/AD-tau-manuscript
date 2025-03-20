import matplotlib.pyplot as plt
import matplotlib.font_manager
plt.rcParams['font.size'] = 12
plt.rcParams['axes.linewidth'] = 2
plt.style.use('seaborn-v0_8-poster')
plt.rcParams['axes.spines.right'] = False
plt.rcParams['axes.spines.top'] = False
import seaborn as sns
import numpy as np
np.set_printoptions(precision=8)
from scipy import stats
import statsmodels.api as sm
import statsmodels.formula.api as smf
import pandas as pd
import sys
import re

def corr_2bp():

	snv_residue_input = '/home/boj924/AD_Tau_PTA/figures/excessive_snv-burden.csv'
	indel_residue_input = '/home/boj924/AD_Tau_PTA/figures/excessive_indel-burden.csv'
	indel_residue_input = '/home/boj924/AD_Tau_PTA/figures/excessive_indel-twobp_burden.csv'

	out_name    = '/home/boj924/AD_Tau_PTA/figures/excessive_2bp'
	snv_residue = pd.read_csv(snv_residue_input, header=0)	
	indel_residue = pd.read_csv(indel_residue_input, header=0)
	all_residue = pd.merge(snv_residue, indel_residue, on='sample', suffixes=('_snv',''))

	fig,ax = plt.subplots(constrained_layout=True,figsize=(4, 4), frameon = False)
	color_palette= { 'AD': 'tab:orange', 'ctrl': 'tab:blue', 'Tau':'tab:red', 'noTau':'tab:pink' }	
	color_map = [ color_palette[group] for group in all_residue['group']]
	plt.scatter(all_residue['burden_residue_snv'], all_residue['burden_residue'], 
			c=color_map, s=20,)
	sns.regplot(x='burden_residue_snv', y='burden_residue', data=all_residue, 
			scatter=False, color=".1") 
	ax.set_xlabel('sSNV excess burden')
	ax.set_ylabel('2bp deletion excess burden')	

	res = stats.spearmanr(all_residue['burden_residue_snv'], all_residue['burden_residue'])
	print(res)
	res = stats.pearsonr(all_residue['burden_residue_snv'], all_residue['burden_residue'])

	plt.savefig(out_name+'_corr.pdf', dpi=200,  bbox_inches='tight')
	plt.show()


def corr():

	snv_residue_input = '/home/boj924/AD_Tau_PTA/figures/excessive_snv-burden.csv'
	indel_residue_input = '/home/boj924/AD_Tau_PTA/figures/excessive_indel-burden.csv'
	out_name    = '/home/boj924/AD_Tau_PTA/figures/excessive'
	snv_residue = pd.read_csv(snv_residue_input, header=0)	
	indel_residue = pd.read_csv(indel_residue_input, header=0)
	all_residue = pd.merge(snv_residue, indel_residue, on='sample', suffixes=('_snv',''))
#	all_residue = all_residue.loc[all_residue['group']!='ctrl']	

	fig,ax = plt.subplots(constrained_layout=True,figsize=(4, 4), frameon = False)
	color_palette= { 'AD': 'tab:orange', 'ctrl': 'tab:blue', 'Tau':'tab:red', 'noTau':'tab:pink' }	
	color_map = [ color_palette[group] for group in all_residue['group']]
	plt.scatter(all_residue['burden_residue_snv'], all_residue['burden_residue'], 
			c=color_map, s=20,)
	sns.regplot(x='burden_residue_snv', y='burden_residue', data=all_residue, 
			scatter=False, color=".1") 
	ax.set_xlabel('sSNV excess burden')
	ax.set_ylabel('indel excess burden')	

#	res = stats.pearsonr(all_residue['burden_residue_snv'], all_residue['burden_residue'])
	res = stats.spearmanr(all_residue['burden_residue_snv'], all_residue['burden_residue'])
	print(res)
	plt.savefig(out_name+'_corr.pdf', dpi=200,  bbox_inches='tight')
	plt.show()


def excess_2bp():

	out_name    = '/home/boj924/AD_Tau_PTA/figures/excessive'
	muttype = 'indel'

	tbp_count = '/home/boj924/AD_Tau_PTA/results/indel_spectrum_count.tab'
	burden_in = '/home/boj924/AD_Tau_PTA/results/tau_AD_ctrl_res.tab'
	metafile_in = '/home/boj924/AD_Tau_PTA/metafiles/all_meta.tab'
	clinic_in   = '/home/boj924/AD_Tau_PTA/metafiles/all_clinic'
	AT8pct_in   = '/home/boj924/AD_Tau_PTA/metafiles/AT8.txt'
	outdir      = '/home/boj924/AD_Tau_PTA/figures/'

	burden = pd.read_csv(burden_in, sep=',', header=0,index_col=0)
	burden = burden.loc[burden['mutation_type']=='indel',]
	tbp_count_df = pd.read_csv(tbp_count, index_col=0, header=0)

	tbp_col = ['2bp_deletion_1','2bp_deletion_2', \
		'2bp_deletion_3','2bp_deletion_4','2bp_deletion_5','2bp_deletion_6+', \
		'2bp_deletion_with_microhomology_1']

	other_col = set(tbp_count_df.columns) - set(tbp_col)

	tbp_pct = tbp_count_df/tbp_count_df.sum(axis=1).values.reshape((-1,1))
	tbp_burden = pd.merge(tbp_pct, burden[['burden']], left_index=True, right_index=True)
	tbp_burden[tbp_pct.columns] = tbp_burden[tbp_pct.columns] * tbp_burden['burden'].values.reshape((-1,1))
	tbp_burden['twobp_burden'] = tbp_burden[tbp_col].sum(axis=1)
	tbp_burden['other_burden'] = tbp_burden[other_col].sum(axis=1)

	color_palette = 'Pastel1'
	metafile = pd.read_csv(metafile_in, header=0,dtype=str)
	clinic = pd.read_csv(clinic_in, header=0, sep='\t',dtype=str)
	metafile = pd.merge(metafile,clinic, left_on='donor',right_on='Case_ID')
	metafile['Age'] = metafile['Age'].astype(float)

	burden_meta = pd.merge(tbp_burden, metafile, left_index=True, right_on='sample')
	burden_meta = burden_meta.loc[burden_meta['sample']!='1995P_201001E3']
	ctrl = burden_meta.loc[burden_meta['group']=='ctrl']

	var = 'twobp_burden'
#	var = 'other_burden'
	md = smf.mixedlm(" %s ~ Age "%var, ctrl, groups=ctrl["donor"])
	mdf = md.fit()
	mlr_results = mdf.summary().tables[1].round(8)    	
	mlr_results['Coef.'] = mlr_results['Coef.'].astype(float)

	intercept = mlr_results.loc['Intercept','Coef.']
	slope = mlr_results.loc['Age','Coef.']

	plot_df = burden_meta
	color_palette= { 'AD': 'tab:orange', 'ctrl': 'tab:blue', 'Tau':'tab:red', 'noTau':'tab:pink' }	
	label = { 'AD': 'AD', 'ctrl': 'Control', 'Tau':'P-tau+', 'noTau':'P-tau-' }
	fig,ax = plt.subplots(constrained_layout=True,figsize=(5,4), frameon = False)

	plot_df['burden_residue'] = plot_df[['Age',var]].apply(lambda row:  row[var] - (row['Age']*slope+intercept), axis=1)	
	grp_plot_df = plot_df.sort_values(['Case_ID', 'group',var])

	color_list = [ color_palette[igrp] for igrp in plot_df['group'] ]	
	sns.boxplot(x='group',y='burden_residue',data=plot_df, ax=ax, palette=color_palette, \
		fill=False,saturation=1, fliersize=0.8, linewidth=2, legend=False)
	sns.stripplot(x='group',y='burden_residue',data = plot_df, ax=ax, palette=color_palette, edgecolor='k', size=8, linewidth=1, legend=False)

	ax.axhline(y=0, c='gray', ls='--', lw=1.5)
	ax.set_ylabel('Excess 2bp deletion burden')
	plt.legend(loc='center left',bbox_to_anchor=(1.01, 0.5), frameon=False)
	plt.savefig(out_name+'_'+muttype+'-'+var+'.pdf', dpi=200,  bbox_inches='tight')
	plt.show()
	
	select1 = burden_meta['group'] == 'Tau' 
	select2 = burden_meta['group'] == 'noTau'  
	res = stats.mannwhitneyu(burden_meta.loc[select1,'burden_residue'], burden_meta.loc[select2, 'burden_residue']) 
	print(res)

	select1 = burden_meta['group'] == 'Tau' 
	select2 = burden_meta['group'] == 'ctrl'  
	res = stats.mannwhitneyu(burden_meta.loc[select1,'burden_residue'], burden_meta.loc[select2, 'burden_residue']) 
	print(res)

	select1 = burden_meta['group'] == 'noTau' 
	select2 = burden_meta['group'] == 'ctrl'  
	res = stats.mannwhitneyu(burden_meta.loc[select1,'burden_residue'], burden_meta.loc[select2, 'burden_residue']) 
	print(res)

	select1 = burden_meta['group'] == 'AD'
	select2 = burden_meta['group'] == 'ctrl'  
	res = stats.mannwhitneyu(burden_meta.loc[select1,'burden_residue'], burden_meta.loc[select2, 'burden_residue']) 
	print(res)


def main_2bp():

	out_name    = '/home/boj924/AD_Tau_PTA/figures/excessive'
	muttype = 'indel'

	tbp_count = '/home/boj924/AD_Tau_PTA/results/indel_spectrum_count.tab'
	burden_in = '/home/boj924/AD_Tau_PTA/results/tau_AD_ctrl_res.tab'
	metafile_in = '/home/boj924/AD_Tau_PTA/metafiles/all_meta.tab'
	clinic_in   = '/home/boj924/AD_Tau_PTA/metafiles/all_clinic'
	AT8pct_in   = '/home/boj924/AD_Tau_PTA/metafiles/AT8.txt'
	outdir      = '/home/boj924/AD_Tau_PTA/figures/'

	burden = pd.read_csv(burden_in, sep=',', header=0,index_col=0)
	burden = burden.loc[burden['mutation_type']=='indel',]
	tbp_count_df = pd.read_csv(tbp_count, index_col=0, header=0)

	tbp_col = ['2bp_deletion_1','2bp_deletion_2', \
		'2bp_deletion_3','2bp_deletion_4','2bp_deletion_5','2bp_deletion_6+', \
		'2bp_deletion_with_microhomology_1']

	tbp_pct = tbp_count_df/tbp_count_df.sum(axis=1).values.reshape((-1,1))
	tbp_burden = pd.merge(tbp_pct, burden[['burden']], left_index=True, right_index=True)
	tbp_burden[tbp_pct.columns] = tbp_burden[tbp_pct.columns] * tbp_burden['burden'].values.reshape((-1,1))
	tbp_burden['twobp_burden'] = tbp_burden[tbp_col].sum(axis=1)

	color_palette = 'Pastel1'
	metafile = pd.read_csv(metafile_in, header=0,dtype=str)
	clinic = pd.read_csv(clinic_in, header=0, sep='\t',dtype=str)
	metafile = pd.merge(metafile,clinic, left_on='donor',right_on='Case_ID')
	metafile['Age'] = metafile['Age_(yrs)'].astype(float)

	burden_meta = pd.merge(tbp_burden, metafile, left_index=True, right_on='sample')
	burden_meta = burden_meta.loc[burden_meta['sample']!='1995P_201001E3']
	ctrl = burden_meta.loc[burden_meta['group']=='ctrl']

	var = 'twobp_burden'
	md = smf.mixedlm(" %s ~ Age "%var, ctrl, groups=ctrl["donor"])
	mdf = md.fit()
	mlr_results = mdf.summary().tables[1].round(8)    	
	mlr_results['Coef.'] = mlr_results['Coef.'].astype(float)

	intercept = mlr_results.loc['Intercept','Coef.']
	slope = mlr_results.loc['Age','Coef.']

	plot_df = burden_meta
	color_palette= { 'AD': 'tab:orange', 'ctrl': 'tab:blue', 'Tau':'tab:red', 'noTau':'tab:pink' }	
	label = { 'AD': 'AD', 'ctrl': 'Control', 'Tau':'P-tau+', 'noTau':'P-tau-' }
	fig,ax = plt.subplots(constrained_layout=True,figsize=(15,3), frameon = False)

	plot_df['burden_residue'] = plot_df[['Age',var]].apply(lambda row:  row[var] - (row['Age']*slope+intercept), axis=1)	
	grp_plot_df = plot_df.sort_values(['Case_ID', 'group', var])
	grp_plot_df.to_csv(out_name+'_'+muttype+'-'+var+'.csv')
	color_list = [ color_palette[igrp] for igrp in plot_df['group'] ]	

	sns.barplot(x='sample',y='burden_residue', hue='group',order=grp_plot_df['sample'], palette=color_palette, data=grp_plot_df, ax=ax, legend=False)
	plt.xticks(rotation=90, fontsize=8)

	ax.set_ylabel('Excess 2bp deletion burden')
	plt.legend(loc='center left',bbox_to_anchor=(1.01, 0.5), frameon=False)
	plt.savefig(out_name+'_'+muttype+'-'+var+'.pdf', dpi=200,  bbox_inches='tight')
	plt.show()


def excess():

	burden_in   = '/home/boj924/AD_Tau_PTA/results/tau_AD_ctrl_res.tab'
	sig_in      = '/home/boj924/AD_Tau_PTA/results/tauADCtrl_aging_contribute.tab'
	metafile_in = '/home/boj924/AD_Tau_PTA/metafiles/all_meta.tab'
	clinic_in   = '/home/boj924/AD_Tau_PTA/metafiles/all_clinic'
	out_name    = '/home/boj924/AD_Tau_PTA/figures/excessive'

	muttype     = 'snv'
	muttype     = 'indel'
	var         = 'burden'
#	var         = 'scaled_Signature_A'	
#	var         = 'scaled_Signature_C'	
		
	color_palette = 'Pastel1'
	burden = pd.read_csv(burden_in, sep=',', header=0,index_col=0)
	burden = burden.loc[burden['mutation_type']==muttype]
	burden.dropna(inplace=True)
	burden['sample'] = burden.index

	sig_tab = pd.read_csv(sig_in, sep=',', header=0, index_col=0)
	sig_tab = sig_tab.astype(float)
	sig_tab = sig_tab/sig_tab.sum(axis=1).values.reshape((-1,1))
	burden_sig = pd.merge(burden, sig_tab, left_index=True, right_index=True)

	cols=['Signature_A','Signature_C']
	scaled_cols=['scaled_' + x for x in cols]
	for i, item in enumerate(cols):
		burden_sig[scaled_cols[i]] = burden_sig.apply(lambda row: row['burden'] * row[item], axis=1 )

	metafile = pd.read_csv(metafile_in, header=0,dtype=str)
	clinic = pd.read_csv(clinic_in, header=0, sep='\t',dtype=str)
	metafile = pd.merge(metafile,clinic, left_on='donor',right_on='Case_ID')
	metafile['Age'] = metafile['Age'].astype(float)
	burden_meta = pd.merge(burden_sig,metafile, on='sample', suffixes=('','_meta'))
	burden_meta = burden_meta.loc[burden_meta['sample']!='1995P_201001E3']

	ctrl = burden_meta.loc[burden_meta['group']=='ctrl']

	md = smf.mixedlm(" %s ~ Age "%var, ctrl, groups=ctrl["donor"])
	mdf = md.fit()
	mlr_results = mdf.summary().tables[1].round(8)    	
	mlr_results['Coef.'] = mlr_results['Coef.'].astype(float)

	intercept = mlr_results.loc['Intercept','Coef.']
	slope = mlr_results.loc['Age','Coef.']

	plot_df = burden_meta
	color_palette= { 'AD': 'tab:orange', 'ctrl': 'tab:blue', 'Tau':'tab:red', 'noTau':'tab:pink' }	
	label = { 'AD': 'AD', 'ctrl': 'Control', 'Tau':'P-tau+', 'noTau':'P-tau-' }
	fig,ax = plt.subplots(constrained_layout=True,figsize=(5,4), frameon = False)

	plot_df['burden_residue'] = plot_df[['Age',var]].apply(lambda row:  row[var] - (row['Age']*slope+intercept), axis=1)	
	grp_plot_df = plot_df.sort_values(['Case_ID', 'group','burden'])

	color_list = [ color_palette[igrp] for igrp in plot_df['group'] ]	
	sns.boxplot(x='group',y='burden_residue',data=plot_df, ax=ax, palette=color_palette, \
			fill=False,saturation=1, fliersize=0.8, linewidth=2, legend=False)
	sns.stripplot(x='group',y='burden_residue',data = plot_df, ax=ax, palette=color_palette, edgecolor='k', size=8, linewidth=1,legend=False)

	ax.axhline(y=0, c='tab:blue', ls='--', lw=1)
	ax.set_ylabel('Excess burden')
	plt.legend(loc='center left',bbox_to_anchor=(1.01, 0.5), frameon=False)
	plt.savefig(out_name+'_'+muttype+'-'+var+'.pdf', dpi=200,  bbox_inches='tight')
	plt.show()
	
	select1 = burden_meta['group'] == 'Tau' 
	select2 = burden_meta['group'] == 'noTau'  
	res = stats.mannwhitneyu(burden_meta.loc[select1,'burden_residue'], burden_meta.loc[select2, 'burden_residue']) 
	print(res)

	select1 = burden_meta['group'] == 'Tau' 
	select2 = burden_meta['group'] == 'ctrl'  
	res = stats.mannwhitneyu(burden_meta.loc[select1,'burden_residue'], burden_meta.loc[select2, 'burden_residue']) 
	print(res)

	select1 = burden_meta['group'] == 'noTau' 
	select2 = burden_meta['group'] == 'ctrl'  
	res = stats.mannwhitneyu(burden_meta.loc[select1,'burden_residue'], burden_meta.loc[select2, 'burden_residue']) 
	print(res)

	select1 = burden_meta['group'] == 'AD'
	select2 = burden_meta['group'] == 'ctrl'  
	res = stats.mannwhitneyu(burden_meta.loc[select1,'burden_residue'], burden_meta.loc[select2, 'burden_residue']) 

	print(np.mean(burden_meta.loc[select1,'burden_residue']) - np.mean(burden_meta.loc[select2, 'burden_residue'])) 
	print(np.median(burden_meta.loc[select1,'burden_residue']) - np.median(burden_meta.loc[select2, 'burden_residue'])) 

	print(res)


def main():

	burden_in   = '/home/boj924/AD_Tau_PTA/results/tau_AD_ctrl_res.tab'
	metafile_in = '/home/boj924/AD_Tau_PTA/metafiles/all_meta.tab'
	clinic_in   = '/home/boj924/AD_Tau_PTA/metafiles/all_clinic'
	out_name    = '/home/boj924/AD_Tau_PTA/figures/excessive'
#	muttype = 'indel'
	muttype = 'snv'
	
	color_palette = 'Pastel1'
	burden = pd.read_csv(burden_in, sep=',', header=0,index_col=0)
	burden = burden.loc[burden['mutation_type']==muttype]
	burden.dropna(inplace=True)
	burden['sample'] = burden.index

	metafile = pd.read_csv(metafile_in, header=0,dtype=str)
	clinic = pd.read_csv(clinic_in, header=0, sep='\t',dtype=str)
	metafile = pd.merge(metafile,clinic, left_on='donor',right_on='Case_ID')
	metafile['Age'] = metafile['Age'].astype(float)
	burden_meta = pd.merge(burden,metafile, on='sample', suffixes=('','_meta'))
	burden_meta = burden_meta.loc[burden_meta['sample']!='1995P_201001E3']

	ctrl = burden_meta.loc[burden_meta['group']=='ctrl']

	var='burden'
	md = smf.mixedlm(" %s ~ Age "%var, ctrl, groups=ctrl["donor"])
	mdf = md.fit()
	mlr_results = mdf.summary().tables[1].round(8)    	
	mlr_results['Coef.'] = mlr_results['Coef.'].astype(float)

	intercept = mlr_results.loc['Intercept','Coef.']
	slope = mlr_results.loc['Age','Coef.']

	plot_df = burden_meta
	color_palette= { 'AD': 'tab:orange', 'ctrl': 'tab:blue', 'Tau':'tab:red', 'noTau':'tab:pink' }	
	label = { 'AD': 'AD', 'ctrl': 'Control', 'Tau':'P-tau+', 'noTau':'P-tau-' }
	fig,ax = plt.subplots(constrained_layout=True,figsize=(15,3), frameon = False)

	plot_df['burden_residue'] = plot_df[['Age','burden']].apply(lambda row:  row['burden'] - (row['Age']*slope+intercept), axis=1)	
	grp_plot_df = plot_df.sort_values(['Case_ID', 'group','burden'])
	grp_plot_df.to_csv(out_name+'_'+muttype+'-'+var+'.csv')

	color_list = [ color_palette[igrp] for igrp in plot_df['group'] ]	

	sns.barplot(x='sample',y='burden_residue', hue='group',order=grp_plot_df['sample'], palette=color_palette, data=grp_plot_df, ax=ax, legend=False)
	plt.xticks(rotation=90, fontsize=8)
	ax.set_yscale('log')

	ax.set_ylabel('Excess burden')
	plt.legend(loc='center left',bbox_to_anchor=(1.01, 0.5), frameon=False)
	plt.savefig(out_name+'_'+muttype+'-'+var+'.pdf', dpi=200,  bbox_inches='tight')
	

if __name__ == "__main__":
#	main()
#	corr()
	excess()
#	excess_2bp()
#	main_2bp()
#	corr_2bp()






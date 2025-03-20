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

def excess():

	burden_in   = '/home/boj924/AD_Tau_PTA/results/tau_AD_ctrl_res.tab'
	sig_in      = '/home/boj924/AD_Tau_PTA/results/tauADCtrl_aging_contribute.tab'
	metafile_in = '/home/boj924/AD_Tau_PTA/metafiles/all_meta.tab'
	clinic_in   = '/home/boj924/AD_Tau_PTA/metafiles/all_clinic'
	out_name    = '/home/boj924/AD_Tau_PTA/figures/excessive_antibody'

	muttype     = 'snv'
#	muttype     = 'indel'
	var         = 'burden'
#	var         = 'scaled_Signature_A'	
		
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
	burden_meta = burden_meta.loc[burden_meta['group'].isin(['ctrl','Tau'])]
	burden_meta['antibody'] = burden_meta['sample'].apply(lambda x: None if len(x.split('_')) < 3 else x.split('_')[1].lower())
	burden_meta['antibody'].fillna('ctrl', inplace=True)

	ctrl = burden_meta.loc[burden_meta['group']=='ctrl']

	md = smf.mixedlm(" %s ~ Age "%var, ctrl, groups=ctrl["donor"])
	mdf = md.fit()
	mlr_results = mdf.summary().tables[1].round(8)    	
	mlr_results['Coef.'] = mlr_results['Coef.'].astype(float)

	intercept = mlr_results.loc['Intercept','Coef.']
	slope = mlr_results.loc['Age','Coef.']

	plot_df = burden_meta
	color_palette= { 'tau': 'brown', 'ctrl': 'tab:blue', 't205':'magenta', 'ht7':'pink' }
	label = { 'tau': 'Ser404', 'ctrl': 'Control', 't205':'T205', 'ht7':'HT7' }
	fig,ax = plt.subplots(constrained_layout=True,figsize=(5,4), frameon = False)

	plot_df['burden_residue'] = plot_df[['Age',var]].apply(lambda row:  row[var] - (row['Age']*slope+intercept), axis=1)	
	grp_plot_df = plot_df.sort_values(['Case_ID', 'group','burden'])
	antibody_order = ['ctrl','tau','t205','ht7']

	color_list = [ color_palette[igrp] for igrp in plot_df['antibody'] ]	
	sns.boxplot(x='antibody',y='burden_residue',data=plot_df, ax=ax, palette=color_palette, order=antibody_order,
			fill=False,saturation=1, fliersize=0.8, linewidth=2, legend=False)
	sns.stripplot(x='antibody',y='burden_residue',data = plot_df, ax=ax, palette=color_palette, 
			order=antibody_order,edgecolor='k', size=8, linewidth=1,legend=False)

	ax.axhline(y=0, c='tab:blue', ls='--', lw=1)
	ax.set_ylabel('Excess burden')
	plt.legend(loc='center left',bbox_to_anchor=(1.01, 0.5), frameon=False)
	plt.savefig(out_name+'_'+muttype+'-'+var+'.pdf', dpi=200,  bbox_inches='tight')
	
	select1 = burden_meta['antibody'] == 'tau' 
	select2 = burden_meta['antibody'] == 'ctrl'  
	res = stats.mannwhitneyu(burden_meta.loc[select1,'burden_residue'], burden_meta.loc[select2, 'burden_residue']) 
	print(res)

	select1 = burden_meta['antibody'] == 't205' 
	select2 = burden_meta['antibody'] == 'ctrl'  
	res = stats.mannwhitneyu(burden_meta.loc[select1,'burden_residue'], burden_meta.loc[select2, 'burden_residue']) 
	print(res)

	select1 = burden_meta['antibody'] == 'ht7' 
	select2 = burden_meta['antibody'] == 'ctrl'  
	res = stats.mannwhitneyu(burden_meta.loc[select1,'burden_residue'], burden_meta.loc[select2, 'burden_residue']) 
	print(res)

	select1 = burden_meta['antibody'] == 'tau'
	select2 = burden_meta['antibody'] == 't205'  
	res = stats.mannwhitneyu(burden_meta.loc[select1,'burden_residue'], burden_meta.loc[select2, 'burden_residue']) 
	print(res)

	select1 = burden_meta['antibody'] == 'tau'
	select2 = burden_meta['antibody'] == 'ht7'  
	res = stats.mannwhitneyu(burden_meta.loc[select1,'burden_residue'], burden_meta.loc[select2, 'burden_residue']) 
	print(res)

def main():

	burden_in   = '/home/boj924/AD_Tau_PTA/results/tau_AD_ctrl_res.tab'
	metafile_in = '/home/boj924/AD_Tau_PTA/metafiles/all_meta.tab'
	clinic_in   = '/home/boj924/AD_Tau_PTA/metafiles/all_clinic'
	out_name    = '/home/boj924/AD_Tau_PTA/figures/excessive'
	muttype = 'indel'
	
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
	plt.show()
	

if __name__ == "__main__":
#	main()
#	corr()
	excess()






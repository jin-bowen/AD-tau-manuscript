import matplotlib.pyplot as plt
import matplotlib.font_manager
plt.rcParams['font.size'] = 10
plt.rcParams['axes.linewidth'] = 1
plt.style.use('seaborn-v0_8-poster')
plt.rcParams['axes.spines.right'] = False
plt.rcParams['axes.spines.top'] = False
import numpy as np
np.set_printoptions(precision=8)
import seaborn as sns
from scipy import stats
import statsmodels.api as sm
import statsmodels.formula.api as smf
import pandas as pd
import sys
import re

def print(*args):
	__builtins__.print(*("%.8f" % a if isinstance(a, float) else a for a in args))


def excess():

	burden_in   = '/home/boj924/AD_Tau_PTA/results/tau_AD_ctrl_res.tab'
	sig_in      = '/home/boj924/AD_Tau_PTA/results/tauADCtrl_aging_contribute.tab'
	metafile_in = '/home/boj924/AD_Tau_PTA/metafiles/all_meta.tab'
	clinic_in   = '/home/boj924/AD_Tau_PTA/metafiles/all_clinic'
	out_name     = '/home/boj924/AD_Tau_PTA/figures/antibody'

	muttype     = 'indel'
	muttype     = 'snv'
	var         = 'burden'	
		
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
	burden_meta = pd.merge(burden_sig,metafile, on='sample', suffixes=('','_meta'))
	burden_meta = burden_meta.loc[burden_meta['group'].isin(['ctrl','Tau'])]
	burden_meta['antibody'] = burden_meta['sample'].apply(lambda x: None if len(x.split('_')) < 3 else x.split('_')[1].lower())
	burden_meta['antibody'].fillna('ctrl', inplace=True)
	burden_meta['Age'] = burden_meta['Age'].astype(float)
#	burden_meta['antibody_group'] = burden_meta['antibody'].apply(lambda x: 'S404' if x=='tau' else( 'ctrl' if x=='ctrl' else 'other'))
	ctrl = burden_meta.loc[burden_meta['group']=='ctrl']

	md = smf.mixedlm(" %s ~ Age "%var, ctrl, groups=ctrl["donor"])
	mdf = md.fit()
	mlr_results = mdf.summary().tables[1].round(8)
	mlr_results['Coef.'] = mlr_results['Coef.'].astype(float)

	intercept = mlr_results.loc['Intercept','Coef.']
	slope = mlr_results.loc['Age','Coef.']

	plot_df = burden_meta
	color_palette= { 'tau': 'pink', 'ctrl': 'tab:blue', 't205':'magenta', 'ht7':'brown' }
	fig,ax = plt.subplots(constrained_layout=True,figsize=(5,4), frameon = False)

	plot_df['burden_residue'] = plot_df[['Age',var]].apply(lambda row:  row[var] - (row['Age']*slope+intercept), axis=1)
	color_list = [ color_palette[igrp] for igrp in plot_df['antibody'] ]
	sns.boxplot(x='antibody',y='burden_residue',data=plot_df, ax=ax, palette=color_palette, \
		fill=False,saturation=1, fliersize=0.8, linewidth=2, legend=False)
	sns.stripplot(x='antibody',y='burden_residue',data = plot_df, ax=ax, palette=color_palette, edgecolor='k', size=8, linewidth=1,legend=False)

	ax.axhline(y=0, c='tab:blue', ls='--', lw=1)
	ax.set_ylabel('Excess burden')
	plt.legend(loc='center left',bbox_to_anchor=(1.01, 0.5), frameon=False)
	plt.savefig(out_name+'_'+muttype+'.pdf', dpi=200,  bbox_inches='tight')
	plt.show()

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

def main():

	burden_in   = '/home/boj924/AD_Tau_PTA/results/tau_AD_ctrl_res.tab'
	sig_in      = '/home/boj924/AD_Tau_PTA/results/tauADCtrl_aging_contribute.tab'
	metafile_in = '/home/boj924/AD_Tau_PTA/metafiles/all_meta.tab'
	clinic_in   = '/home/boj924/AD_Tau_PTA/metafiles/all_clinic'
	outdir      = '/home/boj924/AD_Tau_PTA/figures/'

	muttype     = 'indel'
#	muttype     = 'snv'
	var         = 'burden'	
		
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
	burden_meta = pd.merge(burden_sig,metafile, on='sample', suffixes=('','_meta'))
	burden_meta = burden_meta.loc[burden_meta['group'].isin(['ctrl','Tau'])]
	burden_meta['antibody'] = burden_meta['sample'].apply(lambda x: None if len(x.split('_')) < 3 else x.split('_')[1].lower())
	burden_meta['antibody'].fillna('ctrl', inplace=True)
	burden_meta['Age'] = burden_meta['Age'].astype(float)

	out_name = outdir + 'ctrl_antibody_%s'%var
	md = smf.mixedlm(" %s ~ Age +  C(antibody, Treatment('ctrl')) "%var, burden_meta, groups=burden_meta["donor"])
	mdf = md.fit()
	mlr_results = mdf.summary().tables[1].round(8)
	mlr_results['Coef.'] = mlr_results['Coef.'].astype(float)
	print(mlr_results)
	print("\n")
	print(mdf.pvalues)

	md = smf.mixedlm(" %s ~ Age *  C(antibody, Treatment('ctrl')) "%var, burden_meta, groups=burden_meta["donor"])
	mdf = md.fit()
	mlr_results = mdf.summary().tables[1].round(8)
	mlr_results['Coef.'] = mlr_results['Coef.'].astype(float)

	plot_df = burden_meta
	color_palette= { 'tau': 'pink', 'ctrl': 'tab:blue', 't205':'magenta', 'ht7':'brown' }	
	label = { 'tau': 'Ser404', 'ctrl': 'Control', 't205':'T205', 'ht7':'HT7' }

	fig,ax = plt.subplots(constrained_layout=True,figsize=(6,4), frameon = False)
	for igrp, df_sub in plot_df.groupby('antibody',sort=False):

		xpos=np.arange(57, 91)
		base_intercept = mlr_results.loc['Intercept','Coef.']
		base_slope = mlr_results.loc['Age','Coef.']
		if igrp == 'ctrl':
			intercept = base_intercept
			slope = base_slope
			xpos=np.arange(0, 110)
		else:
			r = re.compile(".*\.%s.*"%igrp)
			sub_status = list(filter(r.match,mlr_results.index.to_list()))[0]
			add_on_intercept = mlr_results.loc[sub_status,'Coef.']
			add_on_slope = 0

			sub_status2 = list(filter(r.match,mlr_results.index.to_list()))[1]
			add_on_slope = mlr_results.loc[sub_status2,'Coef.'] 

			intercept=base_intercept + add_on_intercept
			slope = base_slope + add_on_slope
		ax.plot(xpos, slope * xpos + intercept, c=color_palette[igrp], lw=2, alpha=0.8)
		ax.scatter(x=df_sub['Age'], y=df_sub[var], c=color_palette[igrp], edgecolors='black', \
					label=label[igrp], linewidth=1,s=80)
		ax.set_xlabel('Age')
		ax.set_ylabel('%s per neuron'%muttype)
	plt.legend(loc='center left',bbox_to_anchor=(1.0, 0.5), frameon=False)
	plt.show()
#	plt.savefig(out_name+'_'+muttype+'-'+var+'.pdf', dpi=200,  bbox_inches='tight')

if __name__ == "__main__":
#	main()
	excess()



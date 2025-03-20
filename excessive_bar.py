import matplotlib.pyplot as plt
import matplotlib.font_manager
plt.rcParams['font.size'] = 10
plt.rcParams['axes.linewidth'] = 2
#plt.style.use('seaborn-v0_8-poster')
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

def print(*args):
	__builtins__.print(*("%.8f" % a if isinstance(a, float) else a for a in args))

def main():

	burden_in   = '/home/boj924/AD_Tau_PTA/results/tau_AD_ctrl_res.tab'
	metafile_in = '/home/boj924/AD_Tau_PTA/metafiles/all_meta.tab'
	clinic_in   = '/home/boj924/AD_Tau_PTA/metafiles/all_clinic'
	out_name    = '/home/boj924/AD_Tau_PTA/figures/excessive_bar'
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
	metafile['Age'] = metafile['Age_(yrs)'].astype(float)
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

	plot_df = burden_meta.loc[burden_meta['group'].isin(['Tau','noTau','AD'])]
	color_palette= { 'AD': 'tab:orange', 'ctrl': 'tab:blue', 'Tau':'tab:red', 'noTau':'tab:pink' }	
	label = { 'AD': 'AD', 'ctrl': 'Control', 'Tau':'P-tau+', 'noTau':'P-tau-' }
	fig,ax = plt.subplots(constrained_layout=True,figsize=(8,3), frameon = False)

	plot_df['burden_residue'] = plot_df[['Age','burden']].apply(lambda row:  row['burden'] - (row['Age']*slope+intercept), axis=1)	
	grp_plot_df = plot_df.sort_values(['Case_ID', 'group','burden'])

	color_list = [ color_palette[igrp] for igrp in plot_df['group'] ]	

	sns.stripplot(x='Case_ID',y='burden_residue', hue='group', dodge=True,linewidth=1,  
			palette=color_palette, data=grp_plot_df, ax=ax, legend=False, s=4)
	sns.boxplot(x='Case_ID',y='burden_residue', hue='group', palette=color_palette, 
			linewidth=1, data=grp_plot_df, ax=ax, legend=False, fliersize=0.1, )

	plt.xticks(rotation=90, fontsize=10)
	
	ax.set_xlabel('')
	ax.set_ylabel('Excess burden')
	ax.axhline(y=0,ls='--', lw=1)
	plt.legend(loc='center left',bbox_to_anchor=(1.01, 0.5), frameon=False)
	plt.savefig(out_name+'_'+muttype+'-'+var+'.pdf', dpi=200,  bbox_inches='tight')
	plt.show()
	
	df_reform = plot_df.pivot(index=['sample','donor'], columns='group', values='burden_residue').reset_index()
	print(df_reform)
	#for igrp, grp in df_reform.groupby('donor'):

if __name__ == "__main__":
	main()






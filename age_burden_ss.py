import matplotlib.pyplot as plt
import matplotlib.font_manager
plt.rcParams['font.size'] = 10
plt.rcParams['axes.linewidth'] = 1
plt.style.use('seaborn-v0_8-poster')
plt.rcParams['axes.spines.right'] = False
plt.rcParams['axes.spines.top'] = False
import numpy as np
np.set_printoptions(precision=8)
from scipy import stats
import seaborn as sns
import statsmodels.api as sm
import statsmodels.formula.api as smf
import pandas as pd
import sys
import re


def percent():

	metafile_in = '/home/boj924/Tn5-duplex-calling/tlpk_AD/all_meta.tab'
	clinic_in   = '/home/boj924/Tn5-duplex-calling/tlpk_AD/all_clinic'
	metafile = pd.read_csv(metafile_in, header=0,dtype=str)
	clinic = pd.read_csv(clinic_in, header=0,dtype=str)
	metafile = pd.merge(metafile,clinic, left_on='donor',right_on='Case_ID')
	metafile['Age'] = metafile['Age'].astype(float)

	ss_2bp = pd.read_csv('results/metacs_ss_2bp.csv', index_col=0, header=0)
	ds_2bp = pd.read_csv('results/metacs_ds_2bp_percent.csv', index_col=0, header=0)
	ss_2bp.dropna(inplace=True)
	ds_2bp.dropna(inplace=True)

	print(ss_2bp.head())
	print(metafile.head())
	ss_2bp = pd.merge(ss_2bp, metafile, left_index=True, right_on='sample')
	ds_2bp = pd.merge(ds_2bp, metafile, left_index=True, right_on='sample')

	ss_2bp['type']='ss'
	ds_2bp['type']='ds'
	all_2bp = pd.concat([ss_2bp, ds_2bp])

	color_palette= { 'AD': 'tab:orange', 'ctrl': 'tab:blue', 'Tau':'tab:red', 'noTau':'tab:pink' }	
	label = { 'AD': 'AD', 'ctrl': 'Control', 'Tau':'P-tau+', 'noTau':'P-tau-' }
	fig,ax = plt.subplots(constrained_layout=True,figsize=(5,5), frameon = False)

	order = ['ctrl','AD','Tau','noTau']

	sns.boxplot(x='group',y='percent', data=ss_2bp, ax=ax, palette=color_palette, \
		fill=False,saturation=1, fliersize=0.8, linewidth=2, legend=False, order=order)
	sns.stripplot(x='group',y='percent',data = ss_2bp, ax=ax, palette=color_palette, order=order, \
			edgecolor='k', size=8, linewidth=1,legend=False)
	ax.set_ylim([-0.1, 1.1])
	plt.savefig('figures/2bp_ss.pdf', dpi=200,  bbox_inches='tight')
	plt.show()

	select1 = ss_2bp['group'] == 'AD'
	select2 = ss_2bp['group'] == 'ctrl'
	res = stats.mannwhitneyu(ss_2bp.loc[select1,'percent'], ss_2bp.loc[select2, 'percent'])
	print(res)

	select1 = ss_2bp['group'] == 'Tau'
	select2 = ss_2bp['group'] == 'ctrl'
	res = stats.mannwhitneyu(ss_2bp.loc[select1,'percent'], ss_2bp.loc[select2, 'percent'])
	print(res)

	select1 = ss_2bp['group'] == 'noTau'
	select2 = ss_2bp['group'] == 'ctrl'
	res = stats.mannwhitneyu(ss_2bp.loc[select1,'percent'], ss_2bp.loc[select2, 'percent'])
	print(res)

	select1 = ss_2bp['group'] == 'Tau'
	select2 = ss_2bp['group'] == 'noTau'
	res = stats.mannwhitneyu(ss_2bp.loc[select1,'percent'], ss_2bp.loc[select2, 'percent'])
	print(res)

	fig,ax = plt.subplots(constrained_layout=True,figsize=(5,5), frameon = False)
	sns.boxplot(x='group',y='percent', data=ds_2bp, ax=ax, palette=color_palette, order=order,\
		fill=False,saturation=1, fliersize=0.8, linewidth=2, legend=False)
	sns.stripplot(x='group',y='percent',data = ds_2bp, ax=ax, palette=color_palette, order=order,\
		edgecolor='k', size=8, linewidth=1,legend=False)
	ax.set_ylim([-0.1, 1.1])
	plt.savefig('figures/2bp_ds.pdf', dpi=200,  bbox_inches='tight')
	plt.show()

	select1 = ds_2bp['group'] == 'AD'
	select2 = ds_2bp['group'] == 'ctrl'
	res = stats.mannwhitneyu(ds_2bp.loc[select1,'percent'], ds_2bp.loc[select2, 'percent'])
	print(res)

	select1 = ds_2bp['group'] == 'Tau'
	select2 = ds_2bp['group'] == 'ctrl'
	res = stats.mannwhitneyu(ds_2bp.loc[select1,'percent'], ds_2bp.loc[select2, 'percent'])
	print(res)

	select1 = ds_2bp['group'] == 'noTau'
	select2 = ds_2bp['group'] == 'ctrl'
	res = stats.mannwhitneyu(ds_2bp.loc[select1,'percent'], ds_2bp.loc[select2, 'percent'])
	print(res)

	select1 = ds_2bp['group'] == 'Tau'
	select2 = ds_2bp['group'] == 'noTau'
	res = stats.mannwhitneyu(ds_2bp.loc[select1,'percent'], ds_2bp.loc[select2, 'percent'])
	print(res)

	
def main():

	ds_abs = pd.read_csv('results/All_ss_indel_abs.csv', index_col=0, header=0)
	ds_percent = pd.read_csv('results/All_ss_indel_percent.csv', index_col=0, header=0)

	ds_burden = ds_abs.loc[ds_abs['variable']=='Signature_ssIndel_re_burden']
	md = smf.mixedlm(" value ~ Age *  C(group, Treatment('ctrl'))", ds_burden, groups=ds_burden["donor"])
	mdf = md.fit()
	mlr_results = mdf.summary().tables[1].round(8)
	mlr_results['Coef.'] = mlr_results['Coef.'].astype(float)
	print(mlr_results)
	print("\n")
	print(mdf.pvalues)

	plot_df = ds_burden
	color_palette= { 'AD': 'tab:orange', 'ctrl': 'tab:blue', 'Tau':'tab:red', 'noTau':'tab:pink' }	
	label = { 'AD': 'AD Tau agnostic', 'ctrl': 'Control', 'Tau':'P-tau+', 'noTau':'P-tau-' }
	var='value'
	slope_input='relax'

	fig,ax = plt.subplots(constrained_layout=True,figsize=(8, 4), frameon = False)
	for igrp, df_sub in plot_df.groupby('group',sort=False):

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

			if slope_input == 'relax':
				sub_status2 = list(filter(r.match,mlr_results.index.to_list()))[1]
				add_on_slope = mlr_results.loc[sub_status2,'Coef.'] 

			intercept=base_intercept + add_on_intercept
			slope = base_slope + add_on_slope
		ax.plot(xpos, slope * xpos + intercept, c=color_palette[igrp], lw=2, alpha=0.8)
		ax.scatter(x=df_sub['Age'], y=df_sub[var], c=color_palette[igrp], edgecolors='black', \
					label=label[igrp], linewidth=1,s=80)
		ax.set_xlabel('Age')
		ax.set_ylabel('doule-strand sIndel per neuron')


	plt.legend(loc='center left',bbox_to_anchor=(1.0, 0.5), frameon=False)
	plt.savefig('figures/ss_indel_burden.pdf', dpi=200,  bbox_inches='tight')
	plt.show()



if __name__ == "__main__":
#	main()
	percent()






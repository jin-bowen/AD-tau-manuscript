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
import statsmodels.api as sm
import statsmodels.formula.api as smf
import pandas as pd
import sys
import re

def print(*args):
	__builtins__.print(*("%.8f" % a if isinstance(a, float) else a for a in args))

def main():

	burden_in   = '/home/boj924/AD_Tau_PTA/results/tau_AD_ctrl_res.tab'
	sig_in      = '/home/boj924/AD_Tau_PTA/results/tauADCtrl_aging_contribute.tab'
	metafile_in = '/home/boj924/AD_Tau_PTA/metafiles/all_meta.tab'
	clinic_in   = '/home/boj924/AD_Tau_PTA/metafiles/all_clinic'
	outdir      = '/home/boj924/AD_Tau_PTA/figures/'

	muttype     = 'indel'
	var         = 'burden'	
	slope_input = 'relax'

	muttype     = sys.argv[1]
	var         = sys.argv[2]
	slope_input = sys.argv[3]
	compare     = sys.argv[4]
		
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

	if slope_input != 'relax':
		if compare == 'AD-ctrl':	
			out_name = outdir + 'AD_ctrl_%s'%var
			print("results for AD-ctrl")
			burden_meta_sub = burden_meta.loc[~burden_meta['group'].isin(['Tau','noTau'])]
			md = smf.mixedlm(" %s ~ Age +  C(group, Treatment('ctrl')) "%var, burden_meta_sub, groups=burden_meta_sub["donor"])
			mdf = md.fit()
			mlr_results = mdf.summary().tables[1].round(8)
			mlr_results['Coef.'] = mlr_results['Coef.'].astype(float)
			print(mlr_results)
			print("\n")
			print(mdf.pvalues)

		elif compare == 'tau-notau':	
			out_name = outdir + 'Tau-notau_%s'%var
			print("results for Tau-notau")

			burden_meta_sub = burden_meta.loc[burden_meta['group'].isin(['Tau','noTau'])]
			md = smf.mixedlm(" %s ~ Age +  C(group, Treatment('noTau')) "%var, burden_meta_sub, groups=burden_meta_sub["donor"])
			mdf = md.fit()
			mlr_results = mdf.summary().tables[1].round(8)
			mlr_results['Coef.'] = mlr_results['Coef.'].astype(float)
			print(mlr_results)
			print("\n")

			burden_meta_sub = burden_meta.loc[burden_meta['group'].isin(['Tau','noTau','ctrl'])]
			md = smf.mixedlm(" %s ~ Age +  C(group, Treatment('ctrl')) "%var, burden_meta_sub, groups=burden_meta_sub["donor"])
			mdf = md.fit()
			mlr_results = mdf.summary().tables[1].round(8)
			mlr_results['Coef.'] = mlr_results['Coef.'].astype(float)
			print(mlr_results)
			print("\n")
			print(mdf.pvalues)

		elif compare == 'all':
			out_name = outdir + 'all_%s'%var
			print("results for all-ctrl")
			burden_meta_sub = burden_meta
			md = smf.mixedlm(" %s ~ Age +  C(group, Treatment('ctrl')) "%var, burden_meta_sub, groups=burden_meta_sub["donor"])
			mdf = md.fit()
			mlr_results = mdf.summary().tables[1].round(8)
			mlr_results['Coef.'] = mlr_results['Coef.'].astype(float)
			print(mlr_results)
			print("\n")
			print(mdf.pvalues)

	else:	
		if compare == 'AD-ctrl':	
			out_name = outdir + 'AD_ctrl_%s'%var
			print("results for AD-ctrl")
			burden_meta_sub = burden_meta.loc[~burden_meta['group'].isin(['Tau','noTau'])]
			md = smf.mixedlm(" %s ~ Age *  C(group, Treatment('ctrl')) "%var, burden_meta_sub, groups=burden_meta_sub["donor"])
			mdf = md.fit()
			mlr_results = mdf.summary().tables[1].round(8)
			mlr_results['Coef.'] = mlr_results['Coef.'].astype(float)
			print(mlr_results)
			print("\n")
			print(mdf.pvalues)

		elif compare == 'tau-notau':	
			out_name = outdir + 'Tau-notau_%s'%var
			print("results for Tau-notau")
			burden_meta_sub = burden_meta.loc[burden_meta['group'].isin(['Tau','noTau','ctrl'])]
			print(burden_meta_sub)
			md = smf.mixedlm(" %s ~ Age *  C(group, Treatment('ctrl')) "%var, burden_meta_sub, groups=burden_meta_sub["donor"])
			mdf = md.fit()
			mlr_results = mdf.summary().tables[1].round(8)
			mlr_results['Coef.'] = mlr_results['Coef.'].astype(float)
			print(mlr_results)
			print("\n")
			print(mdf.pvalues)

		elif compare == 'all':
			out_name = outdir + 'all_%s'%var
			print("results for all-ctrl")
			burden_meta_sub = burden_meta
			md = smf.mixedlm(" %s ~ Age *  C(group, Treatment('ctrl')) "%var, burden_meta_sub, groups=burden_meta_sub["donor"])
			mdf = md.fit()
			mlr_results = mdf.summary().tables[1].round(8)
			mlr_results['Coef.'] = mlr_results['Coef.'].astype(float)
			print(mlr_results)
			print("\n")
			print(mdf.pvalues)


	plot_df = burden_meta_sub
	color_palette= { 'AD': 'tab:orange', 'ctrl': 'tab:blue', 'Tau':'tab:red', 'noTau':'tab:pink' }	
	label = { 'AD': 'AD Tau agnostic', 'ctrl': 'Control', 'Tau':'P-tau+', 'noTau':'P-tau-' }

	fig,ax = plt.subplots(constrained_layout=True,figsize=(3,4), frameon = False)
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
		ax.set_ylabel('sSNV per neuron')

	plt.legend(loc='center left',bbox_to_anchor=(1.0, 0.5), frameon=False)
	plt.savefig(out_name+'_'+muttype+'-'+var+'_%s.pdf'%slope_input, dpi=200,  bbox_inches='tight')
	plt.show()

if __name__ == "__main__":
	main()




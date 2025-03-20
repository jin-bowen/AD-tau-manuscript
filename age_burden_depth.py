import matplotlib.pyplot as plt
import matplotlib.font_manager
plt.rcParams['font.size'] = 14
plt.rcParams['axes.linewidth'] = 2
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
	metafile_in = '/home/boj924/AD_Tau_PTA/metafiles/all_meta.tab'
	clinic_in   = '/home/boj924/AD_Tau_PTA/metafiles/all_clinic'
	out_name    = 'AD_ctrl_burden'
	muttype     = 'indel'
	depth_profile_in = '/home/boj924/AD_Tau_PTA/AD_depth_qc/%s_depth_quantile.txt'%muttype
	
	color_palette = 'Pastel1'
	burden = pd.read_csv(burden_in, sep=',', header=0,index_col=0)
	burden = burden.loc[burden['muttype']==muttype]
	burden.dropna(inplace=True)
	burden['sample'] = burden.index

	metafile = pd.read_csv(metafile_in, header=0,dtype=str)
	clinic = pd.read_csv(clinic_in, header=0, sep='\t',dtype=str)
	metafile = pd.merge(metafile,clinic, left_on='donor',right_on='Case_ID')
	metafile['Age'] = metafile['Age_(yrs)'].astype(float)
	burden_meta = pd.merge(burden,metafile, on='sample', suffixes=('','_meta'))

	var='burden'
	print("results for AD-ctrl")
	burden_meta_AD_ctrl = burden_meta.loc[~burden_meta['group'].isin(['Tau','noTau'])]
	md = smf.mixedlm(" %s ~ Age +  C(group, Treatment('ctrl')) "%var, burden_meta_AD_ctrl, groups=burden_meta_AD_ctrl["donor"])
	mdf = md.fit()
	mlr_results = mdf.summary().tables[1].round(8)
	mlr_results['Coef.'] = mlr_results['Coef.'].astype(float)
	print(mlr_results)
	print("\n")
	print(mdf.pvalues)

#	print("results for Tau-notau")
#	burden_meta_Tau = burden_meta.loc[burden_meta['group'].isin(['Tau','noTau'])]
#	md = smf.mixedlm(" %s ~ Age +  C(group, Treatment('noTau')) "%var, burden_meta_Tau, groups=burden_meta_Tau["donor"])
#	mdf = md.fit()
#	mlr_results = mdf.summary().tables[1].round(8)
#	mlr_results['Coef.'] = mlr_results['Coef.'].astype(float)
#	print(mlr_results)
#	print("\n")
#	print(mdf.pvalues)
#
#	print("results for Tau-ctrl")
#	burden_meta_Tau_ctrl = burden_meta.loc[~burden_meta['group'].isin(['AD'])]
#	md = smf.mixedlm(" %s ~ Age +  C(group, Treatment('ctrl')) "%var, burden_meta_Tau_ctrl, groups=burden_meta_Tau_ctrl["donor"])
#	mdf = md.fit()
#	mlr_results = mdf.summary().tables[1].round(8)
#	mlr_results['Coef.'] = mlr_results['Coef.'].astype(float)
#	print(mlr_results)
#	print("\n")
#	print(mdf.pvalues)

#	var='burden'
#	print("results for AD-ctrl")
#	md = smf.mixedlm(" %s ~ Age +  C(group, Treatment('ctrl')) "%var, burden_meta, groups=burden_meta["donor"])
#	mdf = md.fit()
#	mlr_results = mdf.summary().tables[1].round(8)
#	mlr_results['Coef.'] = mlr_results['Coef.'].astype(float)
#	print(mlr_results)
#	print("\n")
#	print(mdf.pvalues)

	depth_profile = pd.read_csv(depth_profile_in, header=0, index_col=0)
	depth_profile['class'] = depth_profile.apply(lambda row: 0 if row['25%'] <10 else 1, axis=1)

	plot_df = burden_meta_AD_ctrl
	color_palette= { 'AD': 'tab:orange', 'ctrl': 'tab:blue', 'Tau':'tab:red', 'noTau':'tab:pink' }	
	label = { 'AD': 'AD', 'ctrl': 'Control', 'Tau':'P-tau+', 'noTau':'P-tau-' }

	color_palette2 = { 0:'black', 1: 'white'}

	fig,ax = plt.subplots(constrained_layout=True,figsize=(6,4), frameon = False)
	for igrp, df_sub in plot_df.groupby('group',sort=False):
		xpos=np.arange(57, 91)
		base_intercept = mlr_results.loc['Intercept','Coef.']
		slope = mlr_results.loc['Age','Coef.']
		if igrp == 'ctrl':
			intercept = base_intercept
			xpos=np.arange(0, 110)
		else:
			r = re.compile(".*%s.*"%igrp)
			sub_status = list(filter(r.match,mlr_results.index.to_list()))[0]
			add_on_intercept = mlr_results.loc[sub_status,'Coef.'] 
			intercept=base_intercept + add_on_intercept

		ax.plot(xpos, slope * xpos + intercept, c=color_palette[igrp], lw=2, alpha=0.8)

		if igrp == 'ctrl':
			ax.scatter(x=df_sub['Age'], y=df_sub[var], c=color_palette[igrp], edgecolors='black', \
					linewidth=2, label=label[igrp])
		else:
			df_sub = pd.merge(df_sub, depth_profile, left_on='sample', right_index=True)
			clist = [ color_palette2[x]  for x  in df_sub['class'] ]
			ax.scatter(x=df_sub['Age'], y=df_sub[var], c=clist, edgecolors='black', \
					linewidth=2)		

		ax.set_xlabel('Age')
		ax.set_ylabel('Burden')

	plt.legend(loc='center left',bbox_to_anchor=(1.01, 0.5), frameon=False)
	plt.savefig(out_name+'_'+muttype+'-'+var+'.png', dpi=200,  bbox_inches='tight')
	plt.show()	
	plt.close()


if __name__ == "__main__":
	main()



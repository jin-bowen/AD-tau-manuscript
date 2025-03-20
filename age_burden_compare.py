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

	lira_in   = '/home/boj924/AD_Tau_PTA/metafiles/control_metafile_2022nature' 

	burden_scan2_in = '/home/boj924/AD_Tau_PTA/results/tau_AD_ctrl_res.tab'
	metafile_in = '/home/boj924/AD_Tau_PTA/metafiles/all_meta.tab'
	clinic_in = '/home/boj924/AD_Tau_PTA/metafiles/all_clinic'
	color_palette = 'Pastel1'
	out_name='caller'

	burden_lira = pd.read_csv(lira_in, sep=',', header=0)
	burden_lira['burden'] = burden_lira['Somatic_rate'] * 5.845001134

	burden_scan2 = pd.read_csv(burden_scan2_in, sep=',', header=0,index_col=0)
	burden_scan2 = burden_scan2.loc[burden_scan2['mutation_type']=='snv']
	burden_scan2.dropna(inplace=True)
	burden_scan2['sample'] = burden_scan2.index
	burden_scan2 = burden_scan2.loc[burden_scan2['sample']!='1995P_201001E3']

	metafile = pd.read_csv(metafile_in, header=0,dtype=str)
	clinic = pd.read_csv(clinic_in, header=0, sep='\t',dtype=str)
	metafile = pd.merge(metafile, clinic, left_on='donor', right_on='Case_ID')
	metafile['Age'] = metafile['Age_(yrs)'].astype(float)
	burden_scan2_meta = pd.merge(burden_scan2,metafile, on='sample', suffixes=('','_meta'))
	burden_scan2_meta['caller'] = 'SCAN2'

	burden_lira['Case_ID'] =  burden_lira['Case_ID'].astype(str)
	burden_lira_meta = pd.merge(burden_lira,metafile[['donor','group']].drop_duplicates(), 
				left_on='Case_ID', right_on='donor', suffixes=('','_meta'))	
	burden_lira_meta['caller'] = 'LiRA'

	print(burden_lira_meta.head())
	burden_lira_meta['Cell_ID'].to_csv('lira.sample')
	burden_scan2_meta['sample'].to_csv('scan2.sample')

	df_scan2 = burden_scan2_meta[['Age','Case_ID','burden','group','caller']]
	df_lira = burden_lira_meta[['Age','Case_ID','burden','group', 'caller']]

	df = pd.concat([df_scan2, df_lira])

	df['Age'] = df['Age'].astype(float)
	df['Case_ID'] = df['Case_ID'].astype(str)
	df['burden'] = df['burden'].astype(float)
	df = df.drop_duplicates()

	group='ctrl'
	df_sub = df.loc[df['group']==group]

	var='burden'
	md = smf.mixedlm(" %s ~ Age * C(caller, Treatment('LiRA'))"%var, df_sub, groups=df_sub["Case_ID"])
	mdf = md.fit()
	mlr_results = mdf.summary().tables[1].round(8)
	mlr_results['Coef.'] = mlr_results['Coef.'].astype(float)
	print(mlr_results)
	print(mdf.pvalues)

	plot_df = df_sub
	color_palette= {  'SCAN2':'tab:green',  'LiRA': 'gray'}	

	fig,ax = plt.subplots(constrained_layout=True,figsize=(7,4), frameon = False)
	for igrp, df_sub in plot_df.groupby('caller',sort=False):

		xpos=np.arange(57,91)
		base_intercept = mlr_results.loc['Intercept','Coef.']
		base_slope = mlr_results.loc['Age','Coef.']
		if igrp == 'LiRA':
			intercept = base_intercept
			slope = base_slope
		else:
			r = re.compile(".*%s.*"%igrp)
			sub_status = list(filter(r.match,mlr_results.index.to_list()))[0]
			add_on_intercept = mlr_results.loc[sub_status,'Coef.'] 

			sub_status2 = list(filter(r.match,mlr_results.index.to_list()))[1]
			add_on_slope = mlr_results.loc[sub_status2,'Coef.']

			intercept=base_intercept + add_on_intercept
			slope = base_slope + add_on_slope 

		ax.plot(xpos, slope * xpos + intercept, c=color_palette[igrp], lw=2, alpha=0.8)
		ax.scatter(x=df_sub['Age'], y=df_sub[var], c=color_palette[igrp], edgecolors='black', \
					label=igrp, linewidth=2)

		ax.set_xlabel('Age')
		ax.set_ylabel('sSNVs per neuron')

	ax.set_xlim([0,100])
	ax.set_ylim([0,4000])

	plt.legend(loc='center left',bbox_to_anchor=(1.01, 0.5), frameon=False)
	plt.savefig(out_name+'_'+group+'.png', dpi=200,  bbox_inches='tight')
	plt.show()	
	plt.close()

if __name__ == "__main__":
	main()



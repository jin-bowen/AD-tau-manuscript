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

def indel_prop():
	
	tbp_count = '/home/boj924/AD_Tau_PTA/results/indel_spectrum_count.tab'
	burden_in = '/home/boj924/AD_Tau_PTA/results/tau_AD_ctrl_res.tab'

	metafile_in = '/home/boj924/AD_Tau_PTA/metafiles/all_meta.tab'
	clinic_in   = '/home/boj924/AD_Tau_PTA/metafiles/all_clinic'
	out_name    = '/home/boj924/AD_Tau_PTA/results/indel'

	burden = pd.read_csv(burden_in, sep=',', header=0,index_col=0)
	burden = burden.loc[burden['mutation_type']=='indel',]
	
	tbp_count_df = pd.read_csv(tbp_count, index_col=0, header=0)

	tbp_col = ['2bp_deletion_1','2bp_deletion_2', \
			'2bp_deletion_3','2bp_deletion_4','2bp_deletion_5','2bp_deletion_6+', \
			'2bp_deletion_with_microhomology_1']

	tbp_pct = tbp_count_df/tbp_count_df.sum(axis=1).values.reshape((-1,1))
	tbp_pct['twobp_burden'] = tbp_pct[tbp_col].sum(axis=1)

	color_palette = 'Pastel1'

	metafile = pd.read_csv(metafile_in, header=0,dtype=str)
	clinic = pd.read_csv(clinic_in, header=0, sep='\t',dtype=str)
	metafile = pd.merge(metafile,clinic, left_on='donor',right_on='Case_ID')
	metafile['Age'] = metafile['Age_(yrs)'].astype(float)

	df = pd.merge(tbp_pct, metafile, left_index=True, right_on='sample')
	var = 'twobp_burden'
	
	color_palette= { 'AD': 'tab:orange', 'ctrl': 'tab:blue', 'Tau':'tab:red', 'noTau':'tab:pink' }	
	label = { 'AD': 'AD', 'ctrl': 'Control', 'Tau':'P-tau+', 'noTau':'P-tau-' }

	df_sub = df.loc[df['group'].isin(['ctrl','AD'])]
	mean_df = df_sub.groupby('group')['twobp_burden'].mean().to_frame().reset_index()

	fig,ax = plt.subplots(constrained_layout=True,figsize=(4,5), frameon = False)
	ax.bar(x=mean_df['group'], height=mean_df['twobp_burden'], color='tab:green')
	ax.set_ylabel('% of 2bp deletion to total indel')

	plt.legend(loc='center left',bbox_to_anchor=(1.01, 0.5), frameon=False)
	plt.savefig(out_name+'_'+var+'.png', dpi=200,  bbox_inches='tight')
	plt.show()	
	plt.close()


def scale_indel():
	
	tbp_count = '/home/boj924/AD_Tau_PTA/results/indel_spectrum_count.tab'
	burden_in = '/home/boj924/AD_Tau_PTA/results/tau_AD_ctrl_res.tab'

	metafile_in = '/home/boj924/AD_Tau_PTA/metafiles/all_meta.tab'
	clinic_in   = '/home/boj924/AD_Tau_PTA/metafiles/all_clinic'
	AT8pct_in   = '/home/boj924/AD_Tau_PTA/metafiles/AT8.txt'
	outdir      = '/home/boj924/AD_Tau_PTA/figures/'
	compare     = 'tau-notau'

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
	var = 'twobp_burden'

	burden_meta_sub = burden_meta.loc[burden_meta['group'].isin(['Tau','noTau','ctrl'])]
	md = smf.mixedlm(" %s ~ Age + C(group, Treatment('ctrl')) "%var, burden_meta_sub, groups=burden_meta_sub["donor"])
	mdf = md.fit()
	mlr_results = mdf.summary().tables[1].round(8)
	mlr_results['Coef.'] = mlr_results['Coef.'].astype(float)
	print(mlr_results)
	print("\n")
	print(mdf.pvalues)

	burden_meta_sub = burden_meta.loc[burden_meta['group'].isin(['Tau','noTau'])]
	md = smf.mixedlm(" %s ~ Age + C(group, Treatment('noTau')) "%var, burden_meta_sub, groups=burden_meta_sub["donor"])
	mdf = md.fit()
	mlr_results = mdf.summary().tables[1].round(8)
	mlr_results['Coef.'] = mlr_results['Coef.'].astype(float)
	print(mlr_results)
	print("\n")
	print(mdf.pvalues)

	if compare == 'AD-ctrl':
		out_name = outdir + 'AD_ctrl_%s'%var
		print("results for AD-ctrl")
		burden_meta_sub = burden_meta.loc[~burden_meta['group'].isin(['Tau','noTau'])]
		md = smf.mixedlm(" %s ~ Age *  C(group, Treatment('ctrl')) "%var, burden_meta_sub, groups=burden_meta_sub["donor"])
		mdf = md.fit()
		mlr_results = mdf.summary().tables[1].round(8)
		mlr_results['Coef.'] = mlr_results['Coef.'].astype(float)

	elif compare == 'tau-notau':
		out_name = outdir + 'Tau-notau_%s'%var
		print("results for Tau-notau")
		burden_meta_sub = burden_meta.loc[burden_meta['group'].isin(['Tau','noTau','ctrl'])]
		md = smf.mixedlm(" %s ~ Age *  C(group, Treatment('ctrl')) "%var, burden_meta_sub, groups=burden_meta_sub["donor"])
		mdf = md.fit()
		mlr_results = mdf.summary().tables[1].round(8)
		mlr_results['Coef.'] = mlr_results['Coef.'].astype(float)

	elif compare == 'all':
		out_name = outdir + 'all_%s'%var
		print("results for all-ctrl")
		burden_meta_sub = burden_meta
		md = smf.mixedlm(" %s ~ Age *  C(group, Treatment('ctrl')) "%var, burden_meta_sub, groups=burden_meta_sub["donor"])
		mdf = md.fit()
		mlr_results = mdf.summary().tables[1].round(8)
		mlr_results['Coef.'] = mlr_results['Coef.'].astype(float)

#	print(mlr_results)
#	print("\n")
#	print(mdf.pvalues)

	plot_df = burden_meta_sub
	color_palette= { 'AD': 'tab:orange', 'ctrl': 'tab:blue', 'Tau':'tab:red', 'noTau':'tab:pink' }
	label = { 'AD': 'AD', 'ctrl': 'Control', 'Tau':'P-tau+', 'noTau':'P-tau-' }
	fig,ax = plt.subplots(constrained_layout=True,figsize=(8,5), frameon = False)
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
			sub_status = list(filter(r.match,mlr_results.index.to_list()))[1]
			add_on_slope = mlr_results.loc[sub_status,'Coef.']
			add_on_intercept = mlr_results.loc[sub_status,'Coef.'] 

			slope = base_slope + add_on_slope
			intercept=base_intercept + add_on_intercept

		ax.plot(xpos, slope * xpos + intercept, c=color_palette[igrp], lw=2, alpha=0.8)
		ax.scatter(x=df_sub['Age'], y=df_sub['twobp_burden'], c=color_palette[igrp], edgecolors='black', \
					label=label[igrp], linewidth=2)

		ax.set_xlabel('Age')
		ax.set_ylabel('2bp deletion burden')

	plt.legend(loc='center left',bbox_to_anchor=(1.01, 0.5), frameon=False)
	plt.savefig(out_name+'_'+var+'.pdf', dpi=300,  bbox_inches='tight')
	plt.show()	
	plt.close()


def scale_ID4():
	
	tbp_count = '/home/boj924/AD_Tau_PTA/results/indel_spectrum_count.tab'
	ID4_count = '/home/boj924/AD_Tau_PTA/results/ID83/Suggested_Solution/Assignment_Solution/Activities/Assignment_Solution_Activities.txt'
	burden_in = '/home/boj924/AD_Tau_PTA/results/tau_AD_ctrl_res.tab'

	metafile_in = '/home/boj924/AD_Tau_PTA/metafiles/all_meta.tab'
	clinic_in   = '/home/boj924/AD_Tau_PTA/metafiles/all_clinic'
	AT8pct_in   = '/home/boj924/AD_Tau_PTA/metafiles/AT8.txt'
	out_name    = '/home/boj924/AD_Tau_PTA/results/indel'

	burden = pd.read_csv(burden_in, sep=',', header=0,index_col=0)
	burden = burden.loc[burden['muttype']=='indel',]
	
	ID4_count_df = pd.read_csv(ID4_count, sep='\t', header=0, index_col=0)
	ID4_count_df.index = ID4_count_df.index.to_series().apply(lambda x: x.strip('_indel'))
	ID4_pct = ID4_count_df/ID4_count_df.sum(axis=1).values.reshape((-1,1))
	ID4_burden = pd.merge(ID4_pct, burden[['burden']], left_index=True, right_index=True)
	ID4_burden[ID4_pct.columns] = ID4_burden[ID4_pct.columns] * ID4_burden['burden'].values.reshape((-1,1))

	tbp_count_df = pd.read_csv(tbp_count, index_col=0, header=0)

	tbp_col = ['2bp_deletion_1','2bp_deletion_2', \
			'2bp_deletion_3','2bp_deletion_4','2bp_deletion_5','2bp_deletion_6+', \
			'2bp_deletion_with_microhomology_1']

	tbp_pct = tbp_count_df/tbp_count_df.sum(axis=1).values.reshape((-1,1))
	tbp_burden = pd.merge(tbp_pct, burden[['burden']], left_index=True, right_index=True)
	tbp_burden[tbp_pct.columns] = tbp_burden[tbp_pct.columns] * tbp_burden['burden'].values.reshape((-1,1))
	tbp_burden['twobp_burden'] = tbp_burden[tbp_col].sum(axis=1)

	burden_df = pd.merge(ID4_burden, tbp_burden, left_index=True, right_index=True,suffixes=('','_y'))
	color_palette = 'Pastel1'

	metafile = pd.read_csv(metafile_in, header=0,dtype=str)
	clinic = pd.read_csv(clinic_in, header=0, sep='\t',dtype=str)
	metafile = pd.merge(metafile,clinic, left_on='donor',right_on='Case_ID')
	metafile['Age'] = metafile['Age_(yrs)'].astype(float)

	df = pd.merge(burden_df, metafile, left_index=True, right_on='sample')
	df = df.loc[df['Case_ID']!='1995']

	var = 'twobp_burden'
#	var = 'ID4'
	print(df)

	plot_df = df
	color_palette= { 'AD': 'tab:orange', 'ctrl': 'tab:blue', 'Tau':'tab:red', 'noTau':'tab:pink' }
	label = { 'AD': 'AD', 'ctrl': 'Control', 'Tau':'P-tau+', 'noTau':'P-tau-' }

	fig,ax = plt.subplots(constrained_layout=True,figsize=(8,5), frameon = False)
	for igrp, df_sub in plot_df.groupby('group',sort=False):

#		xpos=np.arange(57, 91)	
#		base_intercept = results.loc['Intercept','Coef.']
#		slope = results.loc['Age','Coef.']  
#		if igrp == 'ctrl':
#			intercept = base_intercept
#			xpos=np.arange(0, 110)	
#		else:
#			r = re.compile(".*%s.*"%igrp)
#			sub_status = list(filter(r.match,results.index.to_list()))[0]
#			add_on_intercept = results.loc[sub_status,'Coef.'] 
#			intercept=base_intercept + add_on_intercept
#
#		print(igrp, intercept, slope)
#		ax.plot(xpos, slope * xpos + intercept, c=color_palette[igrp], lw=2, alpha=0.8)
		print(df_sub['twobp_burden'], df_sub['Age'])
		print(df_sub.columns)
		ax.scatter(x=df_sub['Age'], y=df_sub['twobp_burden'], c=color_palette[igrp], edgecolors='black', \
					label=label[igrp], linewidth=2)

		ax.set_xlabel('Age')
		ax.set_ylabel('2bp deletion burden')

	plt.legend(loc='center left',bbox_to_anchor=(1.01, 0.5), frameon=False)
	plt.savefig(out_name+'_'+var+'.png', dpi=200,  bbox_inches='tight')
	plt.show()	
	plt.close()


if __name__ == "__main__":
	scale_indel()
#	indel_prop()



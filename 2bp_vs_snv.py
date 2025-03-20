import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib as mpl
import matplotlib.font_manager
plt.rcParams['font.size'] = 10
plt.rcParams['axes.linewidth'] = 1
plt.style.use('fivethirtyeight')
import numpy as np
from scipy import stats
from scipy.interpolate import UnivariateSpline,CubicSpline
from matplotlib.lines import Line2D
import pandas as pd
import sys
import re


def main2():

	deletion_in = '/home/boj924/AD_Tau_PTA/results/All_indel_spectrum.tab'
	burden_in   = '/home/boj924/AD_Tau_PTA/results/burden_meta.res'
	sig_in      = '/home/boj924/AD_Tau_PTA/results/contribute_mat_burden_meta.res'
	metafile_in = '/home/boj924/AD_Tau_PTA/metafiles/all_meta.tab'
	out_name = '2bp_snv'

	deletion = pd.read_csv(deletion_in, sep=',', header=0, index_col=0)
	metafile = pd.read_csv(metafile_in, sep=',', header=0)
	burden = pd.read_csv(burden_in, sep=',', header=0, index_col=0)
	sig_burden = pd.read_csv(sig_in, sep=',', header=0, index_col=0)

	indel_burden = burden.loc[burden['muttype']=='indel']

	deletion_sum = deletion.sum(axis=0)
	deletion = deletion/deletion_sum

	deletion.reset_index(inplace=True)
	deletion_reform = pd.melt(deletion, id_vars='index',value_vars=set(deletion.columns)-set(['index']))

	deletion_burden = pd.merge(deletion_reform, indel_burden, left_on='variable', right_on='sample')
	deletion_burden['burden'] = deletion_burden.apply(lambda x: x['value'] * x['burden'], axis=1)

	select_index = ['2bp_deletion_1','2bp_deletion_2','2bp_deletion_with_microhomology_1']
	deletion_sub = deletion_burden.loc[deletion_burden['index'].isin(select_index)]
	deletion_2bp_sum = deletion_sub.groupby(['donor','sample','group'])[['burden']].sum().reset_index()


	for var in ['scale_sigA_burden', 'scale_sigC_burden']:
		snv_burden = sig_burden[['donor','sample','group',var]]
	
		all_df = pd.merge(deletion_2bp_sum,snv_burden, on='sample', suffixes=('_2bp','_snv'))
		all_df_sub = all_df.loc[all_df['group_2bp']!='ctrl']
		print(all_df_sub)
		mean_df_sub = all_df_sub.groupby(['donor_2bp','group_2bp'])[var,'burden'].mean().reset_index()
		var_df_sub = all_df_sub.groupby(['donor_2bp','group_2bp'])[var,'burden'].std().reset_index()
	
		color_dict= { 'AD-2bp': 'tab:orange', 'Tau-2bp':'tab:red', 'noTau-2bp':'tab:pink',\
				'AD-snv': 'tab:orange', 'Tau-snv':'tab:red', 'noTau-snv':'tab:pink' }		
	
		order=['1353','1647','2208','1456','2207','1995','1828']
	
		fig,ax = plt.subplots(constrained_layout=True,figsize=(10,5),frameon = False)	
		customPalette = sns.set_palette(sns.color_palette(['tab:orange','tab:Red','tab:pink']))
	
		corr_res = pd.DataFrame(columns=['donor','corr','p-val'])
		for igrp, grp in all_df_sub.groupby('donor_2bp'):
			corr = stats.pearsonr(grp['burden'],grp[var])
			corr_res.loc[len(corr_res)]=[igrp,corr[0], corr[1]]
	
		print(corr_res)
		p2 = sns.relplot(data=all_df_sub, x=var, y="burden", col_wrap=4, \
			col="donor_2bp", hue="group_2bp", kind="scatter", s=30, ax=ax)
	
		#plt.ylim([-10,5000])
		#plt.xlim([500,2500])
	
		plt.xlabel(var)
		plt.ylabel('2bp burden')
		ax.legend(loc='center left',bbox_to_anchor=(1.01, 0.5), frameon=False)	
	#	ax.legend(*p1.legend_elements(),loc='center left',bbox_to_anchor=(1.01, 0.5), frameon=False)
		plt.savefig(out_name+var+'.png')	
		plt.show()

def main():

	deletion_in = '/home/boj924/AD_Tau_PTA/results/All_indel_spectrum.tab'
	burden_in   = '/home/boj924/AD_Tau_PTA/results/burden_meta.res'
	sig_in      = '/home/boj924/AD_Tau_PTA/results/contribute_mat_burden_meta.res'
	metafile_in = '/home/boj924/AD_Tau_PTA/metafiles/all_meta.tab'
	out_name = '2bp_snv'

	deletion = pd.read_csv(deletion_in, sep=',', header=0, index_col=0)
	burden = pd.read_csv(burden_in, sep=',', header=0, index_col=0)
	metafile = pd.read_csv(metafile_in, sep=',', header=0)
	sig_burden = pd.read_csv(sig_in, sep=',', header=0, index_col=0)

	indel_burden = burden.loc[burden['muttype']=='indel']

	deletion_sum = deletion.sum(axis=0)
	deletion = deletion/deletion_sum

	deletion.reset_index(inplace=True)
	deletion_reform = pd.melt(deletion, id_vars='index',value_vars=set(deletion.columns)-set(['index']))

	deletion_burden = pd.merge(deletion_reform, indel_burden, left_on='variable', right_on='sample')
	deletion_burden['burden'] = deletion_burden.apply(lambda x: x['value'] * x['burden'], axis=1)

	select_index = ['2bp_deletion_1','2bp_deletion_2','2bp_deletion_with_microhomology_1']
	deletion_sub = deletion_burden.loc[deletion_burden['index'].isin(select_index)]
	deletion_2bp_sum = deletion_sub.groupby(['donor','sample','group'])[['burden']].sum().reset_index()

	snv_burden = burden.loc[burden['muttype']=='snv']
	snv_burden_sub = snv_burden[['donor','sample','group','burden']]

	all_df = pd.merge(deletion_2bp_sum,snv_burden_sub, on='sample', suffixes=('_2bp','_snv'))
	all_df_sub = all_df.loc[all_df['group_2bp']!='ctrl']
	mean_df_sub = all_df_sub.groupby(['donor_2bp','group_2bp'])['burden_snv','burden_2bp'].mean().reset_index()
	var_df_sub = all_df_sub.groupby(['donor_2bp','group_2bp'])['burden_snv','burden_2bp'].std().reset_index()
#	all_df['combine_feature'] = all_df[['group','type']].agg(lambda x: '-'.join(x), axis=1) 

	color_dict= { 'AD-2bp': 'tab:orange', 'Tau-2bp':'tab:red', 'noTau-2bp':'tab:pink',\
			'AD-snv': 'tab:orange', 'Tau-snv':'tab:red', 'noTau-snv':'tab:pink' }		

	order=['1353','1647','2208','1456','2207','1995','1828']

	fig,ax = plt.subplots(constrained_layout=True,figsize=(10,5),frameon = False)	
	customPalette = sns.set_palette(sns.color_palette(['tab:orange','tab:Red','tab:pink']))
#	p1 = sns.scatterplot(mean_df_sub, x='burden_snv',y='burden_2bp', \
#		hue='group_2bp', s=40, palette=customPalette)

	corr_res = pd.DataFrame(columns=['donor','corr','p-val'])
	for igrp, grp in all_df_sub.groupby('donor_2bp'):
		corr = stats.pearsonr(grp['burden_2bp'],grp['burden_snv'])
		corr_res.loc[len(corr_res)]=[igrp,corr[0], corr[1]]

	print(corr_res)
	p2 = sns.relplot(data=all_df_sub, x="burden_snv", y="burden_2bp", col_wrap=4, \
		col="donor_2bp", hue="group_2bp", kind="scatter", s=30, ax=ax)

	#plt.ylim([-10,5000])
	#plt.xlim([500,2500])

	plt.xlabel('sSNV burden')
	plt.ylabel('2bp burden')
	ax.legend(loc='center left',bbox_to_anchor=(1.01, 0.5), frameon=False)	
#	ax.legend(*p1.legend_elements(),loc='center left',bbox_to_anchor=(1.01, 0.5), frameon=False)
	plt.savefig(out_name+'.png')	
	plt.show()

if __name__ == "__main__":
	main2()



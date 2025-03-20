import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib as mpl
import matplotlib.font_manager
plt.rcParams['font.size'] = 10
plt.rcParams['axes.linewidth'] = 1
import numpy as np
from scipy import stats
from scipy.interpolate import UnivariateSpline,CubicSpline
from matplotlib.lines import Line2D
import pandas as pd
import sys
import re

def pct():
	infile='/home/boj924/AD_Tau_PTA/results/ID83/Assignment_Solution/Activities/Assignment_Solution_Activities.txt'
	metafile = '/home/boj924/AD_Tau_PTA/metafiles/all_meta.tab'
	meta_df = pd.read_csv(metafile, header=0)

	df = pd.read_csv(infile, sep='\t', header=0)
	df['Samples'] = df['Samples'].str.strip('_indel')
	df.set_index('Samples', inplace=True)
	df = df/ df.sum(axis=1).values.reshape((-1, 1))
	df['Samples'] = df.index
	df_reform = df.melt(id_vars='Samples', var_name='IDsig', value_name='pct')
	df_reform = pd.merge(df_reform, meta_df, left_on='Samples', right_on='sample')

	for igrp, grp in df_reform.groupby('group'):

		print(len(grp))	
		fig, ax=plt.subplots(constrained_layout=True, figsize=(8, 4) )
		sns.histplot(data=grp, x="Samples", hue="IDsig", weights='pct', \
				palette=colors, multiple="stack", ax=ax)
		ax.set_xlabel('')
		ax.set_ylabel('percentage')
		plt.xticks(rotation=90)
		legend = ax.get_legend()
		legend.set_bbox_to_anchor((1, 1))
		plt.savefig('/home/boj924/AD_Tau_PTA/results/AD_%s_Decompose_pct.pdf'%igrp )
		plt.show()

def snv():

	infile='/home/boj924/AD_Tau_PTA/results/SBS96/Assignment_Solution/Activities/Assignment_Solution_Activities.txt'
#	infile='/home/boj924/AD_Tau_PTA/results/SBS96/Decompose_Solution/Activities/Decompose_Solution_Activities.txt'
	burden='/home/boj924/AD_Tau_PTA/results/tau_AD_ctrl_res.tab'
	burden_df = pd.read_csv(burden, header=0, sep=',')
	burden_df = burden_df.loc[burden_df['mutation_type']=='snv']

	metafile = '/home/boj924/AD_Tau_PTA/metafiles/all_meta.tab'
	meta_df = pd.read_csv(metafile, header=0)

	df = pd.read_csv(infile, sep='\t', header=0)
	df['Samples'] = df['Samples'].str.strip('_snv')
	df.set_index('Samples', inplace=True)
	df = df/ df.sum(axis=1).values.reshape((-1, 1))

	select_sbs = df.sum(axis=0)>0
	df = df.loc[:,select_sbs.values]
	df['Samples'] = df.index
	df_reform = df.melt(id_vars='Samples', var_name='SBSsig', value_name='pct')
	df_reform = pd.merge(df_reform, meta_df, left_on='Samples', right_on='sample')
	df_reform = pd.merge(df_reform, burden_df[['Cell_ID','burden']], left_on='sample', 
				right_on='Cell_ID')	
	df_reform['SBSsig_burden'] = df_reform.apply(lambda row:  row['pct'] * row['burden'], axis=1)

	for igrp, grp in df_reform.groupby('group'):

		fig, ax=plt.subplots(constrained_layout=True, figsize=(8, 4) )
		sns.histplot(data=grp, x="Samples", hue="SBSsig", weights='SBSsig_burden', \
				palette=sbs_color_code, multiple="stack", ax=ax)
		ax.set_xlabel('')
		ax.set_ylabel('Genome-wide burden')
		plt.xticks(rotation=90)
		legend = ax.get_legend()
		legend.set_bbox_to_anchor((1, 1))
		legend.set_frame_on(False)
		ax.set_ylim([0,3000])
		plt.savefig('/home/boj924/AD_Tau_PTA/results/AD_snv_%s_Decompose_count_rescale.pdf'%igrp, bbox_inches='tight' )


def main():

	infile='/n/data1/bwh/pathology/miller/lab/AD_Tau/results_hg19/annovar/Tau_vs_noTau.exonic_intronic.GO.tsv'
	df = pd.read_csv(infile, sep='\t', header=0)
	df['log_corrected_p_Tau'] = df['corrected_p_Tau'].transform(np.log10)
	df['log_corrected_p_Tau'] = df['log_corrected_p_Tau'] * -1
	df['log_corrected_p_noTau'] = df['corrected_p_noTau'].transform(np.log10)
	df['log_corrected_p_noTau'] = df['log_corrected_p_noTau'] * -1
	df['sub_type'] = 'None'

	df['sub_type'] = df.apply(lambda row: 'synapse' if re.search('synap',row['term']) else row['sub_type'], axis=1)
	df['sub_type'] = df.apply(lambda row: 'neuron' if re.search('neuron',row['term']) else row['sub_type'], axis=1)
	df['sub_type'] = df.apply(lambda row: 'kinase' if re.search('kinase',row['term']) else row['sub_type'], axis=1)
	df['sub_type'] = df.apply(lambda row: 'phospho' if re.search('phospho',row['term']) else row['sub_type'], axis=1)
#	df['sub_type'] = df.apply(lambda row: 'None' if row['log_corrected_p_Tau'] <= 2 or row['log_corrected_p_noTau']  <= 2  else row['sub_type'], axis=1)
	df = df.loc[(df['log_corrected_p_Tau'] >= 2) & (df['log_corrected_p_noTau']  >= 2)]

	tau_gene_list = df.sort_values(by='log_corrected_p_Tau', ascending=False)
	notau_gene_list = df.sort_values(by='log_corrected_p_noTau', ascending=False)
	top_gene = tau_gene_list['category'].tolist() +  notau_gene_list['category'].tolist()
	share_gene=set(tau_gene_list['category'].tolist()).intersection(set(notau_gene_list['category'].tolist()))
	
	df_sub = df.loc[df['category'].isin(top_gene),]

	#print(df.loc[df['category'].isin(share_gene),'term'])

	metafile = '/home/boj924/AD_Tau_PTA/metafiles/all_meta.tab'
	meta_df = pd.read_csv(metafile, header=0)

	fig, ax=plt.subplots(constrained_layout=True, figsize=(5,3) )

	color = {'synapse':'gold','neuron':'violet','kinase':'gray','phospho':'gray','None':'gray'}
	sns.scatterplot(data=df_sub, x="log_corrected_p_Tau", y="log_corrected_p_noTau", hue="sub_type", 
			palette=color, edgecolors='k', linewidths=1)
	ax.set_xlabel('P-tau+')
	ax.set_ylabel('P-tau-')
	ax.set_xlim([1.8,8])
	ax.set_ylim([1.8,8])

	xpos = np.arange(1.8, 8)
	ypos = np.arange(1.8, 8)
	ax.plot(xpos,ypos,lw=1, c='k')
	ax.axhline(y=2, lw=1, ls='--', c='k')
	ax.axvline(x=2, lw=1, ls='--', c='k')	
	
	ax.set_aspect('equal', adjustable='box')
	legend = ax.get_legend()
	legend.set_bbox_to_anchor((1, 1))
	legend.set_frame_on(False)
	plt.savefig('/home/boj924/AD_Tau_PTA/results/Tau-noTau_enrich.pdf', bbox_inches='tight' )
	plt.show()

if __name__ == "__main__":
        main()
#	snv()






import matplotlib.pyplot as plt
import matplotlib.font_manager
import seaborn as sns
plt.rcParams['font.size'] = 10
plt.rcParams['axes.linewidth'] = 1
plt.rcParams['axes.spines.right'] = False
plt.rcParams['axes.spines.top'] = False
import numpy as np
import seaborn as sns
import seaborn.objects as so
np.set_printoptions(precision=8)
from scipy import stats
import statsmodels.api as sm
import statsmodels.formula.api as smf
import pandas as pd
import sys
import re


def indel_bin():

	csq_in    = '/home/boj924/AD_Tau_PTA/results/all_indel_list.annovar.exonic_variant_function'
	metafile_in = '/home/boj924/AD_Tau_PTA/metafiles/all_meta.tab'
	outdir    = '/home/boj924/AD_Tau_PTA/figures/'
	clinic_in  = '/home/boj924/AD_Tau_PTA/metafiles/all_clinic'

	csq = pd.read_csv(csq_in, sep='\t', header=None, names=['num','csq','gene','chr','start',\
		'end','ref','alt','chr2','start2','end2','ref2','alt2','score','filter','infor'])
	csq['sample'] = csq['infor'].apply(lambda row: row.split(',')[-1]) 

	metafile = pd.read_csv(metafile_in, header=0,dtype=str)
	clinic = pd.read_csv(clinic_in, header=0, sep='\t', dtype=str)
	metafile = pd.merge(metafile, clinic, left_on='donor', right_on='Case_ID')
	metafile['Age'] = metafile['Age'].astype(float)

	metafile = metafile.loc[(metafile['Age'] >=57 ) & (metafile['Age'] <=91)]

	csq_meta = pd.merge(csq, metafile, left_on='sample', right_on='sample')
	print(csq_meta.loc[csq_meta['group']=='ctrl'].count())

	csq_meta_summary = csq_meta.groupby(['sample','group'])['csq'].value_counts().to_frame(name='count')
	csq_meta_summary.reset_index(inplace=True)
	csq_meta_summary_sub = csq_meta_summary.loc[csq_meta_summary['csq']!='unknown']

	csq_meta_summary_refill = pd.DataFrame()
	for icsq, grp in csq_meta_summary_sub.groupby('csq'):
		merge_grp = pd.merge(grp, metafile, left_on='sample', right_on='sample', 
					how='right', suffixes=('_grp',''))
		merge_grp['csq'] = icsq	
		merge_grp['count'].fillna(0, inplace=True)
		if len(csq_meta_summary_refill) == 0:
			csq_meta_summary_refill = merge_grp
		else:
			csq_meta_summary_refill = pd.concat([csq_meta_summary_refill, merge_grp])

	csq_meta_summary_refill['bin_grp'] = csq_meta_summary_refill['group'].apply(lambda x: 'ctrl' if x=='ctrl' else 'AD')
	print(csq_meta_summary_refill.loc[csq_meta_summary_refill['group']=='ctrl'])

	for icsq, grp in csq_meta_summary_refill.groupby('csq'):
		model = smf.ols('count ~ group', data=grp).fit()
		if len(grp) < 4: continue	
		print(icsq)
		print(model.fvalue, model.f_pvalue)

	for icsq, grp in csq_meta_summary_refill.groupby('csq'):
		model = smf.ols('count ~ bin_grp', data=grp).fit()
		if len(grp) < 4: continue	
		print(icsq)
		print(model.fvalue, model.f_pvalue)

	color_palette= { 'AD': 'tab:orange', 'ctrl': 'tab:blue', 'Tau':'tab:red', 'noTau':'tab:pink' }	
	bin_color_palette= { 'AD': 'magenta', 'ctrl': 'tab:blue'}	

	label = { 'AD': 'AD Tau agnostic', 'ctrl': 'Control', 'Tau':'P-tau+', 'noTau':'P-tau-' }

	sns.set_style("whitegrid")	
	(
	so.Plot(csq_meta_summary_refill, y='count', x="csq", group='group', color='group')
	.add(so.Dot(),so.Agg(), so.Dodge(), )
	.scale(color=color_palette)
	.add(so.Range(),so.Est(errorbar=("se",2)), so.Dodge())
	.add(so.Range(),so.Est(errorbar=("se",2)), so.Dodge())
	.save(outdir+'/AD_indel_csq.pdf',format='pdf', dpi=300, bbox_inches='tight')
	)

	sns.set_style("whitegrid")	
	(
	so.Plot(csq_meta_summary_refill, y='count', x="csq", group='bin_grp', color='bin_grp')
	.add(so.Dot(),so.Agg(), so.Dodge(), )
	.scale(color=bin_color_palette)
	.add(so.Range(),so.Est(errorbar=("se",2)), so.Dodge())
	.add(so.Range(),so.Est(errorbar=("se",2)), so.Dodge())
	.save(outdir+'/AD_indel_csq_binary.pdf',format='pdf', dpi=300, bbox_inches='tight')
	)


	##only sig
	csq_meta_summary_refill_sub = csq_meta_summary_refill.loc[csq_meta_summary_refill['csq'].isin(['frameshift deletion','nonframeshift deletion'])]
	for icsq, grp in csq_meta_summary_refill_sub.groupby('csq'):
		ctrl = grp.loc[grp['group'] == 'ctrl', 'count']
		AD = grp.loc[grp['group'] == 'AD','count']
		Tau = grp.loc[grp['group'] == 'Tau', 'count']
		noTau = grp.loc[grp['group'] == 'noTau', 'count']
		print(stats.ttest_ind(ctrl, AD))
		print(stats.ttest_ind(ctrl, Tau))
		print(stats.ttest_ind(ctrl, noTau))
		print(stats.ttest_ind(Tau, noTau))

	sns.set_style("whitegrid")	
	(
	so.Plot(csq_meta_summary_refill_sub, y='count', x="csq", group='group', color='group')
	.add(so.Dot(),so.Agg(), so.Dodge(), )
	.scale(color=color_palette)
	.add(so.Range(),so.Est(errorbar=("se",2)), so.Dodge())
	.add(so.Range(),so.Est(errorbar=("se",2)), so.Dodge())
	.save(outdir+'/AD_indel_csq_sig.pdf',format='pdf', dpi=300, bbox_inches='tight')
	)


def snv_bin():

	csq_in    = '/home/boj924/AD_Tau_PTA/results/all_snv_list.annovar.exonic_variant_function'
	metafile_in = '/home/boj924/AD_Tau_PTA/metafiles/all_meta.tab'
	clinic_in  = '/home/boj924/AD_Tau_PTA/metafiles/all_clinic'
	outdir    = '/home/boj924/AD_Tau_PTA/figures/'

	csq = pd.read_csv(csq_in, sep='\t', header=None, names=['num','csq','gene','infor'])
	csq['sample'] = csq['infor'].apply(lambda row: row.split(' ')[-1]) 

	metafile = pd.read_csv(metafile_in, header=0,dtype=str)
	clinic = pd.read_csv(clinic_in, header=0, sep='\t', dtype=str)
	metafile = pd.merge(metafile, clinic, left_on='donor', right_on='Case_ID')
	metafile['Age'] = metafile['Age'].astype(float)
#	metafile = metafile.loc[(metafile['Age'] >=57 ) & (metafile['Age'] <=91)]
	csq_meta = pd.merge(csq, metafile, left_on='sample', right_on='sample')

	csq_meta_summary = csq_meta.groupby(['sample','group'])['csq'].value_counts().to_frame(name='count')
	csq_meta_summary.reset_index(inplace=True)
	csq_meta_summary_sub = csq_meta_summary.loc[csq_meta_summary['csq']!='unknown']

	csq_meta_summary_refill = pd.DataFrame()
	for icsq, grp in csq_meta_summary_sub.groupby('csq'):
		merge_grp = pd.merge(grp, metafile, left_on='sample', right_on='sample', 
					how='right', suffixes=('_grp',''))
		merge_grp['csq'] = icsq	
		merge_grp['count'].fillna(0, inplace=True)
		if len(csq_meta_summary_refill) == 0:
			csq_meta_summary_refill = merge_grp
		else:
			csq_meta_summary_refill = pd.concat([csq_meta_summary_refill, merge_grp])

	csq_meta_summary_refill['bin_grp'] = csq_meta_summary_refill['group'].apply(lambda x: 'ctrl' if x=='ctrl' else 'AD')

	for icsq, grp in csq_meta_summary_sub.groupby('csq'):
		model = smf.ols('count ~ group', data=grp).fit()
		if len(grp) < 4: continue
		print(icsq)
		print(model.fvalue, model.f_pvalue)

	for icsq, grp in csq_meta_summary_refill.groupby('csq'):
		model = smf.ols('count ~ bin_grp', data=grp).fit()
		if len(grp) < 4: continue	
		print(icsq)
		print(model.fvalue, model.f_pvalue)

	csq_meta_summary_refill['bin_grp'] = csq_meta_summary_refill['group'].apply(lambda x: 'ctrl' if x=='ctrl' else 'AD')


	color_palette= { 'AD': 'tab:orange', 'ctrl': 'tab:blue', 'Tau':'tab:red', 'noTau':'tab:pink' }	
	bin_color_palette= { 'AD': 'magenta', 'ctrl': 'tab:blue'}
	label = { 'AD': 'AD Tau agnostic', 'ctrl': 'Control', 'Tau':'P-tau+', 'noTau':'P-tau-' }

	sns.set_style("whitegrid")	
	(
	so.Plot(csq_meta_summary_refill, y='count', x="csq", group='group', color='group')
	.add(so.Dot(),so.Agg(), so.Dodge(), )
	.scale(color=color_palette)
	.add(so.Range(),so.Est(errorbar=("se",2)), so.Dodge())
	.add(so.Range(),so.Est(errorbar=("se",2)), so.Dodge())
	.save(outdir+'/AD_snv_csq.pdf',format='pdf', dpi=300, bbox_inches='tight')
	)

	sns.set_style("whitegrid")	
	(
	so.Plot(csq_meta_summary_refill, y='count', x="csq", group='bin_grp', color='bin_grp')
	.add(so.Dot(),so.Agg(), so.Dodge(), )
	.scale(color=bin_color_palette)
	.add(so.Range(),so.Est(errorbar=("se",2)), so.Dodge())
	.add(so.Range(),so.Est(errorbar=("se",2)), so.Dodge())
	.save(outdir+'/AD_snv_csq_binary.pdf',format='pdf', dpi=300, bbox_inches='tight')
	)


if __name__ == "__main__":
#	snv_bin()
	indel_bin()




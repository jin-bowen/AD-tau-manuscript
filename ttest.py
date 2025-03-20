import matplotlib.pyplot as plt
import matplotlib.font_manager
plt.rcParams['font.size'] = 10
plt.rcParams['axes.linewidth'] = 2
plt.style.use('seaborn-v0_8-talk')
plt.rcParams['axes.spines.right'] = False
plt.rcParams['axes.spines.top'] = False
import numpy as np
from scipy import stats
from matplotlib.lines import Line2D
import statsmodels.api as sm
import statsmodels.formula.api as smf
import pandas as pd
import sys
import re

def lmer(df):

	md = smf.mixedlm("burden ~  C(group) * Age", df, groups=df["donor"],			)
	mdf = md.fit()
	print(mdf.summary())

def var():

	burden_in   = sys.argv[1]
	metafile_in = sys.argv[2]
	clinic_in   = sys.argv[3]
	out_name    = sys.argv[4]

	color_palette = 'Pastel1'
	burden = pd.read_csv(burden_in, sep=',', header=0, index_col=0)
	burden.dropna(inplace=True)

	metafile = pd.read_csv(metafile_in, header=0)
	clinic = pd.read_csv(clinic_in, header=0, sep='\t')
	metafile = pd.merge(metafile,clinic, left_on='donor',right_on='Case_ID')
	metafile['Age'] = metafile['Age_(yrs)']

	burden_meta = pd.merge(burden,metafile,left_index=True, right_on='sample')

	color_palette= { 'AD': 'tab:orange', 'ctrl': 'tab:blue', 'Tau':'tab:red', 'noTau':'tab:pink' }	
	cols=['snv','indel']

	burden_meta.to_csv('burden_meta.res')

	burden_meta_grp = burden_meta.groupby('group')
	for col in cols:
		max_value = 0
		fig,ax = plt.subplots(constrained_layout=True, figsize=(10,5), ncols=2,\
				frameon = False, sharex=True, sharey=True)
		for igrp in ['AD','ctrl']:
			df = burden_meta_grp.get_group(igrp)
			df_sub = df.loc[df['muttype']==col]	
			scale_burden = df_sub["burden"]
			if igrp == 'ctrl':
				xpos=np.arange(0, 110)
				res = stats.linregress(df_sub['Age'], scale_burden)
				ax[0].plot(xpos, res.slope * xpos + res.intercept, c='tab:blue', lw=1)
	
			ax[0].scatter(x=df_sub['Age'], y=scale_burden, c=color_palette[igrp], edgecolors='black', \
					label=igrp, linewidth=2)
			max_value = max(max_value, max(scale_burden))
		ax[0].set_xlabel('Age')
		ax[0].set_ylabel('Burden')

		for igrp in ['Tau','noTau','ctrl']:
			df = burden_meta_grp.get_group(igrp)
			df_sub = df.loc[df['muttype']==col]	
			scale_burden = df_sub["burden"]
			if igrp == 'ctrl':
				xpos=np.arange(0, 110)
				res = stats.linregress(df_sub['Age'], scale_burden)
				ax[1].plot(xpos, res.slope * xpos + res.intercept, c='tab:blue', lw=1)
	
			ax[1].scatter(x=df_sub['Age'], y=scale_burden, c=color_palette[igrp], edgecolors='black', \
					label=igrp, linewidth=2)
			max_value = max(max_value, max(scale_burden))

		ax[1].set_xlabel('Age')
		ax[1].set_ylabel('Burden')
		ax[1].set_xlabel('Age')
		plt.title(col)
		plt.ylim([-1,5000])
		plt.legend(loc='center left',bbox_to_anchor=(1.01, 0.5), frameon=False)
		plt.savefig(out_name+'_'+col+'.pdf', dpi=200,  bbox_inches='tight')

def main():

	burden_in   = sys.argv[1]
	metafile_in = sys.argv[2]
	clinic_in   = sys.argv[3]
	out_name    = sys.argv[4]

	muttype='snv'
	color_palette = 'Pastel1'
	burden = pd.read_csv(burden_in, sep=',', header=0, index_col=0)
	burden.dropna(inplace=True)

	metafile = pd.read_csv(metafile_in, header=0)
	clinic = pd.read_csv(clinic_in, header=0, sep='\t')
	metafile = pd.merge(metafile,clinic, left_on='donor',right_on='Case_ID')

	burden_meta = pd.merge(burden,metafile,left_index=True, right_on='sample')

	color_palette= { 'AD': 'tab:orange', 'ctrl': 'tab:blue', 'Tau':'tab:red', 'noTau':'tab:pink' }	
	cols=['snv','indel']

	for col in cols:
		max_value = 0
		fig,ax = plt.subplots(constrained_layout=True,figsize=(8,5), frameon = False)
		burden_meta_sub = burden_meta.loc[burden_meta['muttype']==col]
		for igrp, df in burden_meta_sub.groupby('group'):
		
			scale_burden = df["burden"]
			if igrp == 'ctrl':
				xpos=np.arange(0, 110)
				res = stats.linregress(df['Age'], scale_burden)
				ax.plot(xpos, res.slope * xpos + res.intercept, c='tab:blue', lw=1)
	
			ax.scatter(x=df['Age_(yrs)'], y=scale_burden, c=color_palette[igrp], edgecolors='black', \
					label=igrp, linewidth=2)
			max_value = max(max_value, max(scale_burden))
	
		ax.set_xlabel('Age')
		ax.set_ylabel('Burden')
		ax.set_title(col)
		ax.set_ylim([-1,max_value*1.3])
		ax.set_ylim([-1,5000])

		plt.legend(loc='center left',bbox_to_anchor=(1.01, 0.5), frameon=False)
		plt.savefig(out_name+'_'+col+'.pdf', dpi=200,  bbox_inches='tight')
		plt.show()

if __name__ == "__main__":
	var()



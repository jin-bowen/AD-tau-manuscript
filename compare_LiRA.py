import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib as mpl
import matplotlib.font_manager
plt.rcParams['font.size'] = 10
plt.rcParams['axes.linewidth'] = 2
plt.style.use('seaborn-v0_8-talk')
plt.rcParams['axes.spines.right'] = False
plt.rcParams['axes.spines.top'] = False
import numpy as np
from scipy import stats
import pandas as pd
import sys
import re

def figure2(df, out_name):

	all_burden_type = ['sim.negbin.burden','mean.burden','poisson.burden','negbin.burden']
	df = df.loc[(df['mutation.type']=='indel') & (df['Cell_ID']!='ALZ1647BA9-C')]

	fig,ax = plt.subplots(constrained_layout=True,figsize=(len(all_burden_type)*3,3), frameon = False, ncols=len(all_burden_type))
	for i, burden_type in enumerate(all_burden_type):

		slope, intercept, r_value, p_value, std_err = stats.linregress(df["burden.genome"],\
			df[burden_type])
	
		maxval = max( df[burden_type])
		maxval = max( maxval, max(df["burden.genome"]))
		maxval = maxval*1.3
	
		xpos=np.arange(0, maxval)
		ypos=np.arange(0, maxval)

		print(df[burden_type])
		print(df["burden.genome"])
	
		ax[i].plot(xpos, ypos, c='gray', lw=2)
		ax[i].scatter(x=df[burden_type], y=df["burden.genome"], c='tab:green', edgecolors='black', linewidth=2)
		ax[i].set_xlabel(burden_type)
		ax[i].set_ylabel('Original burden')
		ax[i].set_title('corr=%s'%round(r_value,3), loc='right')
	plt.savefig(out_name+'.png', dpi=200,  bbox_inches='tight' )
#	plt.show()

def figure1(df, out_name):

	all_burden_type = ['burden.genome','sim.negbin.burden','mean.burden','poisson.burden','negbin.burden']
	
	fig,ax = plt.subplots(constrained_layout=True,figsize=(len(all_burden_type)*3,3), 
			sharey=True, sharex=True, frameon = False, ncols=len(all_burden_type))
	for i, burden_type in enumerate(all_burden_type):

		slope, intercept, r_value, p_value, std_err = stats.linregress(df["Somatic_burden"],\
			df[burden_type])
	
		maxval = max( df[burden_type])
		maxval = max( maxval, max(df["Somatic_burden"]))
		maxval = maxval*1.3
	
		xpos=np.arange(0, maxval)
		ypos=np.arange(0, maxval)

	
		ax[i].plot(xpos, ypos, c='gray', lw=2)
		ax[i].scatter(x=df[burden_type], y=df["Somatic_burden"], c='tab:green', edgecolors='black', linewidth=2)
		ax[i].set_xlabel(burden_type)
		ax[i].set_ylabel('LiRA call')
		ax[i].set_title('corr=%s'%round(r_value,3), loc='right')

	plt.xlim([0,5000])
	plt.ylim([0,5000])
	plt.savefig(out_name+'.png', dpi=200,  bbox_inches='tight' )
#	plt.show()

def figure(df, out_name):

	burden_type = 'negbin.burden'
	burden_type = 'burden'
	
	fig,ax = plt.subplots(constrained_layout=True,figsize=(3,3),frameon = False)
		
	slope, intercept, r_value, p_value, std_err = stats.linregress(df["Somatic_burden"],\
		df[burden_type])

	maxval = max( df[burden_type])
	maxval = max( maxval, max(df["Somatic_burden"]))
	maxval = maxval*1.1

	xpos=np.arange(0, maxval)
	ypos=np.arange(0, maxval)

	ax.plot(xpos, ypos, c='gray', lw=2)
	ax.scatter(x=df[burden_type], y=df["Somatic_burden"], c='tab:green', edgecolors='black', linewidth=2)
	ax.set_xlabel(burden_type)
	ax.set_ylabel('LiRA call')
	ax.set_title('corr=%s'%round(r_value,3), loc='right')

	plt.xlim([0,maxval])
	plt.ylim([0,maxval])
	plt.savefig(out_name+'.pdf', dpi=200,  bbox_inches='tight' )
	plt.show()

def main():

#	infile = '/home/boj924/AD_Tau_PTA/results/tau_AD_ctrl_res.tab' 
#	out_name='ctrl'

	infile = '/n/scratch/users/b/boj924/AD_hg19/summary.combined.txt'
	out_name='ADag'

#	infile = '/n/data1/bwh/pathology/miller/lab/Tau_hg19/summary.combined.txt' 
#	out_name='Tau'

	metafile = '/home/boj924/AD_Tau_PTA/AD_case_config/scan2_summary_meta.csv'
	metafile = '/home/boj924/AD_Tau_PTA/metafiles/ADnCtrl_metafile_august'
#	metafile = '/home/boj924/AD_Tau_PTA/metafiles/control_LiRA'

	color_palette = 'Pastel1'

	df = pd.read_csv(infile, sep='\t', header=0)
#	df_snv = df
#	df_snv = df.loc[df['muttype']=='snv']
	df_snv = df.loc[df['mutation.type']=='snv']

	meta_file = pd.read_csv(metafile, sep=',', header=0)
	meta_file['Somatic_burden'] = meta_file['Somatic_rate'] * 5.845001134
	print(meta_file)
	print(df_snv)
	df_merge = pd.merge(meta_file,df_snv,on='Cell_ID')

	figure(df_merge, out_name+'LiRA_Corr')
#	figure1(df_merge, out_name+'LiRA_Corr')
#	figure2(df, out_name+'_compare')

if __name__ == "__main__":
	main()


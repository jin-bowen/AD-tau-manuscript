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
from scipy.interpolate import UnivariateSpline,CubicSpline
from matplotlib.lines import Line2D
from sklearn.metrics import r2_score
import pandas as pd
import sys
import re

def realign():

	infile = sys.argv[1]
	out_name = sys.argv[2]

	color_palette = 'Pastel1'
	df = pd.read_csv(infile, sep=',', header=0)

	for spl in ['Sensitivity', 'muts_per_haploid_Gb']:
		fig,ax = plt.subplots(figsize=(6,3),frameon = False)

		sns.relplot(data=df, x="realign", y=spl, hue="mutation_type",\
			kind="scatter")

		ax.set_title(spl, loc='right')
		ax.set_xlim([-1,2])
		plt.savefig(out_name+'_'+spl+'.png', dpi=200,  bbox_inches='tight' )
		plt.show()


def repeat_test():

	infile = sys.argv[1]
	out_name = sys.argv[2]

	color_palette = 'Pastel1'
	df = pd.read_csv(infile, sep=',', header=0)

	cols= ['somatic.sens','rate.per.gb','ncalls']
	df_var = df.groupby(['cell','muttype'])[cols].var()
	df_mean = df.groupby(['cell','muttype'])[cols].mean()
	df_fano = (df_var.values) / (df_mean.values)
	print(df_fano)
	for spl in cols:
		fig,ax = plt.subplots(figsize=(3,3),frameon = False)

		sns.relplot(data=df, x="cell", y=spl, hue="muttype",\
			kind="scatter")

		ax.set_title(spl, loc='right')
		ax.set_xlim([-1,2])
		plt.savefig(out_name+'_'+spl+'.png', dpi=200,  bbox_inches='tight' )
		plt.show()

def sensitivity(infile,infile2,out_name):

	color_palette = 'Pastel1'
	df = pd.read_csv(infile, sep=',', header=0, index_col=0)
	df2 = pd.read_csv(infile2, sep=',', header=0, index_col=0)
	
	for muttype in ['snv','indel']:

		df_sub = df.loc[df['muttype']==muttype]
		df2_sub = df2.loc[df2['muttype']==muttype]

		for spl in ['somatic.sens','ncalls']:
			fig,ax = plt.subplots(figsize=(4,3), constrained_layout=True)
			ax.scatter(x=np.zeros(len(df_sub)), y=df_sub[spl], label='AD',edgecolors='black', linewidth=2)
			ax.scatter(x=np.ones(len(df2_sub)), y=df2_sub[spl], label='control', edgecolors='black', linewidth=2)
	
			ax.set_title(spl, loc='right')
			ax.set_xlim([-1,2])
			ax.set_xticks([0,1])
			ax.set_xticklabels(['AD','control'])
			ax.set_ylabel(spl)
			plt.savefig(out_name+'_%s'%spl+'_%s'%muttype+'.png', dpi=200,  bbox_inches='tight' )
			plt.show()

if __name__ == "__main__":
	globals()[sys.argv[1]](sys.argv[2],sys.argv[3],sys.argv[4])



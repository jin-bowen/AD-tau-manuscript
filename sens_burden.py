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

	burden_in   = sys.argv[1]
	metafile_in = sys.argv[2]
	clinic_in   = sys.argv[3]
	out_name    = sys.argv[4]
	muttype = sys.argv[5]
	
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

	fig,ax = plt.subplots(constrained_layout=True,figsize=(8,4),ncols=2,
			 frameon = False, sharex=True)
	df_sub=burden_meta	
	color_palette= { 'AD': 'tab:orange', 'ctrl': 'tab:blue', 'Tau':'tab:red', 'noTau':'tab:pink' }
	color_list = [ color_palette[x]  for x in df_sub['group'].values]	
	ax[0].scatter(x=df_sub['somatic.sens'], y=df_sub['burden'], edgecolor='k', c=color_list)
	ax[1].scatter(x=df_sub['somatic.sens'], y=df_sub['ncalls'], edgecolor='k', c=color_list)

	ax[1].set_xlim([0.05, 0.25])
	ax[0].set_xlabel('Sensitivity')
	ax[0].set_xlabel('Sensitivity')

	ax[0].set_ylabel('Burden')
	ax[1].set_ylabel('Raw call')

	plt.legend(loc='center left',bbox_to_anchor=(1.01, 0.5), frameon=False)
	plt.savefig(out_name+'_'+muttype+'-'+'.png', dpi=200,  bbox_inches='tight')
	plt.show()	
	plt.close()

if __name__ == "__main__":
	main()



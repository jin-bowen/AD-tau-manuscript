import matplotlib.pyplot as plt
import matplotlib.font_manager
plt.rcParams['font.size'] = 10
plt.rcParams['axes.linewidth'] = 1
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


def main():

	burden_in   = '/home/boj924/AD_Tau_PTA/results/tau_AD_ctrl_res.tab'
	sig_in      = '/home/boj924/AD_Tau_PTA/results/tauADCtrl_aging_contribute.tab'
	metafile_in = '/home/boj924/AD_Tau_PTA/metafiles/all_meta.tab'
	clinic_in   = '/home/boj924/AD_Tau_PTA/metafiles/all_clinic'
		
	burden = pd.read_csv(burden_in, sep=',', header=0,index_col=0)
#	burden = burden.loc[burden['mutation_type']==mutation_type]
	burden.dropna(inplace=True)
	burden['sample'] = burden.index

	sig_tab = pd.read_csv(sig_in, sep=',', header=0, index_col=0)
	sig_tab = sig_tab.astype(float)
	sig_tab = sig_tab/sig_tab.sum(axis=1).values.reshape((-1,1))
	burden_sig = pd.merge(burden, sig_tab, left_index=True, right_index=True)

	cols=['Signature_A','Signature_C']
	scaled_cols=['scaled_' + x for x in cols]
	for i, item in enumerate(cols):
		burden_sig[scaled_cols[i]] = burden_sig.apply(lambda row: row['burden'] * row[item], axis=1 )

	metafile = pd.read_csv(metafile_in, header=0,dtype=str)
	clinic = pd.read_csv(clinic_in, header=0, sep='\t',dtype=str)
	metafile = pd.merge(metafile,clinic, left_on='donor',right_on='Case_ID')
	metafile['Age'] = metafile['Age'].astype(float)
	burden_meta = pd.merge(burden_sig,metafile, on='sample', suffixes=('','_meta'))
	burden_meta = burden_meta.loc[burden_meta['sample']!='1995P_201001E3']

	burden_meta_snv = burden_meta.loc[burden_meta['mutation_type']=='snv']
	burden_meta_indel = burden_meta.loc[burden_meta['mutation_type']=='indel']

	burden_meta_reform = pd.merge(burden_meta_snv, burden_meta_indel, left_on='sample', right_on='sample', suffixes=('','_indel'))
	burden_meta_reform.to_csv('data_table.csv')
	print(burden_meta_reform)


if __name__ == "__main__":
	main()




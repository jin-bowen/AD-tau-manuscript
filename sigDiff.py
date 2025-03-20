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
from scipy.optimize import leastsq
from sklearn.metrics.pairwise import pairwise_distances
import pandas as pd
import sys
import re

def cos_sim(a,b):

	a=a.reshape(-1)
	b=b.reshape(-1)
	return np.dot(a,b)/(np.linalg.norm(a)*np.linalg.norm(b))

def sigDiff_indel(mat_file):
	
	clinic_file = "/home/boj924/AD_Tau_PTA/metafiles/all_meta.tab"
	burden_file = "/home/boj924/AD_Tau_PTA/results/burden_meta.res"
	
	mut = pd.read_csv(mat_file, sep=',', header=0, index_col=0)
	mut = mut/mut.sum(axis=1).values.reshape((-1,1))
	mut = mut.T
	sbs96 = mut.columns
	sbs96_shift = [ ele+'_shift' for ele in sbs96]

	clinic = pd.read_csv(clinic_file,header=0)
	donor = clinic.loc[clinic['group']=='Tau','donor'].unique()
	clinic_sub = clinic.loc[(clinic['donor'].isin(donor)) & \
				(clinic['group']!='AD') &  \
				(clinic['group']!='ctrl')]

	muttype='indel'
	burden = pd.read_csv(burden_file,header=0, index_col=0)
	burden_df = burden.loc[burden['muttype']==muttype,['sample','burden','Age']]
	
	df = pd.merge(mut,clinic_sub, left_index=True, right_on='sample')
	df = pd.merge(df, burden_df, on='sample')
	df[sbs96] = df[sbs96]/df[sbs96].sum(axis=1).values.reshape((-1,1))	
	df[sbs96] = df[sbs96] * df['burden'].values.reshape((-1,1))
	df_reform = df.groupby(['donor','group'])[sbs96].mean().reset_index()
	df_reform[sbs96_shift] = df_reform.groupby(['donor'])[sbs96].diff(axis=0,periods=-1).fillna(0)
	sub_df_reform = df_reform.loc[df_reform['group']=='Tau',sbs96_shift+['donor']]
	sub_df_reform.set_index('donor', inplace=True)
	sub_df_reform = sub_df_reform.T
	sub_df_reform.index = sbs96
	sub_df_reform.to_csv('tau_indel-diff.mat')
	sub_df_reform_mat = sub_df_reform.values

#	cosmic_raw = pd.read_csv('/home/boj924/AD_Tau_PTA/custome_sig/mutagen_id83.tab', sep=',', header=0, index_col=0)	
#	ID_list=cosmic_raw.columns[1:-1].to_list()
#	indel_transfer_tab = pd.read_csv('/home/boj924/AD_Tau_PTA/custome_sig/indel83_lookup.tab',header=0)
#	cosmic1 = pd.merge(cosmic_raw, indel_transfer_tab, left_index=True, right_on='indeltype1')
#	cosmic2 = pd.merge(cosmic1, sub_df_reform, left_on='indeltype2', right_index=True, how='right')
#	cosmic2.set_index('indeltype2',inplace=True)
#	cosmic = cosmic2[ID_list]
#	cosmic_mat = cosmic.values
#	cosmic.to_csv('/home/boj924/AD_Tau_PTA/custome_sig/mutagen_id83_v2.tab')
#	cosmic_raw = pd.read_csv('/home/boj924/AD_Tau_PTA/custome_sig/COSMIC_v3.4_ID_GRCh37.txt', sep='\t', header=0, index_col=0)	
#	ID_list=cosmic_raw.columns[1:-1].to_list()
#	indel_transfer_tab = pd.read_csv('/home/boj924/AD_Tau_PTA/custome_sig/indel83_lookup.tab',header=0)
#	cosmic1 = pd.merge(cosmic_raw, indel_transfer_tab, left_index=True, right_on='indeltype1')
#	cosmic2 = pd.merge(cosmic1, sub_df_reform, left_on='indeltype2', right_index=True, how='right')
#	cosmic2.set_index('indeltype2',inplace=True)
#	cosmic = cosmic2[ID_list]
#	cosmic_mat = cosmic.values
#	cosmic.to_csv('/home/boj924/AD_Tau_PTA/custome_sig/COSMIC_v3.4_ID_GRCh37_v2.txt')

	cosmic = pd.read_csv('/home/boj924/AD_Tau_PTA/custome_sig/COSMIC_v3.4_ID_mutagen_v2.txt',header=0, index_col=0)	
	cosmic_mat = cosmic.values
	print(cosmic)

	coef_res = np.zeros(shape=(len(cosmic.columns),len(sub_df_reform.columns)))
	cos_sim_score_res = np.zeros(shape=(len(sub_df_reform.columns)))
	for i in range(len(sub_df_reform.columns)):
		x0 = np.ones(shape=(len(cosmic.columns),1))
		def fun(x): return np.matmul(cosmic_mat,x) - sub_df_reform_mat[:,i]

		res, los = leastsq(fun, x0=x0)
		reconstruct = np.matmul(cosmic_mat, res.reshape((-1,1)))
		cos_sim_score = cos_sim(reconstruct,sub_df_reform_mat[:,i])

		coef_res[:,i] = res
		cos_sim_score_res[i] = cos_sim_score	


	order = [1353,1647,2208,1456,2207,1995,1828]
	coef_res_df = pd.DataFrame(data=coef_res.T, columns=cosmic.columns, index=sub_df_reform.columns)

	sns.clustermap(coef_res_df.loc[order],xticklabels=True, cmap='coolwarm', vmin=-2000, vmax=2000)
	plt.show()
	plt.close()

	label_list = [str(ele) for ele in sub_df_reform.columns]
	fig, ax = plt.subplots(figsize=(5,3))
	ax.bar(label_list, cos_sim_score_res)
	ax.set_ylabel('Cosine similarity')
	plt.show()


def sigSum(mat_file):
	
	clinic_file = "/home/boj924/AD_Tau_PTA/metafiles/all_meta.tab"
	burden_file = "/home/boj924/AD_Tau_PTA/results/burden_meta.res"
	
	mut = pd.read_csv(mat_file, sep=',', header=0, index_col=0)
	mut = mut - 0.0001	
	mut = mut/mut.sum(axis=1).values.reshape((-1,1))
	sbs96 = mut.columns.tolist()
	sbs96_shift = [ ele+'_shift' for ele in sbs96]

	clinic = pd.read_csv(clinic_file,header=0)
	donor = clinic.loc[clinic['group']=='Tau','donor'].unique()
	clinic_sub = clinic.loc[(clinic['donor'].isin(donor)) & (clinic['group']!='AD')]
#	clinic_sub = clinic.loc[(clinic['donor'].isin(donor)) & (clinic['group']!='noTau')]

	muttype='snv'
	burden = pd.read_csv(burden_file,header=0, index_col=0)
	burden_df = burden.loc[burden['muttype']==muttype,['sample','burden','Age']]
	
	df = pd.merge(mut,clinic_sub, left_index=True, right_on='sample')
	df = pd.merge(df, burden_df, on='sample')
	df[sbs96] = df[sbs96]/df[sbs96].sum(axis=1).values.reshape((-1,1))	
	df[sbs96] = df[sbs96] * df['burden'].values.reshape((-1,1))
	df_reform = df.groupby(['donor','group'])[sbs96].mean().reset_index()

	Tau_df_reform = df_reform.loc[df_reform['group']=='Tau',sbs96+['donor']]
	noTau_df_reform = df_reform.loc[df_reform['group']=='noTau',sbs96+['donor']]
	print(Tau_df_reform)
	Tau_df_reform.set_index('donor', inplace=True)
	Tau_df_reform = Tau_df_reform.T
	Tau_df_reform.to_csv('tau_snv_mean.mat')

	noTau_df_reform.set_index('donor', inplace=True)
	noTau_df_reform = noTau_df_reform.T
	noTau_df_reform.to_csv('notau_snv_mean.mat')

def sigDiff(mat_file):
	
	clinic_file = "/home/boj924/AD_Tau_PTA/metafiles/all_meta.tab"
	burden_file = "/home/boj924/AD_Tau_PTA/results/burden_meta.res"
	
	mut = pd.read_csv(mat_file, sep=',', header=0, index_col=0)
	mut = mut - 0.0001	
	mut = mut/mut.sum(axis=1).values.reshape((-1,1))
	sbs96 = mut.columns
	sbs96_shift = [ ele+'_shift' for ele in sbs96]

	clinic = pd.read_csv(clinic_file,header=0)
	donor = clinic.loc[clinic['group']=='Tau','donor'].unique()
#	clinic_sub = clinic.loc[(clinic['donor'].isin(donor)) & (clinic['group']!='AD')]
	clinic_sub = clinic.loc[(clinic['donor'].isin(donor)) & (clinic['group']!='noTau')]

	muttype='snv'
	burden = pd.read_csv(burden_file,header=0, index_col=0)
	burden_df = burden.loc[burden['muttype']==muttype,['sample','burden','Age']]
	
	df = pd.merge(mut,clinic_sub, left_index=True, right_on='sample')
	df = pd.merge(df, burden_df, on='sample')
	df[sbs96] = df[sbs96]/df[sbs96].sum(axis=1).values.reshape((-1,1))	
	df[sbs96] = df[sbs96] * df['burden'].values.reshape((-1,1))
	df_reform = df.groupby(['donor','group'])[sbs96].mean().reset_index()
	df_reform[sbs96_shift] = df_reform.groupby(['donor'])[sbs96].diff(axis=0,periods=-1).fillna(0)

	#sub_df_reform = df_reform.loc[df_reform['group']=='Tau',sbs96_shift+['donor']]
	sub_df_reform = df_reform.loc[df_reform['group']=='AD',sbs96_shift+['donor']]

	sub_df_reform.set_index('donor', inplace=True)
	sub_df_reform = sub_df_reform.T
	sub_df_reform.index = sbs96
	sub_df_reform.to_csv('tauAD_diff.mat')
	print(sub_df_reform)
#	sub_df_reform.to_csv('tau_diff.mat')

	sub_df_reform_mat = sub_df_reform.values
	sub_df_reform_mat = sub_df_reform_mat * -1	

	cosmic = pd.read_csv('/home/boj924/AD_Tau_PTA/custome_sig/COSMIC_v3.4_SBS_GRCh37.txt', sep='\t', header=0, index_col=0)	
#	cosmic = pd.read_csv('/home/boj924/AD_Tau_PTA/custome_sig/Mutagen53_COSMIC.tsv', sep='\t', header=0, index_col=0)
#	cosmic = pd.read_csv('/home/boj924/AD_Tau_PTA/custome_sig/Mutagen53_sigA.tsv', sep='\t', header=0, index_col=0)	
	print(cosmic.shape)
	cosmic_mat = cosmic.values

	coef_res = np.zeros(shape=(len(cosmic.columns),len(sub_df_reform.columns)))
	cos_sim_score_res = np.zeros(shape=(len(sub_df_reform.columns)))
	for i in range(len(sub_df_reform.columns)):
		x0 = np.ones(shape=(len(cosmic.columns),1))
		def fun(x): return np.matmul(cosmic_mat,x) - sub_df_reform_mat[:,i]

		res, los = leastsq(fun, x0=x0)
		reconstruct = np.matmul(cosmic_mat, res.reshape((-1,1)))
		cos_sim_score = cos_sim(reconstruct,sub_df_reform_mat[:,i])

		coef_res[:,i] = res
		cos_sim_score_res[i] = cos_sim_score	


	order = [1353,1647,2208,1456,2207,1995,1828]
	coef_res_df = pd.DataFrame(data=coef_res.T, columns=cosmic.columns, index=sub_df_reform.columns)
	print(cos_sim_score_res)
#	sns.heatmap(coef_res_df.loc[order], xticklabels=True, cmap='coolwarm')
	sns.clustermap(coef_res_df.loc[order],xticklabels=True, cmap='coolwarm')
#	fig, ax = plt.subplots()
#	ax.imshow(coef_res_df)
#	ax.set_xticks(np.arange(len(cosmic.columns)), labels=cosmic.columns, )
#	ax.set_yticks(np.arange(len(sub_df_reform.columns)), labels=sub_df_reform.columns)
#	ax.set_xticklabels(ax.get_xticks(), rotation = 90)
	plt.show()

if __name__ == "__main__":

	mat_file = sys.argv[1]
	sigSum(mat_file)
#	sigDiff(mat_file)
	sigDiff_indel(mat_file)


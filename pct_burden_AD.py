import matplotlib.pyplot as plt
import matplotlib.font_manager
plt.rcParams['font.size'] = 10
plt.style.use('seaborn-v0_8-poster')
plt.rcParams['axes.linewidth'] = 2
plt.rcParams['axes.spines.right'] = False
plt.rcParams['axes.spines.top'] = False
import numpy as np
np.set_printoptions(precision=8)
from scipy import stats
from sklearn import preprocessing
import seaborn as sns
import statsmodels.api as sm
import statsmodels.formula.api as smf
import pandas as pd
import sys
import re

def norm_to_zero_one(df):
	return (df - df.min()) * 1.0 / (df.max() - df.min())

def print(*args):
	__builtins__.print(*("%.8f" % a if isinstance(a, float) else a for a in args))

def mda_snv():

	burden_in = '/home/boj924/AD_Tau_PTA/metafiles/AD_mda.csv'
	metafile_in = '/home/boj924/AD_Tau_PTA/metafiles/all_meta.tab'
	clinic_in   = '/home/boj924/AD_Tau_PTA/metafiles/all_clinic'
	AT8pct_in   = '/home/boj924/AD_Tau_PTA/metafiles/quant_sub.csv'
	out_name    = '/home/boj924/AD_Tau_PTA/results/AT8_burden_mda'

	burden = pd.read_csv(burden_in, sep=',', header=0)

	color_palette = 'Pastel1'
	clinic = pd.read_csv(clinic_in, header=0, sep='\t',dtype=str)

	AT8pct = pd.read_csv(AT8pct_in, header=0)
	AT8pct['donor'] = AT8pct['donor'].astype(str)

	burden['Case_ID'] = burden['Case_number'].astype(str)
	clinic['Case_ID'] = clinic['Case_ID'].astype(str)

	df = pd.merge(burden, clinic, on='Case_ID')
	df['Age'] = df['Age_(yrs)']
	df['Age'] = df['Age'].astype(float)
	df = df.loc[df['Tissue_origin']=='PFC']

	var = 'sSNVs'
	df_tau = pd.merge(df, AT8pct, left_on='Case_ID', right_on='donor')
	print(df_tau[['sSNVs','Age','IHC_neuron_pos_percent','Case_ID']])

	md = smf.mixedlm(" %s ~ Age + IHC_neuron_pos_percent "%var, df_tau, groups=df_tau["Case_ID"])
	mdf = md.fit()
	rescale_results = mdf.summary().tables[1]
	rescale_results['Coef.'] = rescale_results['Coef.'].astype(float)
	print(rescale_results)

	fig,ax = plt.subplots(constrained_layout=True,figsize=(3,4))
	g=sns.lmplot(x="IHC_neuron_pos_percent", y=var, data=df_tau, \
			legend=False,scatter_kws={'color':'tab:orange'}, line_kws={'color':'tab:orange'})
	g.set(xlim=(0,30))
	g.set(ylim=(-1000, 6000))
	plt.savefig(out_name+'_'+var+'_combine.pdf', dpi=300,  bbox_inches='tight')
	plt.show()
	plt.close()

def quant_snv():
	
	burden_in = '/home/boj924/AD_Tau_PTA/results/tau_AD_ctrl_res.tab'
	metafile_in = '/home/boj924/AD_Tau_PTA/metafiles/all_meta.tab'
	clinic_in   = '/home/boj924/AD_Tau_PTA/metafiles/all_clinic'
	AT8pct_in   = '/home/boj924/AD_Tau_PTA/metafiles/quant.csv'
	out_name    = '/home/boj924/AD_Tau_PTA/results/AT8_snv-AD'

	vartype = 'snv'
	burden = pd.read_csv(burden_in, sep=',', header=0,index_col=0)
	burden = burden.loc[(burden['mutation_type']==vartype) & (burden.index!='1995P_201001E3'),]
	color_palette = 'Pastel1'

	metafile = pd.read_csv(metafile_in, header=0,dtype=str)
	clinic = pd.read_csv(clinic_in, header=0, sep='\t',dtype=str)
	metafile = pd.merge(metafile,clinic, left_on='donor',right_on='Case_ID')
	metafile['Age'] = metafile['Age_(yrs)'].astype(float)

	AT8pct = pd.read_csv(AT8pct_in,sep=',',header=0, index_col=0)
	AT8pct.index = [str(x) for x in AT8pct.index]
	metainfor = pd.merge(metafile, AT8pct, left_on='Case_ID', right_index=True, how='left')
	df = pd.merge(burden, metainfor, left_index=True, right_on='sample')
	var = 'burden'
	df_tau = df.loc[df['group'].isin(['AD'])]

	df_tau[AT8pct.columns] = df_tau[AT8pct.columns].apply(norm_to_zero_one)

	res = []
	for col in AT8pct.columns:
		print(col)
		df_tau_sub = df_tau[['Age','group',var, col,"Case_ID"]].dropna()
		md = smf.mixedlm(" %s ~ Age + %s "%(var,col), df_tau_sub, groups=df_tau_sub["Case_ID"])
		mdf = md.fit()
		rescale_results = mdf.summary().tables[1]
		rescale_results['Coef.'] = rescale_results['Coef.'].astype(float)
		res.append(rescale_results.loc[col,['Coef.','[0.025','0.975]','P>|z|']].tolist())
	res_df = pd.DataFrame(res, columns=['coef','25','975','pval'], index=AT8pct.columns)
	res_df = res_df.astype(float)

	cols = ['sup_frontal_np','sup_frontal_dp','sup_frontal_neuron_loss','sup_frontal_nft_count','IHC_neuron_pos_percent']
	res_df = res_df.loc[cols]
	print(res_df)

	plt.figure(figsize=(7,2), dpi=300)
	ci = [res_df.iloc[::-1]['coef'] - res_df.iloc[::-1]['25'].values, res_df.iloc[::-1]['975'].values - res_df.iloc[::-1]['coef']]
	plt.errorbar(x=res_df.iloc[::-1]['coef'], y=res_df.iloc[::-1].index.values, xerr=ci,
		    color='black',  capsize=3, linestyle='None', linewidth=1,
		    marker="o", markersize=5, mfc="black", mec="black")
	plt.axvline(x=1, linewidth=0.8, linestyle='--', color='black')
	plt.tick_params(axis='both', which='major', labelsize=8)
	plt.xlabel('coef and 95% Confidence Interval', fontsize=8)
	plt.tight_layout()
	plt.savefig('%s_%s_forest_plot.pdf'%(out_name,var), dpi=200,  bbox_inches='tight')
	plt.show()


def scale_indel():
	
	tbp_count = '/home/boj924/AD_Tau_PTA/results/indel_spectrum_count.tab'
	ID4_count = '/home/boj924/AD_Tau_PTA/results/ID83_Musical_AD/Decompose_Solution/Activities/Decompose_Solution_Activities.txt'
	burden_in = '/home/boj924/AD_Tau_PTA/results/tau_AD_ctrl_res.tab'

	metafile_in = '/home/boj924/AD_Tau_PTA/metafiles/all_meta.tab'
	clinic_in   = '/home/boj924/AD_Tau_PTA/metafiles/all_clinic'
	AT8pct_in   = '/home/boj924/AD_Tau_PTA/metafiles/full_digpath_20241103.csv'
	depth_in    = '/home/boj924/AD_Tau_PTA/metafiles/AD_depth_quantile.txt'

	vartype ='indel'
	out_name    = '/home/boj924/AD_Tau_PTA/results/AT8_%s_burden'%vartype

	depth_profile = pd.read_csv(depth_in, header=0,index_col=0)
	depth_profile['class'] = depth_profile.apply(lambda row: 0 if row['25%'] <10 else 1, axis=1)

	burden = pd.read_csv(burden_in, sep=',', header=0,index_col=0)
	burden = burden.loc[(burden['mutation_type']==vartype) & (burden.index!='1995P_201001E3'),]
	burden.dropna(inplace=True)
	
	metafile = pd.read_csv(metafile_in, header=0,dtype=str)
	clinic = pd.read_csv(clinic_in, header=0, sep='\t',dtype=str)
	metafile = pd.merge(metafile,clinic, left_on='donor',right_on='Case_ID')
	metafile['Age'] = metafile['Age_(yrs)'].astype(float)

	AT8pct = pd.read_csv(AT8pct_in,header=0)
	AT8pct['Case_ID'] = AT8pct['ADRCID'].astype(str)
	metainfor = pd.merge(metafile, AT8pct, on='Case_ID', how='left')

	df = pd.merge(burden, metainfor, left_index=True, right_on='sample')
	var = 'burden'

	df_AD = df.loc[df['group']=='AD',]
	md = smf.mixedlm(" %s ~ Age + IHC_neuron_pos_percent "%var, df_AD, groups=df_AD["Case_ID"])#	mdf = md.fit()
	mdf = md.fit()
	rescale_results = mdf.summary().tables[1]
	rescale_results['Coef.'] = rescale_results['Coef.'].astype(float)

	rescale_var = 'resacle_' + var
	df_AD[rescale_var] = df_AD[var]  - (df_AD['Age'] * rescale_results.loc['Age','Coef.'] + rescale_results.loc['Intercept','Coef.'])

	fig,ax = plt.subplots(constrained_layout=True,figsize=(4,5))
	sns.regplot(x="IHC_neuron_pos_percent", y=rescale_var, data=df_AD, color="tab:orange", \
			ax=ax,  scatter_kws={'linewidths':1, 'edgecolor':'k'});
	ax.set_xlabel('AT8+ %')
	ax.set_ylabel('Age-scaled burden')
	ax.set_xlim([0,30])

	plt.savefig(out_name+'_'+var+'_combine.pdf', dpi=200,  bbox_inches='tight')
	plt.show()


def quant_indel():
	
	tbp_count = '/home/boj924/AD_Tau_PTA/results/indel_spectrum_count.tab'
	ID4_count = '/home/boj924/AD_Tau_PTA/results/ID83_Musical_AD/Decompose_Solution/Activities/Decompose_Solution_Activities.txt'
	burden_in = '/home/boj924/AD_Tau_PTA/results/tau_AD_ctrl_res.tab'

	metafile_in = '/home/boj924/AD_Tau_PTA/metafiles/all_meta.tab'
	clinic_in   = '/home/boj924/AD_Tau_PTA/metafiles/all_clinic'
	AT8pct_in   = '/home/boj924/AD_Tau_PTA/metafiles/quant.csv'
	out_name    = '/home/boj924/AD_Tau_PTA/results/AT8_indel-AD'

	burden = pd.read_csv(burden_in, sep=',', header=0,index_col=0)
	burden = burden.loc[(burden['mutation_type']=='indel') & (burden.index!='1995P_201001E3'),]
	burden.dropna(inplace=True)	

	ID4_count_df = pd.read_csv(ID4_count, sep='\t', header=0, index_col=0)
	ID4_count_df.index = ID4_count_df.index.to_series().apply(lambda x: x.strip('_indel'))
	ID4_pct = ID4_count_df/ID4_count_df.sum(axis=1).values.reshape((-1,1))
	ID4_burden = pd.merge(ID4_pct, burden[['burden']], left_index=True, right_index=True)
	ID4_burden[ID4_pct.columns] = ID4_burden[ID4_pct.columns] * ID4_burden['burden'].values.reshape((-1,1))
	ID4_burden['ID4_rawcount'] = ID4_count_df['ID4']

	tbp_count_df = pd.read_csv(tbp_count, index_col=0, header=0)

	tbp_r = re.compile('2bp_deletion')
	tbp_col = list(filter(tbp_r.search, tbp_count_df.columns))

	del_r = re.compile('deletion')	
	del_col = list(filter(del_r.search, tbp_count_df.columns))

	tbp_pct = tbp_count_df/tbp_count_df.sum(axis=1).values.reshape((-1,1))
	tbp_burden = pd.merge(tbp_pct, burden[['burden']], left_index=True, right_index=True)

	tbp_burden[tbp_pct.columns]  = tbp_burden[tbp_pct.columns] * tbp_burden['burden'].values.reshape((-1,1))
	tbp_burden['twobp_burden']   = tbp_burden[tbp_col].sum(axis=1)
	tbp_burden['twpbp_rawcount'] = tbp_count_df[tbp_col].sum(axis=1)

	tbp_burden['deletion_burden']   = tbp_burden[del_col].sum(axis=1)
	tbp_burden['deletion_rawcount'] = tbp_count_df[del_col].sum(axis=1)

	burden_df = pd.merge(ID4_burden, tbp_burden, left_index=True, right_index=True,suffixes=('','_y'))
	color_palette = 'Pastel1'

	metafile = pd.read_csv(metafile_in, header=0,dtype=str)
	clinic = pd.read_csv(clinic_in, header=0, sep='\t',dtype=str)
	metafile = pd.merge(metafile,clinic, left_on='donor',right_on='Case_ID')
	metafile['Age'] = metafile['Age_(yrs)'].astype(float)

	AT8pct = pd.read_csv(AT8pct_in,sep=',',header=0, index_col=0)
	AT8pct.index = [str(x) for x in AT8pct.index]

	metainfor = pd.merge(metafile, AT8pct, left_on='Case_ID', right_index=True, how='left')
	df = pd.merge(burden_df, metainfor, left_index=True, right_on='sample')
	var = 'burden'
	df_tau = df.loc[df['group'].isin(['AD'])]
	df_tau[AT8pct.columns] = df_tau[AT8pct.columns].apply(norm_to_zero_one)
	print(df_tau['group'])

	res = []
	for col in AT8pct.columns:
		print(col)
		df_tau_sub = df_tau[['Age','group',var, col,"Case_ID"]].dropna()
		md = smf.mixedlm(" %s ~ Age + %s  "%(var,col), df_tau_sub, groups=df_tau_sub["Case_ID"])

		mdf = md.fit()
		rescale_results = mdf.summary().tables[1]
		rescale_results['Coef.'] = rescale_results['Coef.'].astype(float)
		res.append(rescale_results.loc[col,['Coef.','[0.025','0.975]','P>|z|']].tolist())
	res_df = pd.DataFrame(res, columns=['coef','25','975','pval'], index=AT8pct.columns)
	res_df = res_df.astype(float)
	print(res_df)

	cols = ['sup_frontal_np','sup_frontal_dp','sup_frontal_neuron_loss','sup_frontal_nft_count','IHC_neuron_pos_percent']
	res_df = res_df.loc[cols]

	plt.figure(figsize=(7,2), dpi=300)
	ci = [res_df.iloc[::-1]['coef'] - res_df.iloc[::-1]['25'].values, res_df.iloc[::-1]['975'].values - res_df.iloc[::-1]['coef']]
	plt.errorbar(x=res_df.iloc[::-1]['coef'], y=res_df.iloc[::-1].index.values, xerr=ci,
	            color='black',  capsize=3, linestyle='None', linewidth=1,
	            marker="o", markersize=5, mfc="black", mec="black")
	plt.axvline(x=1, linewidth=0.8, linestyle='--', color='black')
	plt.tick_params(axis='both', which='major', labelsize=8)
	plt.xlabel('coef and 95% Confidence Interval', fontsize=8)
	plt.tight_layout()
	plt.savefig('%s_%s_forest_plot.pdf'%(out_name,var), dpi=200,  bbox_inches='tight')
	plt.show()

if __name__ == "__main__":
#	mda_snv()
#	scale_indel()
	quant_indel()
#	quant_snv()




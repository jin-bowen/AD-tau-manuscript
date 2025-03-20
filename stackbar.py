import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib as mpl
import matplotlib.font_manager
plt.rcParams['font.size'] = 10
plt.rcParams['axes.linewidth'] = 1
#plt.style.use('fivethirtyeight')
import numpy as np
from scipy import stats
from scipy.interpolate import UnivariateSpline,CubicSpline
from matplotlib.lines import Line2D
import pandas as pd
import sys
import re

sbs_color_code = {"SBS1": "#ebf5bc",
	"SBS5": "#63d69e",
	"SBS2": "#f8b6b3",
	"SBS13": "#f17fb2",
	"SBS3": "#c4abc4",
	"SBS4": "#bcf2f5",
	"SBS7a": "#b5d7f5",
	"SBS7b": "#9ecef7",
	"SBS7c": "#84bdf0",
	"SBS7d": "#6cb2f0",
	"SBS8": "#dfc4f5",
	"SBS9": "#acf2d0",
	"SBS10a": "#f2aeae",
	"SBS10b": "#f08080",
	"SBS17a": "#d9f7b0",
	"SBS17b": "#8cc63f",
	"SBS40": "#c4c4f5",
	"SBS6": "#faf1dc",
	"SBS14": "#faecca",
	"SBS15": "#fcebc2",
	"SBS20": "#fae4af",
	"SBS21": "#fae1a5",
	"SBS26": "#fcde97",
	"SBS44": "#fad682",
	"SBS89": "darkslateblue",
	"SBS88": "olivedrab",
	"SBS87": "cyan",
	"SBS32": "hotpink",
	"SBS30": "rebeccapurple",
	"SBS25": "lime",
	"SBS19": "forestgreen",
	"SBS27": "#C8C8C8",
	"SBS43": "#C0C0C0",
	"SBS45": "#BEBEBE",
	"SBS46": "#B8B8B8",
	"SBS47": "#B0B0B0",
	"SBS48": "#A9A9A9",
	"SBS49": "#A8A8A8",
	"SBS50": "#A0A0A0",
	"SBS51": "#989898",
	"SBS52": "#909090",
	"SBS53": "#888888",
	"SBS54": "#808080",
	"SBS55": "#787878",
	"SBS56": "#707070",
	"SBS57": "#696969",
	"SBS58": "#686868",
	"SBS59": "#606060",
	"SBS60": "#585858",
	"SBS16": "tab:green",
	"SBS40a":"tab:purple",
	"SBS40b":"tab:red",
	"SBS12": "tab:pink",
	"SBS18": "tab:orange",
	"SBS92": "tab:olive",
	"SBS97": "tab:brown",
	"SBS96": "tab:cyan",
	"SBS39": "deeppink",
	"SBS37": "orangered",
	"SBS42": "blueviolet",
	"SBS22a": "chocolate",
	"SBS33": "darkgreen",}



colors = ["tab:green",
	"tab:purple",
	"tab:pink",
	"tab:orange",
	"tab:olive",
	"tab:brown",
	"tab:red",
	"tab:cyan",
	"deeppink",
	"orangered",
	"blueviolet",
	"chocolate",
	"darkgreen",
	"dodgerblue",
	"mediumvioletred",
	"salmon",
	"magenta",
	"sandybrown",
	"forestgreen",
	"royalblue",
	"orchid",
	"indigo",
	"darkseagreen",
	"blue",
	"palevioletred",
	"darkslateblue",
	"olivedrab",
	"cyan",
	"hotpink",
	"rebeccapurple",
	"lime"]

colors = {"ID1":"tab:green",
	"ID2":"tab:purple",
	"ID3":"tab:pink",
	"ID4":"tab:orange",
	"ID5":"tab:olive",
	"ID6":"tab:brown",
	"ID7":"tab:red",
	"ID8":"tab:cyan",
	"ID9":"deeppink",
	"ID10":"orangered",
	"ID11":"tab:grey",
	"ID12":"chocolate",
	"ID13":"darkgreen",
	"ID14":"dodgerblue",
	"ID15":"mediumvioletred",
	"ID16":"salmon",
	"ID17":"magenta",
	"ID19":"tab:blue",
	"ID83A":"tab:red",
	"ID11b":"tab:grey",
	"ID22":"lime",
	"ID83B":"tab:green",
}


def pct():
	infile='/home/boj924/AD_Tau_PTA/results/ID83/Assignment_Solution/Activities/Assignment_Solution_Activities.txt'
	infile='/home/boj924/AD_Tau_PTA/results/ID83/Decompose_Solution/Activities/Decompose_Solution_Activities.txt'

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

#	infile='/home/boj924/AD_Tau_PTA/results/ID83_Musical_AD/Decompose_Solution/Activities/Decompose_Solution_Activities.txt'
	infile='/home/boj924/AD_Tau_PTA/results/ID83/Decompose_Solution/Activities/Decompose_Solution_Activities.txt'
	infile='/home/boj924/AD_Tau_PTA/results/ID83/De_Novo_Solution/Activities/De_Novo_Activities.txt'

	burden='/home/boj924/AD_Tau_PTA/results/tau_AD_ctrl_res.tab'
	burden_df = pd.read_csv(burden, header=0, sep=',')
	burden_df = burden_df.loc[burden_df['mutation_type']=='indel']
	burden_df = burden_df.loc[burden_df['Cell_ID']!='1995P_201001E3']

	metafile = '/home/boj924/AD_Tau_PTA/metafiles/all_meta.tab'
	meta_df = pd.read_csv(metafile, header=0)

	clinicfile = '/home/boj924/AD_Tau_PTA/metafiles/all_clinic'
	clinic = pd.read_csv(clinicfile, header=0, sep='\t')

	meta_df = pd.merge(meta_df, clinic, left_on='donor', right_on='Case_ID')

	df = pd.read_csv(infile, sep='\t', header=0)
	df['Samples'] = df['Samples'].str.strip('_indel')
	df.set_index('Samples', inplace=True)
	df = df/ df.sum(axis=1).values.reshape((-1, 1))
	df['Samples'] = df.index
	df_reform = df.melt(id_vars='Samples', var_name='IDsig', value_name='pct')
	df_reform = pd.merge(df_reform, meta_df, left_on='Samples', right_on='sample')
	df_reform = pd.merge(df_reform, burden_df[['Cell_ID','burden']], left_on='sample', 
				right_on='Cell_ID')	
	df_reform['IDsig_burden'] = df_reform.apply(lambda row:  row['pct'] * row['burden'], axis=1)

	df_reform['group'] = pd.Categorical(df_reform['group'], ['ctrl','AD','Tau','noTau'])
	df_reform['donor'] = pd.Categorical(df_reform['donor'], ['1353','1647','2208','4556','5222','2207','1456','1995','1828'])
	df_reform['class'] = df_reform.apply(lambda row: 'ctrl'  if row['group'] == 'ctrl'  else 'AD', axis=1)
	df_reform['class'] = pd.Categorical(df_reform['class'], ['ctrl','AD'])

	df_reform = df_reform.sort_values(by=['class', 'Age','group'])

	fig, ax=plt.subplots(constrained_layout=True, figsize=(16,4) )
	sns.histplot(data=df_reform, x="Samples", hue="IDsig", weights='IDsig_burden', \
			palette=colors, multiple="stack", ax=ax, lw=0)
	ax.set_xlabel('')
	ax.set_ylabel('Genome-wide burden')
	ax.set_ylim([0, 8000])
	plt.xticks(rotation=90)
	legend = ax.get_legend()
	legend.set_bbox_to_anchor((1, 1))
	legend.set_frame_on(False)
#	plt.savefig('/home/boj924/AD_Tau_PTA/results/all_Decompose_count_rescale_cosmic.pdf', dpi=300)
	plt.savefig('/home/boj924/AD_Tau_PTA/results/denovo_decompose.pdf', dpi=300)

	plt.show()


if __name__ == "__main__":
        main()
#	snv()



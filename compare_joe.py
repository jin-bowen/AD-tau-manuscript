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

def figure1(df, df2, out_name):

	index=set(df.index).intersection(df2.index)
	df  = df.loc[index]
	df2 = df2.loc[index]
	slope, intercept, r_value, p_value, std_err = stats.linregress(df["rate.per.gb"],\
				df2["rate.per.gb"])

	fig,ax = plt.subplots(constrained_layout=True,figsize=(3,3), frameon = False)

	maxval = max( df["rate.per.gb"])
	maxval = max( maxval, max(df2["rate.per.gb"]))
	maxval = maxval*1.3

	xpos=np.arange(0, maxval)
	ypos=np.arange(0, maxval)

	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)
	ax.plot(xpos, ypos, c='gray', lw=2)
	ax.scatter(x=df["rate.per.gb"], y=df2["rate.per.gb"], c='tab:green', \
		edgecolors='black', linewidth=2)
	ax.set_title('mutation rate(per gbp)\ncorr=%s'%round(r_value,3), loc='right')
	ax.set_xlabel('gatk3 gvcf joint call')
	ax.set_ylabel("Joe's call")

	plt.savefig(out_name+'.png', dpi=200,  bbox_inches='tight' )
	plt.show()

def figure2(df, df2, out_name):

	df2.dropna(inplace=True)
	fig,ax = plt.subplots(constrained_layout=True,figsize=(8,5), frameon = False)

	res1 = stats.linregress(df['Age_(yrs)'], df["burden"])
	res2 = stats.linregress(df2['Age_(yrs)'], df2["burden"])
	xpos=np.arange(0, 110)

	ax.scatter(x=df2['Age_(yrs)'], y=df2["burden"], c='tab:blue', edgecolors='black', \
			label="Joe's call", linewidth=2)

	ax.scatter(x=df['Age_(yrs)'], y=df["burden"], c='tab:red', edgecolors='black', \
			label='gatk3 gvcf joint call', linewidth=2)

	xpos1=np.arange(0, 110)
#	ax.plot(xpos, res2.slope * xpos + res2.intercept, c='tab:blue', lw=1)
	xpos2=np.arange(50, 110)

	ax.set_xlabel('Age')
	ax.set_ylabel('Burden')

	maxval = max( df["burden"])
	maxval = max( maxval, max(df2["burden"]))
	maxval = maxval*1.3
	ax.set_ylim([0,maxval])
#	ax.set_ylim([0,1000])

	plt.legend(loc='center left',bbox_to_anchor=(1.01, 0.5), frameon=False)
	plt.savefig(out_name+'.png', dpi=200,  bbox_inches='tight')
	plt.show()

def main():

	infile   = sys.argv[1]
	infile2  = sys.argv[2]
	metafile = sys.argv[3]
	out_name = sys.argv[4]

	for muttype in ['indel','snv']:
#	for muttype in ['indel']:
		color_palette = 'Pastel1'
		df = pd.read_csv(infile, sep='\t', header=0)
	
		df = df.loc[df['mutation.type']==muttype]
		df['burden'] = df['rate.per.gb'] * 5.845001134
		df['Case_ID'] = df['Case_ID'].astype(str)
	
		meta_file = pd.read_csv(metafile, sep='\t', header=0)
		meta_file['Case_ID'] = meta_file['Case_ID'].astype(str)
		df_merge = pd.merge(meta_file,df,on='Case_ID')
		df_merge.set_index('Cell_ID',inplace=True)
	
		df2 = pd.read_csv(infile2, header=0, index_col=0)
		df2 = df2.loc[df2['muttype']==muttype]
		df2['burden'] = df2['rate.per.gb'] * 5.845001134
		df2['Case_ID'] = df2.index.to_series().apply(lambda x: re.findall(r'\d{3,4}', x)[0]).astype(str)
		df2['Cell_ID'] = df2.index.to_series()
		
		df_merge2 = pd.merge(meta_file, df2, on='Case_ID')
		df_merge2.set_index('Cell_ID',inplace=True)
	
		figure1(df_merge, df_merge2, out_name+'_%s'%muttype+'_RateCorr')
		figure2(df_merge, df_merge2, out_name+'_%s'%muttype+'_BurdenCompare')

if __name__ == "__main__":
	main()



import pandas as pd
import gseapy as gp
import matplotlib.pyplot as plt
import matplotlib as mpl
plt.rcParams['font.serif'] = ['Arial']
plt.rcParams['font.size'] = 10
plt.rcParams['axes.linewidth'] = 1
from matplotlib.lines import Line2D
from scipy import sparse, stats
import numpy as np
import scipy as sp
import re
import sys

def main1():
	
	# number of cell sampled
	pathway_in    = sys.argv[1]
	out           = sys.argv[2]

	pathway_raw = pd.read_csv(pathway_in, header=0, index_col=0)
	pathway = pathway_raw.T

#	human = gp.get_library_name()
#	r = re.compile('^GO.')
#	human_go = list(filter(r.match, human))

	human_go = ['MSigDB_Hallmark_2020']
	for i, ipathway in pathway.iterrows():
		gene_cluster = ipathway[ipathway!=0].index.tolist()
		enr = gp.enrichr(gene_list=gene_cluster, # or gene_list=glist
			organism='human',
			gene_sets=human_go, # kegg is a dict object
			background=None, #"hsapiens_gene_ensembl",
			outdir=None,
			verbose=True)
		print(i,len(gene_cluster))
		print(enr.results.head(5)[['Gene_set','Term','Adjusted P-value']])

def main():
	
	# number of cell sampled
	gene_set_in = sys.argv[1]
	out         = sys.argv[2]

	gene_set = np.loadtxt(gene_set_in, dtype='str')

	human = gp.get_library_name()
	r = re.compile('^GO')
	human_go = list(filter(r.match, human))
	#human_go = ['GO_Biological_Process_2023',
	#	'GO_Cellular_Component_2023',
	#	'GO_Molecular_Function_2023', ]
		#'MSigDB_Hallmark_2020',\
		#'Allen_Brain_Atlas_10x_scRNA_2021']
	enr = gp.enrichr(gene_list=list(gene_set), # or gene_list=glist
		organism='human',
		gene_sets=human_go, # kegg is a dict object
		background=None, #"hsapiens_gene_ensembl",
		outdir=None,
		verbose=True)
	print(enr.results.head(10)[['Gene_set','Term','Adjusted P-value','Genes']])

if __name__ == "__main__":
	main()



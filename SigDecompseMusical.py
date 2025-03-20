import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt
from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc as matGen
from SigProfilerMatrixGenerator.scripts import MutationMatrixGenerator as mutGen
import seaborn as sns
import matplotlib as mpl
import pandas as pd
import time
import scipy as sp
import pickle
import musical
import sys

def utility():

	output_dir="/home/boj924/AD_Tau_PTA/results"
	signatures = output_dir+"/ID83/All_Solutions/ID83_2_Signatures/Signatures/ID83_S2_Signatures.txt"
	activities=output_dir+"/ID83/All_Solutions/ID83_2_Signatures/Activities/ID83_S2_NMF_Activities.txt"
	samples=output_dir+"/ID83/Samples.txt"
	output=output_dir+""
	
	decomp.spa_analyze(signatures=signatures,samples=samples, output=output, \
		signature_database='/home/boj924/AD_Tau_PTA/custome_sig/COSMIC_v3.4_ID_GRCh37.txt', \
		collapse_to_SBS96 = False, genome_build="GRCh37")
	
	
	signatures = output_dir+"/SBS96/All_Solutions/SBS96_2_Signatures/Signatures/SBS96_S2_Signatures.txt"
	activities=output_dir+"/SBS96/All_Solutions/SBS96_2_Signatures/Activities/SBS96_S2_NMF_Activities.txt"
	samples=output_dir+"/SBS96/Samples.txt"
	output=output_dir+""
	
	decomp.spa_analyze(signatures=signatures,samples=samples, output=output, \
		signature_database='/home/boj924/AD_Tau_PTA/Mutagen53_COSMIC.tsv',
		collapse_to_SBS96 = False, genome_build="GRCh37")
	
	decomp.spa_analyze(signatures=signatures,samples=samples, output=output, \
			collapse_to_SBS96 = True, genome_build="GRCh37")

def main():


	matrices = matGen.SigProfilerMatrixGeneratorFunc("test", \
			"GRCh37", "/home/boj924/PrimTau/results/indel_with_AD/",\
			plot=False, exome=False, bed_file=None, \
			chrom_based=False, tsb_stat=False, \
			seqInfo=False, cushion=100)

	reform_mtx = pd.DataFrame(list(matrices.values())[8])

	model = musical.DenovoSig(reform_mtx, 
			min_n_components=1, # Minimum number of signatures to test
			max_n_components=10, # Maximum number of signatures to test
			init='random', # Initialization method
			method='mvnmf', # mvnmf or nmf
			n_replicates=20, # Number of mvnmf/nmf replicates to run per n_components
			ncpu=10, # Number of CPUs to use
			max_iter=100000, # Maximum number of iterations for each mvnmf/nmf run
			bootstrap=True, # Whether or not to bootstrap X for each run
			tol=1e-8, # Tolerance for claiming convergence of mvnmf/nmf
			verbose=1, # Verbosity of output
			normalize_X=False) 
	model.fit()
	with open('/home/boj924/PrimTau/results/denovo_indel_with_AD.pkl', 'wb') as f:
		pickle.dump(model, f, pickle.HIGHEST_PROTOCOL)


#	with open('/home/boj924/PrimTau/results/denovo_indel_with_AD.pkl', 'rb') as f:
#		model = pickle.load(f)
	model.plot_selection()
	plt.show()

	fig = musical.sigplot_bar(model.W, sig_type='Indel83')
	plt.savefig('/home/boj924/PrimTau/results/denovo_sig.pdf')

	catalog = musical.load_catalog('COSMIC_v3p1_Indel')
#	catalog = musical.load_catalog('MuSiCal_v4_Indel_WGS')

	W_catalog = catalog.W

	thresh_grid = np.array([
		0.0001, 0.0002, 0.0005,
		0.001, 0.002, 0.005,
		0.01, 0.02, 0.05,
		0.1, 0.2, 0.5,
		1., 2., 5.])

	model.assign_grid(W_catalog, 
		method_assign='likelihood_bidirectional', # Method for performing matching and refitting
		thresh_match_grid=thresh_grid, # Grid of threshold for matchinng
		thresh_refit_grid=thresh_grid, # Grid of threshold for refitting
		thresh_new_sig=0.0, # De novo signatures with reconstructed cosine similarity below this threshold will be considered novel
		connected_sigs=False, # Whether or not to force connected signatures to co-occur
		clean_W_s=False)

	model.validate_grid(validate_n_replicates=1, # Number of simulation replicates to perform for each grid point
		grid_selection_method='distance')# Method for selecting the best grid point

	with open('/home/boj924/PrimTau/results/denovo_indel_with_AD_COSMIC.pkl', 'wb') as f:
		pickle.dump(model, f, pickle.HIGHEST_PROTOCOL)

	H_s = model.H_s
	print(W_s.columns.tolist())
	
if __name__ == "__main__":
	main()






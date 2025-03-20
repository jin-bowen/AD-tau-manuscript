from SigProfilerExtractor import sigpro as sig
from SigProfilerAssignment import Analyzer as Analyze
from SigProfilerAssignment import decomposition as decomp
import pandas as pd
import sys

def utility():

	output_dir="/home/boj924/PrimTau/results"
	signatures = output_dir+"/ID83_ADPrimaryTau/All_Solutions/ID83_2_Signatures/Signatures/ID83_S2_Signatures.txt"
	activities=output_dir+"/ID83_ADPrimaryTau/All_Solutions/ID83_2_Signatures/Activities/ID83_S2_NMF_Activities.txt"
	samples=output_dir+"/ID83_ADPrimaryTau/Samples.txt"
	output="/home/boj924/PrimTau/"
	
	decomp.spa_analyze(signatures=signatures,samples=samples, output=output, \
		signature_database='/home/boj924/AD_Tau_PTA/custome_sig/MuSiCal_IDsig.csv', \
		collapse_to_SBS96 = False, genome_build="GRCh37")
	
#	decomp.spa_analyze(signatures=signatures,samples=samples, output=output, \
#		signature_database='/home/boj924/AD_Tau_PTA/custome_sig/COSMIC_v3.4_ID_GRCh37.txt', \
#		collapse_to_SBS96 = False, genome_build="GRCh37")
#		
#	signatures = output_dir+"/SBS96/All_Solutions/SBS96_2_Signatures/Signatures/SBS96_S2_Signatures.txt"
#	activities=output_dir+"/SBS96/All_Solutions/SBS96_2_Signatures/Activities/SBS96_S2_NMF_Activities.txt"
#	samples=output_dir+"/SBS96/Samples.txt"
#	output=output_dir+""
#	
#	decomp.spa_analyze(signatures=signatures,samples=samples, output=output, \
#		signature_database='/home/boj924/AD_Tau_PTA/Mutagen53_COSMIC.tsv',
#		collapse_to_SBS96 = False, genome_build="GRCh37")
#	
#	decomp.spa_analyze(signatures=signatures,samples=samples, output=output, \
#			collapse_to_SBS96 = True, genome_build="GRCh37")


def main():
	mut_type = sys.argv[1]
	
	if mut_type == 'snv':
		sig.sigProfilerExtractor("vcf", \
			"/home/boj924/AD_Tau_PTA/results", \
			"/home/boj924/AD_Tau_PTA/results/snv", \
			context_type='SBS96', \
			reference_genome="GRCh37", \
			minimum_signatures=1, maximum_signatures=10, \
			nmf_replicates=200,\
			cpu=-1)
	else:	
		sig.sigProfilerExtractor("vcf", \
			"/home/boj924/AD_Tau_PTA/results", \
			"/home/boj924/AD_Tau_PTA/results/indel", \
			context_type='ID83', \
			reference_genome="GRCh37", \
			minimum_signatures=1, maximum_signatures=10, \
			nmf_replicates=200,\
			cpu=-1)

	
if __name__ == "__main__":
	utility()
#	main()



args=commandArgs(trailingOnly=TRUE)

library(MutationalPatterns)
library(ggplot2)
library(reshape2)
library(BSgenome)

dir = '/home/boj924/AD_Tau_PTA/results/all_snv'

vcf_files <- list.files(c(dir),pattern = "*_snv.vcf", full.names=TRUE)
outname='tauADctrl'

ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
library(ref_genome, character.only = TRUE)

getSampleName = function(x){ gsub('_snv.vcf','', basename(x))  }
sample_names = lapply(vcf_files,getSampleName )
sample_names = unlist(sample_names)

grl <- read_vcfs_as_granges(vcf_files, sample_names, ref_genome)
snv_grl <- get_mut_type(grl, type = "snv")
muts <- mutations_from_vcf(grl[[1]])
context <- mut_context(grl[[1]], ref_genome)
type_context <- type_context(grl[[1]], ref_genome)
type_occurrences <- mut_type_occurrences(grl, ref_genome)

mut_mat <- mut_matrix(vcf_list = grl, ref_genome = ref_genome)
mut_mat <- mut_mat + 0.0001
saveRDS(mut_mat, file = "tau_mut.rds")

## combined sig
#sig_file='/n/data1/bch/genetics/lee/lan/CTE/mut_sig/aging.mutation_pattern.3sigs.tsv'
#aging_signatures=as.matrix(read.delim(sig_file))
#chem_signatures <- read.table("/home/boj924/AD_Tau_PTA/custome_sig/Mutagen53_sub_signature.tsv",sep="\t",header=T,check.names=FALSE)
#sigA = aging_signatures[,2] /sum(aging_signatures[,2])
#all_sig = as.matrix(cbind(chem_signatures[2:54], 'SigA'=sigA))


# transcription bias
# BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")
# BiocManager::install("AnnotationDbi")
# BiocManager::install("GenomicFeatures")
library("TxDb.Hsapiens.UCSC.hg19.knownGene")
genes_hg19 <- GenomicFeatures::genes(TxDb.Hsapiens.UCSC.hg19.knownGene)
strand <- mut_strand(grl[[1]], genes_hg19)
mut_mat_s  <- mut_matrix_stranded(grl, ref_genome, genes_hg19)

meta=read.csv('/home/boj924/AD_Tau_PTA/metafiles/all_meta.tab')
getGrpName = function(x){ meta[meta$sample==x,'group']  }
grpname = lapply(  colnames(mut_mat_s), getGrpName )
grpname =unlist(grpname)

strand_counts <- strand_occurrences(mut_mat_s, by = grpname)
strand_bias <- strand_bias_test(strand_counts)
strand_bias_strict <- strand_bias_test(strand_counts,fdr_cutoffs=0.05) 
fig <- plot_strand(strand_counts, mode = "relative")
ggsave("/home/boj924/AD_Tau_PTA/results/transcription_snv.pdf", plot = fig, width = 8, height = 4, dpi = 300)

fig <- plot_strand_bias(strand_bias_strict)
ggsave("/home/boj924/AD_Tau_PTA/results/transcription_bias_snv.pdf", plot = fig, width = 8, height = 4, dpi = 300)


cols=c("group","newvar","relative_contribution")
strand_counts$newvar = paste(strand_counts$type, strand_counts$strand, sep="_")
strand_counts_reshape = reshape(as.data.frame(strand_counts[,cols]), timevar="group", idvar="newvar", v.names="relative_contribution", direction = "wide")
rownames(strand_counts_reshape) = strand_counts_reshape$newvar
chisq.test(strand_counts_reshape[c(2,3,4,5)])







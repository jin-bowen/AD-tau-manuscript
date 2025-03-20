args=commandArgs(trailingOnly=TRUE)

library(MutationalPatterns)
library(ggplot2)
library(reshape2)
library(BSgenome)
library(dplyr)

source('/home/boj924/AD_Tau_PTA/analysis_scripts/specturm.R')
source('/home/boj924/Tn5-duplex-calling/qc_script/plot_spectrum_mod.R')

#dir = '/home/boj924/Tn5-duplex-calling/results/thres_test'
#outname='thres_test'

dir = '/home/boj924/AD_Tau_PTA/results/all_snv'
outname='all_PTA'

dir = '/home/boj924/AD_Tau_PTA/results/all_snv'
outname='all_PTA'

####compare different subset
vcf_files <- list.files(c(dir),pattern = "*_snv.vcf", full.names=TRUE)
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

library(stringr)
meta_df = data.frame('index'=rownames(type_occurrences))
#meta_df[c('sample', 'thres')] <- str_split_fixed(meta_df$index, "[.]", 2)
#meta_df[c('af', 'thres')] <- str_split_fixed(meta_df$thres, "_", 2)

meta <- read.csv('/home/boj924/AD_Tau_PTA/metafiles/all_meta.tab', header=TRUE)
meta_all = merge(meta_df, meta, by.x="index", by.y="sample")
meta_sub = meta_all

type_occurrences_sub = type_occurrences[meta_sub$index,]
getGrpName = function(x){ meta_sub[meta_sub$index==x,'group']  }
grpname = lapply(rownames(type_occurrences_sub), getGrpName )
grpname = unlist(grpname)

p <- plot_spectrum_mod(type_occurrences_sub, by=grpname, legend = TRUE)
ggsave(paste0("/home/boj924/AD_Tau_PTA/figures/",outname,"_ind.pdf"), plot = p, width = 8, height = 4, dpi = 300)

mut_mat <- mut_matrix(vcf_list = grl, ref_genome = ref_genome)
mut_mat_sub = mut_mat
grpname = lapply(colnames(mut_mat_sub), getGrpName )
colnames(mut_mat_sub) = grpname

ctrl_res = rowSums( mut_mat_sub[,colnames(mut_mat_sub) %in% c('ctrl')] )
ctrl_res = ctrl_res/sum(ctrl_res)
AD_res = rowSums( mut_mat_sub[,colnames(mut_mat_sub) %in% c('AD')] )
AD_res = AD_res/sum(AD_res)
Tau_res = rowSums( mut_mat_sub[,colnames(mut_mat_sub) %in% c('Tau')] )
Tau_res = Tau_res/sum(Tau_res)
noTau_res = rowSums( mut_mat_sub[,colnames(mut_mat_sub) %in% c('noTau')] )
noTau_res = noTau_res/sum(noTau_res)

combine_res = t(rbind(ctrl_res,AD_res, noTau_res, Tau_res))
fig = plot_snv_profile(combine_res)
ggsave(paste0("/home/boj924/AD_Tau_PTA/figures/",outname,"_spectrum.pdf"), plot = fig, width = 8, height = 4, dpi = 300)

rdsfilename='/home/boj924/Tn5-duplex-calling/results/pta_ctrl_spectrum.rds'
saveRDS(combine_res, rdsfilename)


###############################meta-cs
dir = '/home/boj924/Tn5-duplex-calling/results/ctrl_curve_NV_vcf'
outname='ctrl_metacs'

####compare different subset
vcf_files <- list.files(c(dir),pattern = "*_NV.vcf", full.names=TRUE)
ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
library(ref_genome, character.only = TRUE)
getSampleName = function(x){ gsub('_NV.vcf','', basename(x))  }
sample_names = lapply(vcf_files,getSampleName )
sample_names = unlist(sample_names)

grl <- read_vcfs_as_granges(vcf_files, sample_names, ref_genome)
snv_grl <- get_mut_type(grl, type = "snv")
muts <- mutations_from_vcf(grl[[1]])
context <- mut_context(grl[[1]], ref_genome)
type_context <- type_context(grl[[1]], ref_genome)
type_occurrences <- mut_type_occurrences(grl, ref_genome)

library(stringr)
meta_df = data.frame('index'=rownames(type_occurrences))
#meta_df[c('sample', 'thres')] <- str_split_fixed(meta_df$index, "[.]", 2)
#meta_df[c('af', 'thres')] <- str_split_fixed(meta_df$thres, "_", 2)

mut_mat <- mut_matrix(vcf_list = grl, ref_genome = ref_genome)
ctrl_res = rowSums( mut_mat[,colnames(mut_mat)] )
ctrl_res = ctrl_res/sum(ctrl_res)

fig = plot_snv_profile(ctrl_res)
ggsave(paste0("/home/boj924/AD_Tau_PTA/figures/",outname,"_spectrum.pdf"), plot = fig, width = 8, height = 4, dpi = 300)

rdsfilename='/home/boj924/Tn5-duplex-calling/results/ctrl_curve_spectrum.rds'
saveRDS(ctrl_res, rdsfilename)


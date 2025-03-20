args=commandArgs(trailingOnly=TRUE)

library(MutationalPatterns)
library(ggplot2)
library(BSgenome)
library(reshape2)
library(dplyr)

dir =   '/home/boj924/AD_Tau_PTA/results/all_indel'
outname='tauADCtrl'
vcf_files <- list.files(c(dir),pattern = "*_indel.vcf", full.names=TRUE)

ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
library(ref_genome, character.only = TRUE)

getSampleName = function(x){ gsub('_indel.vcf','', basename(x))  }
sample_names = lapply(vcf_files,getSampleName )
sample_names = unlist(sample_names)
grl <- read_vcfs_as_granges(vcf_files, sample_names, ref_genome, type = "all")
indel_grl <- get_mut_type(grl, type = "indel")
indel_grl <- get_indel_context(indel_grl, ref_genome)
indel_counts <- count_indel_contexts(indel_grl)
twobp_row = ifelse(grepl("^2bp", rownames(indel_counts)), TRUE, FALSE)

twobp_count = colSums(indel_counts[twobp_row,])
all_indel_count = colSums(indel_counts)
twobp_percent = twobp_count/all_indel_count
write.csv(twobp_percent, '/home/boj924/AD_Tau_PTA/results/2bp_percent.csv')
#indel_pct =  prop.table(indel_counts, margin=2)


meta=read.csv('/home/boj924/AD_Tau_PTA/metafiles/all_meta.tab')
getGrpName = function(x){ meta[meta$sample==x,'group']  }
grpname = lapply( sample_names, getGrpName )
grpname =unlist(grpname)
colnames(indel_counts) = grpname

ctrl_res = rowSums( indel_counts[,colnames(indel_counts) %in% c('ctrl')] )
ctrl_res = ctrl_res/sum(ctrl_res)
AD_res = rowSums( indel_counts[,colnames(indel_counts) %in% c('AD')] )
AD_res = AD_res/sum(AD_res)
Tau_res = rowSums( indel_counts[,colnames(indel_counts) %in% c('Tau')] )
Tau_res = Tau_res/sum(Tau_res)
noTau_res = rowSums( indel_counts[,colnames(indel_counts) %in% c('noTau')] )
noTau_res = noTau_res/sum(noTau_res)

combine_res = t(rbind(ctrl_res,AD_res, noTau_res, Tau_res))
fig = plot_indel_contexts(combine_res, same_y = TRUE, condensed =T)
ggsave("/home/boj924/AD_Tau_PTA/results/ADindel.pdf", plot = fig, width = 15, height = 5, dpi = 300)


sig_file='/home/boj924/AD_Tau_PTA/custome_sig/MuSiCal_IDsig.csv'
signatures=read.delim(sig_file)
meta_file=read.delim('/home/boj924/AD_Tau_PTA/custome_sig/indel83_lookup.tab', sep=',')

signature_meta = left_join(signatures, meta_file, by=c('MutationType'='indeltype1'))
signatures_select = signature_meta[,c('ID4','ID22')]
row.names(signatures_select) = signature_meta$indeltype2
source('/home/boj924/AD_Tau_PTA/analysis_scripts/specturm.R')
p=plot_indel_profile(signatures_select)
ggsave("/home/boj924/AD_Tau_PTA/results/indel_sig.pdf", plot = p, width = 15, height = 5, dpi = 300)



ID_signature <- read.delim("/home/boj924/AD_Tau_PTA/results/ID83_Musical_AD/Decompose_Solution/Signatures/Decompose_Solution_Signatures.txt")
meta_file=read.delim('/home/boj924/AD_Tau_PTA/custome_sig/indel83_lookup.tab', sep=',')
signature_meta = left_join(ID_signature, meta_file, by=c('MutationType'='indeltype1'))
signatures_select = signature_meta[,c('ID4','ID5','ID8','ID9','ID11b','ID19','ID22')]
row.names(signatures_select) = signature_meta$indeltype2
saveRDS(signatures_select,'/home/boj924/AD_Tau_PTA/results/indel_musical_sig.rds')


ID_signature <- read.delim("custome_sig/MuSiCal_IDsig.csv", sep='\t')
meta_file=read.delim('/home/boj924/AD_Tau_PTA/custome_sig/indel83_lookup.tab', sep=',')
signature_meta = left_join(ID_signature, meta_file, by=c('MutationType'='indeltype1'))
row.names(signature_meta) = signature_meta$indeltype2
signature_meta = select(signature_meta, -c(MutationType,indeltype2))
saveRDS(signature_meta,'/home/boj924/AD_Tau_PTA/results/indel_musical_all.rds')


ID_signature <- read.delim("/home/boj924/AD_Tau_PTA/results/ID83/Decompose_Solution/Signatures/Decompose_Solution_Signatures.txt")
meta_file=read.delim('/home/boj924/AD_Tau_PTA/custome_sig/indel83_lookup.tab', sep=',')
signature_meta = left_join(ID_signature, meta_file, by=c('MutationType'='indeltype1'))
signatures_select = signature_meta[,c('ID3','ID4','ID5','ID8','ID11','ID83A')]
row.names(signatures_select) = signature_meta$indeltype2
saveRDS(signatures_select,'/home/boj924/AD_Tau_PTA/results/indel_cosmic_sig.rds')
source('/home/boj924/AD_Tau_PTA/analysis_scripts/specturm.R')
p=plot_indel_profile(signatures_select)
ggsave("/home/boj924/AD_Tau_PTA/results/indel_sig.pdf", plot = p, width = 8, height = 5, dpi = 300)


ID_signature <- read.delim("/home/boj924/AD_Tau_PTA/results/ID83/De_Novo_Solution/Signatures/De_Novo_Signatures.txt")
meta_file=read.delim('/home/boj924/AD_Tau_PTA/custome_sig/indel83_lookup.tab', sep=',')
signature_meta = left_join(ID_signature, meta_file, by=c('MutationType'='indeltype1'))
signatures_select = signature_meta[,c('ID83B','ID83A')]
row.names(signatures_select) = signature_meta$indeltype2
source('/home/boj924/AD_Tau_PTA/analysis_scripts/specturm.R')
p=plot_indel_profile(signatures_select)
ggsave("/home/boj924/AD_Tau_PTA/results/indel_denovo_sig.pdf", plot = p, width = 8, height = 5, dpi = 300)









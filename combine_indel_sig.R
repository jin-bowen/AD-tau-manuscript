library(MutationalPatterns)
library(ggplot2)
library(reshape2)
library(BSgenome)
library(dplyr)
library(tidyr)
library(forcats)

source('/home/boj924/AD_Tau_PTA/analysis_scripts/specturm.R')
source('/home/boj924/Tn5-duplex-calling/qc_script/plot_spectrum_mod.R')

#### Fit SCAN2 calls to META-CS ds/ssSNV signatures ####
decompose = read.delim("/home/boj924/AD_Tau_PTA/results/ID83_Musical_AD/Decompose_Solution/Activities/Decompose_Solution_Activities.txt")
rownames(decompose) = decompose$Samples
decompose = subset(decompose, select = -Samples)
decompose_percent =  decompose / rowSums(decompose)

assignID_signature <- read.delim("/home/boj924/AD_Tau_PTA/results/ID83_Musical_AD/Decompose_Solution/Signatures/Decompose_Solution_Signatures.txt")
denovoID_signature <- read.delim("/home/boj924/AD_Tau_PTA/results/ID83_Musical_AD/De_Novo_Solution/Signatures/De_Novo_Signatures.txt")
id4_prop = mean(decompose_percent$ID4)
id22_prop = mean(decompose_percent$ID22)

combined_sig = as.data.frame(assignID_signature$ID4 * id4_prop + assignID_signature$ID22 * id22_prop)
rownames(combined_sig) = assignID_signature$MutationType
colnames(combined_sig) = 'value'
combined_sig = combined_sig / colSums(combined_sig)


meta_file=read.delim('/home/boj924/AD_Tau_PTA/custome_sig/indel83_lookup.tab', sep=',')
combined_sig$MutationType = factor(rownames(combined_sig), level=rownames(combined_sig))
signature_meta = left_join(combined_sig, meta_file, by=c('MutationType'='indeltype1'))
signatures_select = as.data.frame(signature_meta[,c('value')])
row.names(signatures_select) = signature_meta$indeltype2
p=plot_indel_profile(signatures_select)


p=plot_indel_profile(combined_sig)
cosmic_signature <- read.delim("/home/boj924/AD_Tau_PTA/results/ID83/Decompose_Solution/Signatures/Decompose_Solution_Signatures.txt")

library(lsa)
cosine(combined_sig$value, cosmic_signature$ID83A)
cosine(cosmic_signature$ID4, cosmic_signature$ID83A)
cosine(assignID_signature$ID4, cosmic_signature$ID83A)







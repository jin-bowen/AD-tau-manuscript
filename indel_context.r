args=commandArgs(trailingOnly=TRUE)

library(MutationalPatterns)
library(ggplot2)
library(reshape2)
library(dplyr)
library(vcfR)
library(dplyr)
library(Biostrings)
library(data.table)

dir =   '/home/boj924/AD_Tau_PTA/results/all_indel'
outname='tauADCtrl'
vcf_files <- list.files(c(dir),pattern = "*_indel.vcf", full.names=TRUE)
getSampleName = function(x){ gsub('_indel.vcf','', basename(x))  }

ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
library(ref_genome, character.only = TRUE)
library(Biostrings)
require(BSgenome.Hsapiens.1000genomes.hs37d5)
genome <- BSgenome.Hsapiens.1000genomes.hs37d5

indel_df = data.frame()
for (file in vcf_files) {
	vcf <- read.vcfR(file)
	sample=getSampleName(file)
	temp = as.data.frame(getFIX(vcf))
	if (ncol(temp) != 7){
		temp = as.data.frame(t(temp))
	}
	temp$sample = sample
	indel_df = rbind(indel_df, temp)
}


#n=10
#indel_anno_df = as.data.table(indel_df)
#indel_anno_df$POS = as.numeric(indel_anno_df$POS)
#indel_2bp_df = indel_anno_df[nchar(REF) - nchar(ALT) == 2,]
#indel_2bp_df[nchar(REF) - nchar(ALT) == 2, deca := as.character(getSeq(BSgenome.Hsapiens.1000genomes.hs37d5, GRanges(seqnames=CHROM, ranges=IRanges(start=POS-3, end=POS+6))))]
#indel_2bp_df$dn = substr(indel_2bp_df$REF,2,3)
#
#will.revcomp = c('GG', 'AA', 'TG', 'AG', 'GA', 'AC')
#indel_2bp_df[,rev_deca := ifelse(dn %in% will.revcomp, reverseComplement(DNAStringSet(deca)), DNAStringSet(deca))]
#
#metafile = '~/AD_Tau_PTA/metafiles/all_meta.tab'
#metadata <- read.table(metafile,sep=",", header = T,  colClasses = c("character", "character", "character"))
#
#indel_2bp_df_meta = inner_join(indel_2bp_df, metadata, by = c("sample" = "sample"))
#
#
#for( grp in c('ctrl','AD','Tau','noTau')){
#indel_2bp_df_sub = indel_2bp_df_meta[indel_2bp_df_meta$group==grp,]
#temp = strsplit(indel_2bp_df_sub$rev_deca, split='')
#temp_mat = matrix(unlist(temp), ncol=n, byrow=TRUE)
#
#occurance_map = data.frame()
#for(i in 1:n) {
#	temp_vec = table(temp_mat[,i])
#	occurance_map = rbind(occurance_map, temp_vec)
#}
#colnames(occurance_map) = c('a','c','g','t')
#proportion <- function(x){
#   rs <- sum(x);
#   return(x / rs);
#}
#
#library(seqLogo)
##create position weight matrix
#mef2 <- apply(occurance_map, 1, proportion)
#mef2 <- makePWM(mef2)
#pdf(file=paste0("/home/boj924/AD_Tau_PTA/figures/motif_",grp, ".pdf"))
#seqLogo(mef2)
#dev.off()
#}


muts = as.data.table(indel_df)
muts$POS = as.numeric(muts$POS)
muts[, dn := substr(REF, 2, 3)]

# note this only makes a 10-mer context around 2bp deletions
muts[nchar(REF) - nchar(ALT) == 2, deca := as.character(getSeq(BSgenome.Hsapiens.1000genomes.hs37d5, GRanges(seqnames=CHROM, ranges=IRanges(start=POS-3, end=POS+6))))]
muts[nchar(REF) - nchar(ALT) == 2, deca.rc := ifelse(dn %in% will.revcomp, reverseComplement(DNAStringSet(deca)), DNAStringSet(deca))]

muts.tmp <- muts[nchar(REF)-nchar(ALT)==2]
metafile = '~/AD_Tau_PTA/metafiles/all_meta.tab'
metadata <- read.table(metafile,sep=",", header = T,  colClasses = c("character", "character", "character"))
muts.tmp = inner_join(muts.tmp, metadata, by = c("sample" = "sample"))

require(ggplot2)
require(ggseqlogo)
ggseqlogo(split(muts.tmp[mut!= '2:Del:R:0',, deca.rc], muts.tmp[, group]), ncol=4) +
    annotate('segment', x = 4.5, xend=6.5, y=1.3, yend=1.3, size=2) +
    annotate('text', x=5.5, y=1.45, label='Deletion') +
    xlab("Local deletion position")+
    ylim(0,2)
ggsave(dev=pdf, file='/home/boj924/AD_Tau_PTA/figures/motif_10mer_deletion_context.pdf', 
		height = 4 , width = 16)



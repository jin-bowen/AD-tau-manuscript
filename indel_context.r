args=commandArgs(trailingOnly=TRUE)

library(MutationalPatterns)
library(ggplot2)
library(BSgenome)
library(reshape2)
library(dplyr)
library(vcfR)

source('/home/boj924/AD_Tau_PTA/analysis_scripts/IndelClass.R')
getSequenceContext <- function(genome, position, chr, offsetL= 10, offsetR=50){
	  sequence <- Biostrings::getSeq(genome, 
				GRanges(seqnames=chr,ranges=IRanges(start=position-offsetL, end=position+offsetR)))
	  context_string <- as.character(sequence)
	  sequenceContext <- list(sequence=sequence,context_string=context_string)
}

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

indel_anno_df = data.frame()
for(i in 1:nrow(indel_df)) {
	row <- indel_df[i,]
	pos = as.numeric(row$POS)
	seq_obj = getSequenceContext(genome = genome,position=pos, chr=row$CHROM,
                        offsetL=6, offsetR=6)
	seq = seq_obj$context_string
	type = attribution_of_indels(genome = genome,in_CHROM=row$CHROM, in_POS=pos, 
			in_REF=row$REF, in_ALT=row$ALT) 
	row$seq = seq
	row$type = type
	indel_anno_df = rbind(indel_anno_df, row)
}

n=13
write.csv(indel_anno_df, '/home/boj924/AD_Tau_PTA/results/indel_type.csv', row.names = FALSE)
library(dplyr)
indel_anno_df = read.csv('/home/boj924/AD_Tau_PTA/results/indel_type.csv')
twobp_row = ifelse(grepl("^2:Del", indel_anno_df$type), TRUE, FALSE)
indel_2bp_df = indel_anno_df[twobp_row,]

metafile = '~/AD_Tau_PTA/metafiles/all_meta.tab'
metadata <- read.table(metafile,sep=",", header = T,  colClasses = c("character", "character", "character"))

indel_2bp_df_meta = inner_join(indel_2bp_df, metadata, by = c("sample" = "sample"))

for( grp in c('ctrl','AD','Tau','noTau')){
indel_2bp_df_sub = indel_2bp_df_meta[indel_2bp_df_meta$group==grp,]
temp = strsplit(indel_2bp_df_sub$seq, split='')
temp_mat = matrix(unlist(temp), ncol=n, byrow=TRUE)

occurance_map = data.frame()
for(i in 1:n) {
	temp_vec = table(temp_mat[,i])
	occurance_map = rbind(occurance_map, temp_vec)
}
colnames(occurance_map) = c('a','c','g','t')
proportion <- function(x){
   rs <- sum(x);
   return(x / rs);
}

library(seqLogo)
#create position weight matrix
mef2 <- apply(occurance_map, 1, proportion)
mef2 <- makePWM(mef2)
pdf(file=paste0("/home/boj924/AD_Tau_PTA/figures/motif_",grp, ".pdf"))
seqLogo(mef2)
dev.off()
}





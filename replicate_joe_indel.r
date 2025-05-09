args=commandArgs(trailingOnly=TRUE)

library(MutationalPatterns)
library(ggplot2)
library(reshape2)
library(dplyr)
library(vcfR)
library(dplyr)
library(Biostrings)
library(data.table)

inhouse_indel_df = read.table('/home/boj924/AD_Tau_PTA/results/indel_type.csv', sep=",", header = T)
metafile = '~/AD_Tau_PTA/metafiles/all_meta.tab'
metadata <- read.table(metafile,sep=",", header = T,  colClasses = c("character", "character", "character"))
indel_2bp_df_meta = inner_join(inhouse_indel_df, metadata, by = c("sample" = "sample"))

will.revcomp = c('GG', 'AA', 'TG', 'AG', 'GA', 'AC')

require(ggplot2)
require(ggseqlogo)
require(BSgenome.Hsapiens.1000genomes.hs37d5)

temp_df = indel_2bp_df_meta
temp_df = as.data.table(temp_df)
temp_df$POS = as.numeric(temp_df$POS)
temp_df[, dn := substr(REF, 2, 3)]

temp_df[nchar(REF) - nchar(ALT) == 2, deca := as.character(getSeq(BSgenome.Hsapiens.1000genomes.hs37d5, GRanges(seqnames=CHROM, ranges=IRanges(start=POS-3, end=POS+6))))]
temp_df[nchar(REF) - nchar(ALT) == 2, deca.rc := ifelse(dn %in% will.revcomp, reverseComplement(DNAStringSet(deca)), DNAStringSet(deca))]

temp_df.tmp <- temp_df[nchar(REF)-nchar(ALT)==2]

ggseqlogo(temp_df.tmp[(temp_df.tmp$type!='2:Del:R:0') & (temp_df.tmp$group!='ctrl'), deca.rc]) +
    annotate('segment', x = 4.5, xend=6.5, y=1.3, yend=1.3, size=2) +
    annotate('text', x=5.5, y=1.45, label='Deletion') +
    xlab("Local deletion position")+
    ylim(0,2)
ggsave(dev=pdf, file='/home/boj924/AD_Tau_PTA/figures/ID-AD_motif_AD.pdf', 
		height = 4 , width = 4)

ggseqlogo(temp_df.tmp[(temp_df.tmp$type=='2:Del:R:0') & (temp_df.tmp$group!='ctrl'), deca.rc]) +
    annotate('segment', x = 4.5, xend=6.5, y=1.3, yend=1.3, size=2) +
    annotate('text', x=5.5, y=1.45, label='Deletion') +
    xlab("Local deletion position")+
    ylim(0,2)
ggsave(dev=pdf, file='/home/boj924/AD_Tau_PTA/figures/2bp0_motif_AD.pdf', 
		height = 4 , width = 4)


ggseqlogo(split(temp_df.tmp[temp_df.tmp$type!= '2:Del:R:0',deca.rc], temp_df.tmp[temp_df.tmp$type!= '2:Del:R:0', group]), ncol=4) +
    annotate('segment', x = 4.5, xend=6.5, y=1.3, yend=1.3, size=2) +
    annotate('text', x=5.5, y=1.45, label='Deletion') +
    xlab("Local deletion position")+
    ylim(0,2)
ggsave(dev=pdf, file='/home/boj924/AD_Tau_PTA/figures/ID-AD_motif_split.pdf', 
		height = 4 , width = 16)

ggseqlogo(split(temp_df.tmp[temp_df.tmp$type== '2:Del:R:0',deca.rc], temp_df.tmp[temp_df.tmp$type== '2:Del:R:0', group]), ncol=4) +
    annotate('segment', x = 4.5, xend=6.5, y=1.3, yend=1.3, size=2) +
    annotate('text', x=5.5, y=1.45, label='Deletion') +
    xlab("Local deletion position")+
    ylim(0,2)
ggsave(dev=pdf, file='/home/boj924/AD_Tau_PTA/figures/2bp0_motif_split.pdf', 
		height = 4 , width = 16)

library(data.table)
signature = read.table('/home/boj924/AD_Tau_PTA/results/ID83_Musical_AD/Assignment_Solution/Signatures/Assignment_Solution_Signatures.txt', sep="\t", header = T)
row.names(signature) = signature$MutationType

twobp_row = ifelse(grepl("2:Del", rownames(signature)), TRUE, FALSE)
twobp_list = rownames(signature)[twobp_row]
signature_sub = signature[twobp_list,c('ID4','ID22')]
signature_sub$ID4 = signature_sub$ID4 / sum(signature_sub$ID4)
signature_sub$ID22 = signature_sub$ID22 / sum(signature_sub$ID22)

pinkpeak = colSums(signature_sub[c(-1),])
other = signature_sub[c('2:Del:R:0'),]

df = rbind(pinkpeak,other)
df$type = c('pinkpeak','2:Del:R:0')
long <- melt(df, id.vars = c('type'))

vysg = ggplot(long, aes(x="", y=value, fill=type, )) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0)
vysg<-vysg+facet_wrap(~ variable)
ggsave(dev=pdf, file='/home/boj924/AD_Tau_PTA/figures/2bp_piechart.pdf',
                height = 4 , width = 8)








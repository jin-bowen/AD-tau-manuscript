library(goseq)
library(ggplot2)
library(reshape2)

set.seed(1000)
permutation_round=100
numDEInCat_threshold=10
numInCat_threshold=1000

my.GOseq <- function(input_hit,input_length,output_file,permutation)
{
	pwf=nullp(input_hit,bias.data=input_length)
	GO=goseq(pwf,"hg19","geneSymbol")
	GO_filtered=GO[GO$numDEInCat>=numDEInCat_threshold&GO$numInCat<=numInCat_threshold,]
	GO_filtered$over_represented_pvalue_corrected=p.adjust(GO_filtered$over_represented_pvalue,method="fdr")
	if(!permutation)
	{
		#print(head(GO_filtered))
		write.table(GO_filtered,file=output_file,quote=F,sep="\t",row.names=F)
	}
	GO
}

my.indexDiff <- function(go,input1,input2)
{
	a=which(input1$category==go)/nrow(input1)
	b=which(input2$category==go)/nrow(input2)
	if(length(a)==0)
	{
		a=1
	}
	if(length(b)==0)
	{
		b=1
	}
	a-b
}

my.fisherTest <- function(a,b,c,d)
{
	if(a/c>=b/d)
	{
		fisher.test(matrix(c(a,b,c-a,d-b),nrow=2),alternative="greater")$p.value
	}
	else
	{
		fisher.test(matrix(c(a,b,c-a,d-b),nrow=2),alternative="less")$p.value
	}
}


cols <- c("Tau" = "#D62728", "noTau" = "#FFC0CB")
workdir='/n/data1/bwh/pathology/miller/lab/AD_Tau/results_hg19/'
metafile = '~/AD_Tau_PTA/metafiles/all_meta.tab'
snv_mutation_file='/n/data1/bwh/pathology/miller/lab/AD_Tau/results_hg19/annovar/snv_list.annovar.variant_function'
indel_mutation_file='/n/data1/bwh/pathology/miller/lab/AD_Tau/results_hg19/annovar/snv_list.annovar.variant_function'

# read metafile
metadata_full <- read.table(metafile,sep=",", header = T,  colClasses = c("character", "character", "character"))
metadata  = metadata_full[metadata_full$group %in% c('Tau','noTau'),]
head(metadata)

library(tidyr)
library(stringr)
# read in snv data
snv_mutation <- read.table(snv_mutation_file)
colnames(snv_mutation) <- c("region", "gene", "chr", "start", "end", "ref", "alt", "donor", "single_cell_ID")
snv_mutation = snv_mutation %>% separate(donor, c('donor', 'batch'))

snv_mutation$gene <- str_remove(snv_mutation$gene, "\\(.*\\)$") # remove mutations corresponding to more than one genes
snv_mutation <- snv_mutation[!str_detect(snv_mutation$gene, ","), ]
snv_mutation_meta = dplyr::left_join(snv_mutation, metadata, by=c("single_cell_ID" = "sample"), suffix = c(".x", ""),)
snv_mutation_tau = snv_mutation_meta[snv_mutation_meta$group=='Tau',]
snv_mutation_notau = snv_mutation_meta[snv_mutation_meta$group=='noTau',]

gene_length=read.delim("/n/data1/bwh/pathology/miller/lab/ref/hg19.refGene.length",header=T, sep='\t')
gene_length_deduped=gene_length[!duplicated(gene_length$gene),]

Tau_exonintron_ready=gene_length_deduped$gene %in% snv_mutation_tau$gene[snv_mutation_tau$region %in% c("exonic","intronic")]
names(Tau_exonintron_ready)=gene_length_deduped$gene
Tau_raw=my.GOseq(Tau_exonintron_ready,gene_length_deduped$merged,"/n/data1/bwh/pathology/miller/lab/AD_Tau/results_hg19/annovar/Tau.exonic_intronic.GO.tsv",FALSE)

noTau_exonintron_ready=gene_length_deduped$gene %in% snv_mutation_notau$gene[snv_mutation_notau$region %in% c("exonic","intronic")]
names(noTau_exonintron_ready)=gene_length_deduped$gene
noTau_raw=my.GOseq(noTau_exonintron_ready,gene_length_deduped$merged,"/n/data1/bwh/pathology/miller/lab/AD_Tau/results_hg19/annovar/noTau.exonic_intronic.GO.tsv",FALSE)

Tau_go_list=as.character(Tau_raw$category[Tau_raw$over_represented_pvalue<0.01&Tau_raw$numDEInCat>=numDEInCat_threshold&Tau_raw$numInCat<=numInCat_threshold])
go_list=unique(c(Tau_go_list))
noTau_go_list=as.character(noTau_raw$category[noTau_raw$over_represented_pvalue<0.01&noTau_raw$numDEInCat>=numDEInCat_threshold&noTau_raw$numInCat<=numInCat_threshold])
go_list=unique(c(noTau_go_list))

go_list=unique(c(Tau_go_list,noTau_go_list))
gene_list=unique(c(names(Tau_exonintron_ready)[Tau_exonintron_ready],names(noTau_exonintron_ready)[noTau_exonintron_ready]))
gene_map=getgo(gene_list,"hg19","geneSymbol")
idx_diff=sapply(go_list,function(x){my.indexDiff(x,Tau_raw[Tau_raw$numDEInCat>=numDEInCat_threshold&Tau_raw$numInCat<=numInCat_threshold,],noTau_raw[noTau_raw$numDEInCat>=numDEInCat_threshold&noTau_raw$numInCat<=numInCat_threshold,])})

idx_diff_permutation=numeric(0)
for(i in 1:permutation_round)
{
	Tau_exonintron_ready=gene_length_deduped$gene %in% snv_mutation_tau$gene[snv_mutation_tau$region %in% c("exonic","intronic")]
	Tau_exonintron_ready=sample(Tau_exonintron_ready)
	names(Tau_exonintron_ready)=gene_length_deduped$gene
	Tau_permutation=my.GOseq(Tau_exonintron_ready,gene_length_deduped$merged,NA,TRUE)
	noTau_exonintron_ready=gene_length_deduped$gene %in% snv_mutation_notau$gene[snv_mutation_notau$region %in% c("exonic","intronic")]
	noTau_exonintron_ready=sample(noTau_exonintron_ready)
	names(noTau_exonintron_ready)=gene_length_deduped$gene
	noTau_permutation=my.GOseq(noTau_exonintron_ready,gene_length_deduped$merged,NA,TRUE)
	
	idx_diff_permutation=rbind(idx_diff_permutation,sapply(go_list,function(x){my.indexDiff(x,Tau_permutation[Tau_permutation$numDEInCat>=numDEInCat_threshold&Tau_permutation$numInCat<=numInCat_threshold,],noTau_permutation[noTau_permutation$numDEInCat>=numDEInCat_threshold&noTau_permutation$numInCat<=numInCat_threshold,])}))
}
write.table(idx_diff_permutation,file="/n/data1/bwh/pathology/miller/lab/AD_Tau/results_hg19/annovar/Tau-noTau.exonic_intronic.raw.tsv",quote=F,sep="\t",row.names=F)

idx_diff_permutation=read.delim("/n/data1/bwh/pathology/miller/lab/AD_Tau/results_hg19/annovar/Tau-noTau.exonic_intronic.raw.tsv",header=T)
go_table=data.frame(go_list,
					sapply(go_list,function(x){sum(snv_mutation_tau$gene[which(snv_mutation_tau$region %in% c("exonic","intronic"))] %in% names(gene_map)[grep(x,gene_map)],na.rm=T)}),
					sapply(go_list,function(x){sum(snv_mutation_notau$gene[which(snv_mutation_notau$region %in% c("exonic","intronic"))] %in% names(gene_map)[grep(x,gene_map)],na.rm=T)}),
					sum(snv_mutation_tau$region %in% c("exonic","intronic"),na.rm=T),
					sum(snv_mutation_notau$region %in% c("exonic","intronic"),na.rm=T),
					sapply(go_list,function(x){Tau_raw$numDEInCat[Tau_raw$category==x]}),
					sapply(go_list,function(x){noTau_raw$numDEInCat[noTau_raw$category==x]}),
					sapply(go_list,function(x){Tau_raw$numInCat[Tau_raw$category==x]}),
					sapply(go_list,function(x){Tau_raw$over_represented_pvalue[Tau_raw$category==x]}),
					sapply(go_list,function(x){noTau_raw$over_represented_pvalue[noTau_raw$category==x]}),
					sapply(go_list,function(x){Tau_raw$term[Tau_raw$category==x]}),
					sapply(go_list,function(x){Tau_raw$ontology[Tau_raw$category==x]}))
colnames(go_table)=c("category","hit_Tau","hit_noTau","total_Tau","total_noTau","numDE_Tau","numDE_noTau","numInCat","p_Tau","p_noTau","term","ontology")
go_table$corrected_p_Tau=p.adjust(go_table$p_Tau,method="fdr")
go_table$corrected_p_noTau=p.adjust(go_table$p_noTau,method="fdr")
go_table$p_permutation=sapply(1:nrow(go_table),function(x){if(idx_diff[x]<0){sum(idx_diff[x]>idx_diff_permutation[,x])/length(idx_diff_permutation[,x])}else{sum(idx_diff[x]<idx_diff_permutation[,x])/length(idx_diff_permutation[,x])}})
go_table$corrected_p_permutation=p.adjust(go_table$p_permutation,method="fdr")
write.table(go_table,file="/n/data1/bwh/pathology/miller/lab/AD_Tau/results_hg19/annovar/Tau_vs_noTau.exonic_intronic.GO.tsv",quote=F,sep="\t",row.names=F)


go_table = read.delim("/n/data1/bwh/pathology/miller/lab/AD_Tau/results_hg19/annovar/Tau_vs_noTau.exonic_intronic.GO.tsv", header=T)
ready=melt(go_table[go_table$corrected_p_Tau<0.01&go_table$corrected_p_noTau<0.01,c(1,11,13,14)])
colnames(ready)=c("GO","Term","Variable","Pvalue")
ready$Term=factor(ready$Term,levels=rev(unique(ready$Term)))
ready$Tissue=factor(sapply(ready$Variable,function(x){if(grepl("noTau",x)){"noTau"}else if(grepl("Tau",x)){"Tau"}}),levels=c("noTau","Tau"))
ready$Pconvert=-log10(ready$Pvalue)

noTau_ready = ready[ready$Tissue=='Tau',]
Tau_ready   = ready[ready$Tissue=='noTau',]
reform = merge(x=Tau_ready, y=noTau_ready, by='GO')
reform$diff = reform$Pconvert.y -  reform$Pconvert.x
reform = reform[order(-reform$diff), ]

reform_top10_term = head(reform[order(-reform$Pconvert.x),],n=10)
reform_top10_term = reform_top10_term$Term.y
ready_top10 = ready[ready$Term %in% reform_top10_term,]

p <- ggplot(ready_top10,aes(x=Term,y=Pconvert,fill=Tissue))
p + geom_col(position="dodge") + geom_hline(yintercept=-log10(0.01),linetype="dashed") + scale_fill_discrete(name="") + xlab("") + ylab("-log10(FDR-adjusted P-value)") + coord_flip() + theme(legend.position="bottom")

pdf("/n/data1/bwh/pathology/miller/lab/AD_Tau/results_hg19/annovar/Tau_vs_noTau.exonic_intronic.abs.pdf",height=5,width=8)

p <- ggplot(ready_top10,aes(x=Term,y=Pconvert,fill=factor(Tissue,levels=c("Tau","noTau"))))
p + geom_col(position="dodge",color="black") + geom_hline(yintercept=-log10(0.01),linetype="dashed") + 
	scale_fill_manual(values=c("#d62728","#D87AB2"),name="") +
	xlab("") + ylab("-log10(FDR-adjusted P-value)") + 
	coord_flip() + theme_classic() + theme(legend.position="bottom",text=element_text(size=16))

dev.off()




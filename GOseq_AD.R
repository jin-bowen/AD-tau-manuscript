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


cols <- c("AD" = "#F57F20", "ctrl" = "#2278B5")
workdir='/n/data1/bwh/pathology/miller/lab/AD_AD/results_hg19/'
metafile = '~/AD_AD_PTA/metafiles/all_meta.tab'
snv_mutation_file='/n/data1/bwh/pathology/miller/lab/AD_AD/results_hg19/annovar/snv_list.annovar.variant_function'
indel_mutation_file='/n/data1/bwh/pathology/miller/lab/AD_AD/results_hg19/annovar/snv_list.annovar.variant_function'

# read metafile
metadata_full <- read.table(metafile,sep=",", header = T,  colClasses = c("character", "character", "character"))
metadata  = metadata_full[metadata_full$group %in% c('AD','ctrl'),]
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
snv_mutation_AD = snv_mutation_meta[snv_mutation_meta$group=='AD',]
snv_mutation_ctrl = snv_mutation_meta[snv_mutation_meta$group=='ctrl',]

gene_length=read.delim("/n/data1/bwh/pathology/miller/lab/ref/hg19.refGene.length",header=T, sep='\t')
gene_length_deduped=gene_length[!duplicated(gene_length$gene),]

AD_exonintron_ready=gene_length_deduped$gene %in% snv_mutation_AD$gene[snv_mutation_AD$region %in% c("exonic","intronic")]
names(AD_exonintron_ready)=gene_length_deduped$gene
AD_raw=my.GOseq(AD_exonintron_ready,gene_length_deduped$merged,"/n/data1/bwh/pathology/miller/lab/AD_AD/results_hg19/annovar/AD.exonic_intronic.GO.tsv",FALSE)

ctrl_exonintron_ready=gene_length_deduped$gene %in% snv_mutation_ctrl$gene[snv_mutation_ctrl$region %in% c("exonic","intronic")]
names(ctrl_exonintron_ready)=gene_length_deduped$gene
ctrl_raw=my.GOseq(ctrl_exonintron_ready,gene_length_deduped$merged,"/n/data1/bwh/pathology/miller/lab/AD_AD/results_hg19/annovar/ctrl.exonic_intronic.GO.tsv",FALSE)

AD_go_list=as.character(AD_raw$category[AD_raw$over_represented_pvalue<0.01&AD_raw$numDEInCat>=numDEInCat_threshold&AD_raw$numInCat<=numInCat_threshold])
go_list=unique(c(AD_go_list))
ctrl_go_list=as.character(ctrl_raw$category[ctrl_raw$over_represented_pvalue<0.01&ctrl_raw$numDEInCat>=numDEInCat_threshold&ctrl_raw$numInCat<=numInCat_threshold])
go_list=unique(c(ctrl_go_list))

go_list=unique(c(AD_go_list,ctrl_go_list))
gene_list=unique(c(names(AD_exonintron_ready)[AD_exonintron_ready],names(ctrl_exonintron_ready)[ctrl_exonintron_ready]))
gene_map=getgo(gene_list,"hg19","geneSymbol")
idx_diff=sapply(go_list,function(x){my.indexDiff(x,AD_raw[AD_raw$numDEInCat>=numDEInCat_threshold&AD_raw$numInCat<=numInCat_threshold,],ctrl_raw[ctrl_raw$numDEInCat>=numDEInCat_threshold&ctrl_raw$numInCat<=numInCat_threshold,])})

idx_diff_permutation=numeric(0)
for(i in 1:permutation_round)
{
	AD_exonintron_ready=gene_length_deduped$gene %in% snv_mutation_AD$gene[snv_mutation_AD$region %in% c("exonic","intronic")]
	AD_exonintron_ready=sample(AD_exonintron_ready)
	names(AD_exonintron_ready)=gene_length_deduped$gene
	AD_permutation=my.GOseq(AD_exonintron_ready,gene_length_deduped$merged,NA,TRUE)
	ctrl_exonintron_ready=gene_length_deduped$gene %in% snv_mutation_ctrl$gene[snv_mutation_ctrl$region %in% c("exonic","intronic")]
	ctrl_exonintron_ready=sample(ctrl_exonintron_ready)
	names(ctrl_exonintron_ready)=gene_length_deduped$gene
	ctrl_permutation=my.GOseq(ctrl_exonintron_ready,gene_length_deduped$merged,NA,TRUE)
	
	idx_diff_permutation=rbind(idx_diff_permutation,sapply(go_list,function(x){my.indexDiff(x,AD_permutation[AD_permutation$numDEInCat>=numDEInCat_threshold&AD_permutation$numInCat<=numInCat_threshold,],ctrl_permutation[ctrl_permutation$numDEInCat>=numDEInCat_threshold&ctrl_permutation$numInCat<=numInCat_threshold,])}))
}
write.table(idx_diff_permutation,file="/n/data1/bwh/pathology/miller/lab/AD_AD/results_hg19/annovar/AD-ctrl.exonic_intronic.raw.tsv",quote=F,sep="\t",row.names=F)

idx_diff_permutation=read.delim("/n/data1/bwh/pathology/miller/lab/AD_AD/results_hg19/annovar/AD-ctrl.exonic_intronic.raw.tsv",header=T)
go_table=data.frame(go_list,
					sapply(go_list,function(x){sum(snv_mutation_AD$gene[which(snv_mutation_AD$region %in% c("exonic","intronic"))] %in% names(gene_map)[grep(x,gene_map)],na.rm=T)}),
					sapply(go_list,function(x){sum(snv_mutation_ctrl$gene[which(snv_mutation_ctrl$region %in% c("exonic","intronic"))] %in% names(gene_map)[grep(x,gene_map)],na.rm=T)}),
					sum(snv_mutation_AD$region %in% c("exonic","intronic"),na.rm=T),
					sum(snv_mutation_ctrl$region %in% c("exonic","intronic"),na.rm=T),
					sapply(go_list,function(x){AD_raw$numDEInCat[AD_raw$category==x]}),
					sapply(go_list,function(x){ctrl_raw$numDEInCat[ctrl_raw$category==x]}),
					sapply(go_list,function(x){AD_raw$numInCat[AD_raw$category==x]}),
					sapply(go_list,function(x){AD_raw$over_represented_pvalue[AD_raw$category==x]}),
					sapply(go_list,function(x){ctrl_raw$over_represented_pvalue[ctrl_raw$category==x]}),
					sapply(go_list,function(x){AD_raw$term[AD_raw$category==x]}),
					sapply(go_list,function(x){AD_raw$ontology[AD_raw$category==x]}))
colnames(go_table)=c("category","hit_AD","hit_ctrl","total_AD","total_ctrl","numDE_AD","numDE_ctrl","numInCat","p_AD","p_ctrl","term","ontology")
go_table$corrected_p_AD=p.adjust(go_table$p_AD,method="fdr")
go_table$corrected_p_ctrl=p.adjust(go_table$p_ctrl,method="fdr")
go_table$p_permutation=sapply(1:nrow(go_table),function(x){if(idx_diff[x]<0){sum(idx_diff[x]>idx_diff_permutation[,x])/length(idx_diff_permutation[,x])}else{sum(idx_diff[x]<idx_diff_permutation[,x])/length(idx_diff_permutation[,x])}})
go_table$corrected_p_permutation=p.adjust(go_table$p_permutation,method="fdr")
write.table(go_table,file="/n/data1/bwh/pathology/miller/lab/AD_AD/results_hg19/annovar/AD_vs_ctrl.exonic_intronic.GO.tsv",quote=F,sep="\t",row.names=F)


go_table = read.delim("/n/data1/bwh/pathology/miller/lab/AD_AD/results_hg19/annovar/AD_vs_ctrl.exonic_intronic.GO.tsv", header=T)
ready=melt(go_table[go_table$corrected_p_AD<0.01&go_table$corrected_p_ctrl<0.01,c(1,11,13,14)])
colnames(ready)=c("GO","Term","Variable","Pvalue")
ready$Term=factor(ready$Term,levels=rev(unique(ready$Term)))
ready$Tissue=factor(sapply(ready$Variable,function(x){if(grepl("ctrl",x)){"ctrl"}else if(grepl("AD",x)){"AD"}}),levels=c("ctrl","AD"))
ready$Pconvert=-log10(ready$Pvalue)

ctrl_ready = ready[ready$Tissue=='AD',]
AD_ready   = ready[ready$Tissue=='ctrl',]
reform = merge(x=AD_ready, y=ctrl_ready, by='GO')
reform$diff = reform$Pconvert.y -  reform$Pconvert.x
reform = reform[order(-reform$diff), ]

reform_top10_term = head(reform[order(-reform$Pconvert.x),],n=10)
reform_top10_term = reform_top10_term$Term.y
ready_top10 = ready[ready$Term %in% reform_top10_term,]

p <- ggplot(ready_top10,aes(x=Term,y=Pconvert,fill=Tissue))
p + geom_col(position="dodge") + geom_hline(yintercept=-log10(0.01),linetype="dashed") + scale_fill_discrete(name="") + xlab("") + ylab("-log10(FDR-adjusted P-value)") + coord_flip() + theme(legend.position="bottom")

pdf("/n/data1/bwh/pathology/miller/lab/AD_AD/results_hg19/annovar/AD_vs_ctrl.exonic_intronic.abs.pdf",height=5,width=8)

p <- ggplot(ready_top10,aes(x=Term,y=Pconvert,fill=factor(Tissue,levels=c("AD","ctrl"))))
p + geom_col(position="dodge",color="black") + geom_hline(yintercept=-log10(0.01),linetype="dashed") + 
	scale_fill_manual(values=c("#d62728","#D87AB2"),name="") +
	xlab("") + ylab("-log10(FDR-adjusted P-value)") + 
	coord_flip() + theme_classic() + theme(legend.position="bottom",text=element_text(size=16))

dev.off()




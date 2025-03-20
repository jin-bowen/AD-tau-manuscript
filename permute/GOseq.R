##Adapted from August's script

options(stringsAsFactors = FALSE)

# BiocManager::install("goseq")
# BiocManager::install("org.Hs.eg.db")
httr::set_config(httr::config(ssl_verifypeer = FALSE))
library(goseq)
library(ggplot2)
library(reshape2)
library(dplyr)

set.seed(1000)
permutation_round=1000
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

gene_length=read.delim("~/lan/ref/hg19_refGene.length.tsv",header=F)
colnames(gene_length)=c("Gene_symbol","Transcript","Gene_length","Exon_length")
gene_length_deduped=gene_length[!duplicated(gene_length$Gene_symbol),]

workdir='/n/data1/bwh/pathology/miller/lab/AD_Tau_PTA/results_hg19/'
metafile = '~/AD_Tau_PTA/metafiles/all_meta.tab'
snv_mutation_file='/n/data1/bwh/pathology/miller/lab/AD_Tau_PTA/results_hg19/annovar/snv_list.annovar.variant_function'
indel_mutation_file='/n/data1/bwh/pathology/miller/lab/AD_Tau_PTA/results_hg19/annovar/snv_list.annovar.variant_function'



#### SNV ####
# read in metadata
metadata <- read.table('/home/boj924/AD_Tau_PTA/metafiles/all_meta.tab', sep=',', header = T)
num_Tau =length(unique(metadata[metadata$group=='Tau',]$sample))
num_noTau =length(unique(metadata[metadata$group=='noTau',]$sample))

# read in snv data
snv_mutation <- read.table(snv_mutation_file)
colnames(snv_mutation) <- c("region", "gene", "chr", "start", "end", "ref", "alt", "donor", "single_cell_ID")
snv_mutation <- inner_join(snv_mutation, metadata, by = c("single_cell_ID" = "sample"))

tau_exon_ready= snv_mutation$gene[snv_mutation$region %in% c("exonic") & (snv_mutation$group == "Tau")]
notau_exon_ready= snv_mutation$gene[snv_mutation$region %in% c("exonic") & (snv_mutation$group == "noTau")]

CTE_go_list=as.character(CTE_raw$category[CTE_raw$over_represented_pvalue<0.01&CTE_raw$numDEInCat>=numDEInCat_threshold&CTE_raw$numInCat<=numInCat_threshold])
normal_go_list=as.character(normal_raw$category[normal_raw$over_represented_pvalue<0.01&normal_raw$numDEInCat>=numDEInCat_threshold&normal_raw$numInCat<=numInCat_threshold])
# saveRDS(CTE_go_list, "snv/CTE_go_list.exonic_intronic.rds")
# saveRDS(normal_go_list, "snv/normal_go_list.exonic_intronic.rds")
saveRDS(CTE_go_list, "snv/CTE_go_list.exonic.rds")
saveRDS(normal_go_list, "snv/normal_go_list.exonic.rds")
go_list=unique(c(CTE_go_list,normal_go_list))
# gene_list=unique(c(names(CTE_exonintron_ready)[CTE_exonintron_ready],names(normal_exonintron_ready)[normal_exonintron_ready]))
gene_list=unique(c(names(CTE_exon_ready)[CTE_exon_ready],names(normal_exon_ready)[normal_exon_ready]))
gene_map=getgo(gene_list,"hg19","geneSymbol")
idx_diff=sapply(go_list,function(x){
  my.indexDiff(x, 
               CTE_raw[CTE_raw$numDEInCat>=numDEInCat_threshold&CTE_raw$numInCat<=numInCat_threshold,],
               normal_raw[normal_raw$numDEInCat>=numDEInCat_threshold&normal_raw$numInCat<=numInCat_threshold,])
  })

idx_diff_permutation=numeric(0)
for(i in 1:permutation_round)
{
	# # exon + intron
  # CTE_exonintron_ready=gene_length_deduped$Gene_symbol %in% merged$gene[merged$region %in% c("exonic","intronic") & merged$Clinical=="CTE"]
  # CTE_exonintron_ready=sample(CTE_exonintron_ready)
  # names(CTE_exonintron_ready)=gene_length_deduped$Gene_symbol
  # CTE_permutation=my.GOseq(CTE_exonintron_ready,gene_length_deduped$Gene_length,NA,TRUE)
  # normal_exonintron_ready=gene_length_deduped$Gene_symbol %in% merged$gene[merged$region %in% c("exonic","intronic")& merged$Clinical=="none"]
  # normal_exonintron_ready=sample(normal_exonintron_ready)
  # names(normal_exonintron_ready)=gene_length_deduped$Gene_symbol
  # normal_permutation=my.GOseq(normal_exonintron_ready,gene_length_deduped$Gene_length,NA,TRUE)

  # exon only
  CTE_exon_ready=gene_length_deduped$Gene_symbol %in% merged$gene[merged$region %in% c("exonic") & merged$Clinical=="CTE"]
  CTE_exon_ready=sample(CTE_exon_ready)
  names(CTE_exon_ready)=gene_length_deduped$Gene_symbol
  CTE_permutation=my.GOseq(CTE_exon_ready,gene_length_deduped$Gene_length,NA,TRUE)
  normal_exon_ready=gene_length_deduped$Gene_symbol %in% merged$gene[merged$region %in% c("exonic")& merged$Clinical=="none"]
  normal_exon_ready=sample(normal_exon_ready)
  names(normal_exon_ready)=gene_length_deduped$Gene_symbol
  normal_permutation=my.GOseq(normal_exon_ready,gene_length_deduped$Gene_length,NA,TRUE)

	  
  idx_diff_permutation=rbind(idx_diff_permutation,sapply(go_list,function(x){my.indexDiff(x,CTE_permutation[CTE_permutation$numDEInCat>=numDEInCat_threshold&CTE_permutation$numInCat<=numInCat_threshold,],normal_permutation[normal_permutation$numDEInCat>=numDEInCat_threshold&normal_permutation$numInCat<=numInCat_threshold,])}))
}
# write.table(idx_diff_permutation,file="snv/CTE_vs_normal.exonic_intronic.raw.tsv",quote=F,sep="\t",row.names=F)
write.table(idx_diff_permutation,file="snv/CTE_vs_normal.exonic.raw.tsv",quote=F,sep="\t",row.names=F)

# idx_diff_permutation=read.delim("snv/CTE_vs_normal.exonic_intronic.raw.tsv",header=T)
idx_diff_permutation=read.delim("snv/CTE_vs_normal.exonic.raw.tsv",header=T)
go_table=data.frame(go_list,
                    # # exon + intron
                    # sapply(go_list,function(x){sum(merged$gene[which(merged$region %in% c("exonic","intronic") & merged$Clinical=="CTE")] %in% names(gene_map)[grep(x,gene_map)],na.rm=T)}),
                    # sapply(go_list,function(x){sum(merged$gene[which(merged$region %in% c("exonic","intronic") & merged$Clinical=="none")] %in% names(gene_map)[grep(x,gene_map)],na.rm=T)}),
                    # sum(merged$region %in% c("exonic","intronic") & merged$Clinical=="CTE",na.rm=T),
                    # sum(merged$region %in% c("exonic","intronic") & merged$Clinical=="none",na.rm=T),
                    
                    # exon only
                    sapply(go_list,function(x){sum(merged$gene[which(merged$region %in% c("exonic") & merged$Clinical=="CTE")] %in% names(gene_map)[grep(x,gene_map)],na.rm=T)}),
                    sapply(go_list,function(x){sum(merged$gene[which(merged$region %in% c("exonic") & merged$Clinical=="none")] %in% names(gene_map)[grep(x,gene_map)],na.rm=T)}),
                    sum(merged$region %in% c("exonic") & merged$Clinical=="CTE",na.rm=T),
                    sum(merged$region %in% c("exonic") & merged$Clinical=="none",na.rm=T),
                    
                    sapply(go_list,function(x){CTE_raw$numDEInCat[CTE_raw$category==x]}),
                    sapply(go_list,function(x){normal_raw$numDEInCat[normal_raw$category==x]}),
                    sapply(go_list,function(x){CTE_raw$numInCat[CTE_raw$category==x]}),
                    sapply(go_list,function(x){CTE_raw$over_represented_pvalue[CTE_raw$category==x]}),
                    sapply(go_list,function(x){normal_raw$over_represented_pvalue[normal_raw$category==x]}),
                    sapply(go_list,function(x){CTE_raw$term[CTE_raw$category==x]}),
                    sapply(go_list,function(x){CTE_raw$ontology[CTE_raw$category==x]}))
colnames(go_table)=c("category","hit_CTE","hit_normal","total_CTE","total_normal","numDE_CTE","numDE_normal","numInCat","p_CTE","p_normal","term","ontology")
go_table$corrected_p_CTE=p.adjust(go_table$p_CTE,method="fdr")
go_table$corrected_p_normal=p.adjust(go_table$p_normal,method="fdr")
#go_table$p_fisher=sapply(1:nrow(go_table),function(x){my.fisherTest(go_table$hit_CTE[x],go_table$hit_normal[x],go_table$total_CTE[x],go_table$total_normal[x])})
go_table$p_permutation=sapply(1:nrow(go_table),function(x){
  if(idx_diff[x]<0){
    sum(idx_diff[x]>idx_diff_permutation[,x])/length(idx_diff_permutation[,x])
  }else{
    sum(idx_diff[x]<idx_diff_permutation[,x])/length(idx_diff_permutation[,x])
  }})
go_table$corrected_p_permutation=p.adjust(go_table$p_permutation,method="fdr")
# write.table(go_table,file="snv/CTE_vs_normal.exonic_intronic.GO.tsv",quote=F,sep="\t",row.names=F)
write.table(go_table,file="snv/CTE_vs_normal.exonic.GO.tsv",quote=F,sep="\t",row.names=F)



#### indel ####
mutation <- readRDS("/n/data1/bch/genetics/lee/shulin/CTE/results/CTE.scan2.indel.summary.mutation_site.variant_function.rds")
merged 
# > unique(mutation$region)
# [1] "intergenic"          "intronic"            "UTR3"                "downstream"         
# [5] "UTR5"                "ncRNA_intronic"      "upstream"            "exonic"             
# [9] "ncRNA_exonic"        "ncRNA_splicing"      "splicing"            "upstream;downstream"
# [13] "exonic;splicing"    

# CTE_exonintron_ready=gene_length_deduped$Gene_symbol %in% merged$gene[merged$region %in% c("exonic","intronic","exonic;splicing") & merged$Clinical=="CTE"]
# names(CTE_exonintron_ready)=gene_length_deduped$Gene_symbol
# CTE_raw=my.GOseq(CTE_exonintron_ready,gene_length_deduped$Gene_length,"indel/CTE.exonic_intronic.GO.tsv",FALSE)
CTE_exon_ready=gene_length_deduped$Gene_symbol %in% merged$gene[merged$region %in% c("exonic","exonic;splicing") & merged$Clinical=="CTE"]
names(CTE_exon_ready)=gene_length_deduped$Gene_symbol
CTE_raw=my.GOseq(CTE_exon_ready,gene_length_deduped$Gene_length,"indel/CTE.exonic.GO.tsv",FALSE)

# normal_exonintron_ready=gene_length_deduped$Gene_symbol %in% merged$gene[merged$region %in% c("exonic","intronic","exonic;splicing")& merged$Clinical=="none"]
# names(normal_exonintron_ready)=gene_length_deduped$Gene_symbol
# normal_raw=my.GOseq(normal_exonintron_ready,gene_length_deduped$Gene_length,"indel/normal.exonic_intronic.GO.tsv",FALSE)
normal_exon_ready=gene_length_deduped$Gene_symbol %in% merged$gene[merged$region %in% c("exonic","exonic;splicing")& merged$Clinical=="none"]
names(normal_exon_ready)=gene_length_deduped$Gene_symbol
normal_raw=my.GOseq(normal_exon_ready,gene_length_deduped$Gene_length,"indel/normal.exonic.GO.tsv",FALSE)

CTE_go_list=as.character(CTE_raw$category[CTE_raw$over_represented_pvalue<0.01&CTE_raw$numDEInCat>=numDEInCat_threshold&CTE_raw$numInCat<=numInCat_threshold])
normal_go_list=as.character(normal_raw$category[normal_raw$over_represented_pvalue<0.01&normal_raw$numDEInCat>=numDEInCat_threshold&normal_raw$numInCat<=numInCat_threshold])
# saveRDS(CTE_go_list, "indel/CTE_go_list.exonic_intronic.rds")
# saveRDS(normal_go_list, "indel/normal_go_list.exonic_intronic.rds")
saveRDS(CTE_go_list, "indel/CTE_go_list.exonic.rds")
saveRDS(normal_go_list, "indel/normal_go_list.exonic.rds")
go_list=unique(c(CTE_go_list,normal_go_list))
# gene_list=unique(c(names(CTE_exonintron_ready)[CTE_exonintron_ready],names(normal_exonintron_ready)[normal_exonintron_ready]))
gene_list=unique(c(names(CTE_exon_ready)[CTE_exon_ready],names(normal_exon_ready)[normal_exon_ready]))
gene_map=getgo(gene_list,"hg19","geneSymbol")
idx_diff=sapply(go_list,function(x){
  my.indexDiff(x, 
               CTE_raw[CTE_raw$numDEInCat>=numDEInCat_threshold&CTE_raw$numInCat<=numInCat_threshold,],
               normal_raw[normal_raw$numDEInCat>=numDEInCat_threshold&normal_raw$numInCat<=numInCat_threshold,])
})

idx_diff_permutation=numeric(0)
for(i in 1:permutation_round)
{
  # # exon + intron
  # CTE_exonintron_ready=gene_length_deduped$Gene_symbol %in% merged$gene[merged$region %in% c("exonic","intronic","exonic;splicing") & merged$Clinical=="CTE"]
  # CTE_exonintron_ready=sample(CTE_exonintron_ready)
  # names(CTE_exonintron_ready)=gene_length_deduped$Gene_symbol
  # CTE_permutation=my.GOseq(CTE_exonintron_ready,gene_length_deduped$Gene_length,NA,TRUE)
  # normal_exonintron_ready=gene_length_deduped$Gene_symbol %in% merged$gene[merged$region %in% c("exonic","intronic","exonic;splicing")& merged$Clinical=="none"]
  # normal_exonintron_ready=sample(normal_exonintron_ready)
  # names(normal_exonintron_ready)=gene_length_deduped$Gene_symbol
  # normal_permutation=my.GOseq(normal_exonintron_ready,gene_length_deduped$Gene_length,NA,TRUE)
  
  # exon only
  CTE_exon_ready=gene_length_deduped$Gene_symbol %in% merged$gene[merged$region %in% c("exonic","exonic;splicing") & merged$Clinical=="CTE"]
  CTE_exon_ready=sample(CTE_exon_ready)
  names(CTE_exon_ready)=gene_length_deduped$Gene_symbol
  CTE_permutation=my.GOseq(CTE_exon_ready,gene_length_deduped$Gene_length,NA,TRUE)
  normal_exon_ready=gene_length_deduped$Gene_symbol %in% merged$gene[merged$region %in% c("exonic","exonic;splicing")& merged$Clinical=="none"]
  normal_exon_ready=sample(normal_exon_ready)
  names(normal_exon_ready)=gene_length_deduped$Gene_symbol
  normal_permutation=my.GOseq(normal_exon_ready,gene_length_deduped$Gene_length,NA,TRUE)
  
  idx_diff_permutation=rbind(idx_diff_permutation,sapply(go_list,function(x){my.indexDiff(x,CTE_permutation[CTE_permutation$numDEInCat>=numDEInCat_threshold&CTE_permutation$numInCat<=numInCat_threshold,],normal_permutation[normal_permutation$numDEInCat>=numDEInCat_threshold&normal_permutation$numInCat<=numInCat_threshold,])}))
}
# write.table(idx_diff_permutation,file="indel/CTE_vs_normal.exonic_intronic.raw.tsv",quote=F,sep="\t",row.names=F)
write.table(idx_diff_permutation,file="indel/CTE_vs_normal.exonic.raw.tsv",quote=F,sep="\t",row.names=F)

# idx_diff_permutation=read.delim("indel/CTE_vs_normal.exonic_intronic.raw.tsv",header=T)
idx_diff_permutation=read.delim("indel/CTE_vs_normal.exonic.raw.tsv",header=T)
go_table=data.frame(go_list,
                    # # exon + intron
                    # sapply(go_list,function(x){sum(merged$gene[which(merged$region %in% c("exonic","intronic","exonic;splicing") & merged$Clinical=="CTE")] %in% names(gene_map)[grep(x,gene_map)],na.rm=T)}),
                    # sapply(go_list,function(x){sum(merged$gene[which(merged$region %in% c("exonic","intronic","exonic;splicing") & merged$Clinical=="none")] %in% names(gene_map)[grep(x,gene_map)],na.rm=T)}),
                    # sum(merged$region %in% c("exonic","intronic","exonic;splicing") & merged$Clinical=="CTE",na.rm=T),
                    # sum(merged$region %in% c("exonic","intronic","exonic;splicing") & merged$Clinical=="none",na.rm=T),

                    # exon only
                    sapply(go_list,function(x){sum(merged$gene[which(merged$region %in% c("exonic","exonic;splicing") & merged$Clinical=="CTE")] %in% names(gene_map)[grep(x,gene_map)],na.rm=T)}),
                    sapply(go_list,function(x){sum(merged$gene[which(merged$region %in% c("exonic","exonic;splicing") & merged$Clinical=="none")] %in% names(gene_map)[grep(x,gene_map)],na.rm=T)}),
                    sum(merged$region %in% c("exonic","exonic;splicing") & merged$Clinical=="CTE",na.rm=T),
                    sum(merged$region %in% c("exonic","exonic;splicing") & merged$Clinical=="none",na.rm=T),

                    sapply(go_list,function(x){CTE_raw$numDEInCat[CTE_raw$category==x]}),
                    sapply(go_list,function(x){normal_raw$numDEInCat[normal_raw$category==x]}),
                    sapply(go_list,function(x){CTE_raw$numInCat[CTE_raw$category==x]}),
                    sapply(go_list,function(x){CTE_raw$over_represented_pvalue[CTE_raw$category==x]}),
                    sapply(go_list,function(x){normal_raw$over_represented_pvalue[normal_raw$category==x]}),
                    sapply(go_list,function(x){CTE_raw$term[CTE_raw$category==x]}),
                    sapply(go_list,function(x){CTE_raw$ontology[CTE_raw$category==x]}))
colnames(go_table)=c("category","hit_CTE","hit_normal","total_CTE","total_normal","numDE_CTE","numDE_normal","numInCat","p_CTE","p_normal","term","ontology")
go_table$corrected_p_CTE=p.adjust(go_table$p_CTE,method="fdr")
go_table$corrected_p_normal=p.adjust(go_table$p_normal,method="fdr")
#go_table$p_fisher=sapply(1:nrow(go_table),function(x){my.fisherTest(go_table$hit_CTE[x],go_table$hit_normal[x],go_table$total_CTE[x],go_table$total_normal[x])})
go_table$p_permutation=sapply(1:nrow(go_table),function(x){
  if(idx_diff[x]<0){
    sum(idx_diff[x]>idx_diff_permutation[,x])/length(idx_diff_permutation[,x])
  }else{
    sum(idx_diff[x]<idx_diff_permutation[,x])/length(idx_diff_permutation[,x])
  }})
go_table$corrected_p_permutation=p.adjust(go_table$p_permutation,method="fdr")
# write.table(go_table,file="indel/CTE_vs_normal.exonic_intronic.GO.tsv",quote=F,sep="\t",row.names=F)
write.table(go_table,file="indel/CTE_vs_normal.exonic.GO.tsv",quote=F,sep="\t",row.names=F)

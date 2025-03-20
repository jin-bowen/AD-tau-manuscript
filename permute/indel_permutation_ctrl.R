library(ggplot2)
library(stringr)
library(dplyr)
library(reshape2)
library(tidyverse)
library(tidyr)

calculate_pvalue <- function(x ,y){
  data <- data.frame(x = as.numeric(as.character(x)), y = y)
  lm_result <- summary(lm(y ~ x, data = data))
  pvalue <- lm_result$coefficients["x", "Pr(>|t|)"]
  if(pvalue == 0){
    pvalue <- "<2e-16"
  }else{
    pvalue <- format(pvalue, scientific = T, digits = 4)
  }
  pvalue
}

cols <- c("Tau" = "#D62728", "noTau" = "#FFC0CB")
workdir='/n/data1/bwh/pathology/miller/lab/ctrl_hg19/'
metafile = '~/AD_Tau_PTA/metafiles/all_meta.tab'
indel_mutation_file='/n/data1/bwh/pathology/miller/lab/ctrl_hg19/annovar/indel_list.annovar.variant_function'
indel_permutation_file='/n/data1/bwh/pathology/miller/lab/ctrl_hg19/permute/annovar/perms_by_cell_indel.combined.variant_function'
group_num = 10

# read metafile
metadata <- read.table(metafile,sep=",", header = T,  colClasses = c("character", "character", "character"))

# read in indel data
indel_mutation <- read.table(indel_mutation_file)
indel_mutation <- indel_mutation[,c(1:5, 11,12,15)]
colnames(indel_mutation) <- c("region", "gene", "chr", "start", "end", "ref", "alt", "sample_infor")
indel_mutation = indel_mutation %>% separate(sample_infor, sep=',', c('donor', 'single_cell_ID')) %>% separate(donor, sep='_', c('donor', 'batch'))

genic_indel_mutation <- indel_mutation[indel_mutation$region %in% c("exonic", "exonic;splicing",
                                                  "intronic", "splicing",
                                                  "UTR3", "UTR5", "UTR5;UTR3"),]

genic_indel_mutation$gene <- str_remove(genic_indel_mutation$gene, "\\(.*\\)$") # remove mutations corresponding to more than one genes
genic_indel_mutation <- genic_indel_mutation[!str_detect(genic_indel_mutation$gene, ","), ]

genic_indel_mutation_meta = dplyr::left_join(genic_indel_mutation, metadata, by=c("single_cell_ID" = "sample"), suffix = c(".x", ""),)

indel_mutation_num <- data.frame(table(genic_indel_mutation_meta$gene, genic_indel_mutation_meta$group))
colnames(indel_mutation_num) <- c("gene", "group", "mut_number")

rdsfilename=paste0(workdir,'indel_mutation.rds')
saveRDS(indel_mutation_num, rdsfilename)

permutation_file = indel_permutation_file
permutation <- read.table(permutation_file, header = F)
permutation <- permutation[,c(1:7, 15)]
colnames(permutation) <- c("region", "gene", "chr", "start", "end", "ref", "alt", "perm.id")
permutation = permutation %>% separate(perm.id, sep=';', c('perm.id', 'muttype'))
genic_permutation <- permutation[permutation$region %in% c("exonic", "exonic;splicing",
				"intronic", "splicing", "UTR3", "UTR5", "UTR5;UTR3"),]

genic_permutation$gene <- str_remove(genic_permutation$gene, "\\(.*\\)$") 
genic_permutation <- genic_permutation[!str_detect(genic_permutation$gene, ","), ]
saveRDS(genic_permutation, paste0(workdir,'genic_indel_mutation.rds'))


permutation_num <- data.frame(table(genic_permutation$gene, genic_permutation$perm.id))
colnames(permutation_num) <- c("gene", "perm.id", "permutation_number")

indel_mutation_perm = right_join(indel_mutation_num, permutation_num, by = c("gene"))
indel_mutation_perm$mut_number[is.na(indel_mutation_perm$mut_number)] <- 0
indel_mutation_perm$permutation_number[is.na(indel_mutation_perm$permutation_number)] <- 0
indel_mutation_perm$group=indel_mutation_perm$clinic
saveRDS(indel_mutation_perm, paste0(workdir,'indel_mutation_perm.rds'))

#########################################################################################
group_num = 10
gene_expression_dir='/n/data1/bwh/pathology/miller/lab/AD_Tau/results_hg19/'
rdsfilename=paste0(gene_expression_dir,'cte_expr_level.rds')
expr_level=readRDS(rdsfilename)

expr_level_mutation <- inner_join(expr_level, indel_mutation_perm, by = c("gene" = "gene"))

average_mut_num <- expr_level_mutation %>% group_by(decile, perm.id) %>% summarise("mutation_number" = sum(mut_number),
                                                                                   "permutation_number" = sum(permutation_number))
average_mut_num$enrichment_ratio <- average_mut_num$mutation_number/average_mut_num$permutation_number
saveRDS(average_mut_num, paste0(workdir, 'indel_average_mut_num.rds'))


average_mut_num_summary <- average_mut_num %>% group_by(decile) %>% summarise("average_enrichment" = mean(enrichment_ratio),
                                                             "sd_enrichment" = sd(enrichment_ratio),
                                                             "mutation_num" = sum(mutation_number),
                                                             "permutation_num" = sum(permutation_number))

#ggplot(average_mut_num_summary, aes(x = decile, y = average_enrichment)) +
#  geom_line(size=2) +
#  geom_point(size=2, colour="black") +
#  scale_color_manual(values = cols) +
#  geom_errorbar(aes(ymin = average_enrichment-sd_enrichment, ymax = average_enrichment+sd_enrichment), width = 0.2) +
#  theme_classic() +
#  ylim(c(0, 2)) +
#  labs(x = "Gene expression levels",
#       y = "Mutation enrichment ratio \n (observed/expected)",
#       title = "Total sSNVs")
#ggsave(paste0(workdir,"indel_enrichment_sc.pdf"), height = 4, width = 6)



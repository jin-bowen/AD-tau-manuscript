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
snv_mutation_file='/n/data1/bwh/pathology/miller/lab/ctrl_hg19/annovar/snv_list.annovar.variant_function'
indel_mutation_file='/n/data1/bwh/pathology/miller/lab/ctrl_hg19/annovar/indel_list.annovar.variant_function'
snv_permutation_file='/n/data1/bwh/pathology/miller/lab/ctrl_hg19/permute/annovar/perms_by_cell_snv.combined.variant_function'
indel_permutation_file='/n/data1/bwh/pathology/miller/lab/ctrl_hg19/permute/annovar/perms_by_cell_indel.combined.variant_function'
group_num = 10

# read metafile
metadata <- read.table(metafile,sep=",", header = T,  colClasses = c("character", "character", "character"))

library(tidyr)
# read in snv data
snv_mutation <- read.table(snv_mutation_file)
colnames(snv_mutation) <- c("region", "gene", "chr", "start", "end", "ref", "alt", "donor", "single_cell_ID")

genic_snv_mutation <- snv_mutation[snv_mutation$region %in% c("exonic", "exonic;splicing",
                                                  "intronic", "splicing",
                                                  "UTR3", "UTR5", "UTR5;UTR3"),]

genic_snv_mutation$gene <- str_remove(genic_snv_mutation$gene, "\\(.*\\)$") # remove mutations corresponding to more than one genes
genic_snv_mutation <- genic_snv_mutation[!str_detect(genic_snv_mutation$gene, ","), ]

library(dplyr)
genic_snv_mutation_meta = dplyr::left_join(genic_snv_mutation, metadata, by=c("single_cell_ID" = "sample"), suffix = c(".x", ""),)

snv_mutation_num <- data.frame(table(genic_snv_mutation_meta$gene, genic_snv_mutation_meta$group))
colnames(snv_mutation_num) <- c("gene", "group", "mut_number")

rdsfilename=paste0(workdir,'snv_mutation.rds')
saveRDS(snv_mutation_num, rdsfilename)

permutation_file = snv_permutation_file
permutation <- read.table(permutation_file, header = F)
permutation <- permutation[,c(1:7, 15)]
colnames(permutation) <- c("region", "gene", "chr", "start", "end", "ref", "alt", "perm.id")
genic_permutation <- permutation[permutation$region %in% c("exonic", "exonic;splicing",
				"intronic", "splicing", "UTR3", "UTR5", "UTR5;UTR3"),]

genic_permutation$gene <- str_remove(genic_permutation$gene, "\\(.*\\)$") 
genic_permutation <- genic_permutation[!str_detect(genic_permutation$gene, ","), ]
saveRDS(genic_permutation, paste0(workdir,'genic_snv_mutation.rds'))


permutation_num <- data.frame(table(genic_permutation$gene, genic_permutation$perm.id))
colnames(permutation_num) <- c("gene", "perm.id", "permutation_number")

snv_mutation_perm = right_join(snv_mutation_num, permutation_num, by = c("gene"))
snv_mutation_perm$mut_number[is.na(snv_mutation_perm$mut_number)] <- 0
snv_mutation_perm$permutation_number[is.na(snv_mutation_perm$permutation_number)] <- 0
snv_mutation_perm$group=snv_mutation_perm$clinic
saveRDS(snv_mutation_perm, paste0(workdir,'snv_mutation_perm.rds'))

#########################################################################################
group_num = 10
gene_expression_dir='/n/data1/bwh/pathology/miller/lab/AD_Tau/results_hg19/'
rdsfilename=paste0(gene_expression_dir,'cte_expr_level.rds')
expr_level=readRDS(rdsfilename)

expr_level_mutation <- inner_join(expr_level, snv_mutation_perm, by = c("gene" = "gene"))

average_mut_num <- expr_level_mutation %>% group_by(decile, perm.id) %>% summarise("mutation_number" = sum(mut_number),
                                                                                   "permutation_number" = sum(permutation_number))
average_mut_num$enrichment_ratio <- average_mut_num$mutation_number/average_mut_num$permutation_number

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
#ggsave(paste0(workdir,"snv_enrichment_sc.pdf"), height = 4, width = 6)

################################################################################################
############################################################################################
ref_genome="BSgenome.Hsapiens.UCSC.hg19"
library(ref_genome, character.only = T)
chr_orders <- c(paste0("chr", 1:22), "chrX", "chrY", "chrM")
library(MutationalPatterns)
library(GenomicRanges)

average_mut_num["sigA_mutation"] <- NA
average_mut_num["sigC_mutation"] <- NA

decile_grange_list <- GRangesList()
for (i in 1:group_num){
	gene_name <- unique(expr_level_mutation[expr_level_mutation$decile == i, "gene"])
	genic_mutation_decile <- genic_snv_mutation_meta[(genic_snv_mutation_meta$gene %in% gene_name),]

	genic_mutation_decile <- genic_mutation_decile[order(genic_mutation_decile$start),]
	genic_mutation_decile <- genic_mutation_decile[order(genic_mutation_decile$chr),]

	decile_grange_list[[paste0("decile_", i)]] <- GRanges(seqnames = paste0("chr", genic_mutation_decile$chr),
				ranges = IRanges(start = genic_mutation_decile$start,
						end = genic_mutation_decile$end),
						ref = genic_mutation_decile$ref,
						alt = genic_mutation_decile$alt)
}
chr_length <- seqlengths(Hsapiens)
seqlengths(decile_grange_list) <- chr_length[names(seqlengths(decile_grange_list))]
seqlevels(decile_grange_list) <- seqlevels(decile_grange_list)[order(factor(seqlevels(decile_grange_list), levels = chr_orders))]
genome(decile_grange_list) = 'hg19'
mut_mat <- mut_matrix(decile_grange_list, ref_genome = ref_genome)
aging_signature_all <- readRDS("/n/data1/bch/genetics/lee/shulin/CTE/ref/aging_signature_all.rds")
aging_signature_all <- aging_signature_all[,-2] # remove sig B
aging_signature_all <- aging_signature_all[row.names(mut_mat),]

aging_signature_fit <- fit_to_signatures(mut_mat, aging_signature_all)
aging_signature_contribution <- apply(aging_signature_fit$contribution, 2, function(x){x/sum(x)})
colnames(aging_signature_contribution) <- gsub("decile_", "", colnames(aging_signature_contribution))
average_mut_num$sigA_mutation <- aging_signature_contribution[1, unlist(lapply(colnames(aging_signature_contribution), rep, 1000))]
average_mut_num$sigC_mutation <- aging_signature_contribution[2, unlist(lapply(colnames(aging_signature_contribution), rep, 1000))]

rdsfilename=paste0(workdir,'snv_average_mut_num.rds')
saveRDS(average_mut_num, rdsfilename)

############################################################################
genic_permutation = readRDS(paste0(workdir,'genic_snv_mutation.rds'))
aging_signature_contribution_summary <- c()
for (perm_round in 1:1000){
	decile_grange_list <- GRangesList()
	for (i in 1:group_num){
		gene_name <- unique(expr_level_mutation[expr_level_mutation$decile == i, "gene"])
		genic_permutation_decile <- genic_permutation[(genic_permutation$gene %in% gene_name)  &
								(genic_permutation$perm.id == perm_round), ]
		genic_permutation_decile <- genic_permutation_decile[order(genic_permutation_decile$start),]
		genic_permutation_decile <- genic_permutation_decile[order(genic_permutation_decile$chr),]
		decile_grange_list[[paste0("decile_", i)]] <- GRanges(seqnames = genic_permutation_decile$chr,
					ranges = IRanges(start = genic_permutation_decile$start,
							end = genic_permutation_decile$end),
							ref = genic_permutation_decile$ref,
							alt = genic_permutation_decile$alt)
	}
	chr_length <- seqlengths(Hsapiens)
	seqlengths(decile_grange_list) <- chr_length[names(seqlengths(decile_grange_list))]
	seqlevels(decile_grange_list) <- seqlevels(decile_grange_list)[order(factor(seqlevels(decile_grange_list), levels = chr_orders))]
	genome(decile_grange_list) = 'hg19'
	mut_mat <- mut_matrix(decile_grange_list, ref_genome = ref_genome)
	aging_signature_all <- readRDS("/n/data1/bch/genetics/lee/shulin/CTE/ref/aging_signature_all.rds")
	aging_signature_all <- aging_signature_all[,-2] # remove sig B
	aging_signature_all <- aging_signature_all[row.names(mut_mat),]

	aging_signature_fit <- fit_to_signatures(mut_mat, aging_signature_all)
	aging_signature_contribution <- apply(aging_signature_fit$contribution, 2, function(x){x/sum(x)})

	aging_signature_contribution <- data.frame(t(aging_signature_contribution[1:2,]))
	rownames(aging_signature_contribution) <- gsub("decile_", "", rownames(aging_signature_contribution))
	colnames(aging_signature_contribution) <- c("sigA_permutation", "sigC_permutation")
	aging_signature_contribution[,"perm.id"] <- perm_round	
	aging_signature_contribution[,"decile"] <- rownames(aging_signature_contribution)

	aging_signature_contribution_summary <- rbind(aging_signature_contribution_summary, aging_signature_contribution)
	if(perm_round %% 100 == 0){print(perm_round)}
}
saveRDS(aging_signature_contribution_summary, paste0(workdir,'snv_permutation_sig.rds'))

############################################################################
#average_mut_num <- readRDS(paste0(workdir,'snv_permutation_sig.rds'))
#average_mut_num$decile <- as.factor(average_mut_num$decile)
#average_mut_num$perm.id <- as.integer(as.character(average_mut_num$perm.id))
#
#average_mut_num_rescale <- average_mut_num %>% group_by(decile, group, perm.id) %>% summarise( 
#							mutation_number = sum(mutation_number),
#							sigA_mutation = sum(sigA_mutation),
#							sigC_mutation = sum(sigC_mutation),
#							sigA_permutation = sum(sigA_permutation),
#							sigC_permutation = sum(sigC_permutation),
#							permutation_number = sum(permutation_number))
#
#
#average_mut_num_rescale["sigA_enrichment"] <-
#  (average_mut_num_rescale$sigA_mutation * average_mut_num_rescale$mutation_number) / (average_mut_num_rescale$sigA_permutation * average_mut_num_rescale$permutation_number)
#average_mut_num_rescale[is.infinite(average_mut_num_rescale$sigA_enrichment), "sigA_enrichment"] <- 1
#average_mut_num_rescale["sigC_enrichment"] <-
#  (average_mut_num_rescale$sigC_mutation * average_mut_num_rescale$mutation_number) / (average_mut_num_rescale$sigC_permutation * average_mut_num_rescale$permutation_number)
#average_mut_num_rescale[is.infinite(average_mut_num_rescale$sigC_enrichment), "sigC_enrichment"] <- 1

##############################################################
#average_mut_num_summary_sigA <- average_mut_num_rescale %>% group_by(decile, group) %>% summarise("average_enrichment" = mean(sigA_enrichment),
#                                                                                        "sd_enrichment" = sd(sigA_enrichment))
#pvalue_clinical_sigA <- average_mut_num %>% group_by(group) %>%
#  summarise("pvalue" = calculate_pvalue(decile, sigA_enrichment))
#
#ggplot(average_mut_num_summary_sigA, aes(x = decile, y = average_enrichment, group = group, color = group)) +
#  geom_line(size=2) +
#  geom_point(size=2, colour="black") +
#  scale_color_manual(values = cols) +
#  geom_errorbar(aes(ymin = average_enrichment-sd_enrichment, ymax = average_enrichment+sd_enrichment), width = 0.2) +
#  ylim(c(0, 2)) +
#  theme_classic() +
#  labs(x = "Gene expression levels",
#       y = "Mutation enrichment ratio \n (observed/expected)",
#       title = "Signature A")
#ggsave(paste0(workdir,"signatureA_enrich_analysis.pdf"), height = 4, width = 6)
#
#############################################################
#average_mut_num_summary_sigC <- average_mut_num_rescale %>% group_by(decile, group) %>% summarise("average_enrichment" = mean(sigC_enrichment),
#                                                                                        "sd_enrichment" = sd(sigC_enrichment))
#pvalue_clinical_sigC <- average_mut_num %>% group_by(group) %>%
#  summarise("pvalue" = calculate_pvalue(decile, sigC_enrichment))
#
#ggplot(average_mut_num_summary_sigC, aes(x = decile, y = average_enrichment, group = group, color = group)) +
#  geom_line(size=2) +
#  geom_point(size=2, colour="black") +
#  scale_color_manual(values = cols) +
#  geom_errorbar(aes(ymin = average_enrichment-sd_enrichment, ymax = average_enrichment+sd_enrichment), width = 0.2) +
#  theme_classic() +
#  labs(x = "Gene expression levels",
#       y = "Mutation enrichment ratio \n (observed/expected)",
#       title = "Signature C")
#ggsave(paste0(workdir,"signatureC_enrich_analysis.pdf"), height = 4, width = 6)




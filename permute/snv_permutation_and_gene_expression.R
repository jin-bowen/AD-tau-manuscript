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

cols <- c("Tau" = "#D62A28", "noTau" = "#D87AB2")
workdir='/n/data1/bwh/pathology/miller/lab/AD_Tau/results_hg19/'
metafile = '~/AD_Tau_PTA/metafiles/all_meta.tab'
snv_mutation_file='/n/data1/bwh/pathology/miller/lab/AD_Tau/results_hg19/annovar/snv_list.annovar.variant_function'
indel_mutation_file='/n/data1/bwh/pathology/miller/lab/AD_Tau/results_hg19/annovar/snv_list.annovar.variant_function'
permutation_tau_file='/n/data1/bwh/pathology/miller/lab/AD_Tau/results_hg19/Permute_Tau/annovar/perms_by_cell_snv.combined.variant_function'
permutation_notau_file='/n/data1/bwh/pathology/miller/lab/AD_Tau/results_hg19/Permute_noTau/annovar/perms_by_cell_snv.combined.variant_function'
group_num = 10

# read metafile
metadata_full <- read.table(metafile,sep=",", header = T,  colClasses = c("character", "character", "character"))
metadata  = metadata_full[metadata_full$group %in% c('Tau','noTau'),]
head(metadata)

library(tidyr)
# read in snv data
snv_mutation <- read.table(snv_mutation_file)
colnames(snv_mutation) <- c("region", "gene", "chr", "start", "end", "ref", "alt", "donor", "single_cell_ID")
snv_mutation = snv_mutation %>% separate(donor, c('donor', 'batch'))

genic_snv_mutation <- snv_mutation[snv_mutation$region %in% c("exonic", "exonic;splicing",
                                                  "intronic", "splicing",
                                                  "UTR3", "UTR5", "UTR5;UTR3"),]

genic_snv_mutation$gene <- str_remove(genic_snv_mutation$gene, "\\(.*\\)$") # remove mutations corresponding to more than one genes
genic_snv_mutation <- genic_snv_mutation[!str_detect(genic_snv_mutation$gene, ","), ]
head(genic_snv_mutation)

library(dplyr)
genic_snv_mutation_meta = dplyr::left_join(genic_snv_mutation, metadata, by=c("single_cell_ID" = "sample"), suffix = c(".x", ""),)
nrow(genic_snv_mutation_meta)
nrow(unique(genic_snv_mutation_meta))

snv_mutation_num <- data.frame(table(genic_snv_mutation_meta$gene, genic_snv_mutation_meta$group))
colnames(snv_mutation_num) <- c("gene", "group", "mut_number")

rdsfilename=paste0(workdir,'Tau_noTau_snv_mutation.rds')
saveRDS(snv_mutation_num, rdsfilename)
snv_mutation_num=readRDS(paste0(workdir,'Tau_noTau_snv_mutation.rds'))


snv_mutation_perm = data.frame()
# read in permutation data
genic_permutation_all <- c()
for (group in c('Tau','noTau')){

	if (group=='Tau'){
		permutation_file= permutation_tau_file
	}else if (group=='noTau'){
		permutation_file= permutation_notau_file
	}
	
	permutation <- read.table(permutation_file, header = F)
	permutation <- permutation[,c(1:7, 15)]
	colnames(permutation) <- c("region", "gene", "chr", "start", "end", "ref", "alt", "perm.id")
	genic_permutation <- permutation[permutation$region %in% c("exonic", "exonic;splicing",
	                                                           "intronic", "splicing",
	                                                           "UTR3", "UTR5", "UTR5;UTR3"),]

		
	genic_permutation$gene <- str_remove(genic_permutation$gene, "\\(.*\\)$") 
	genic_permutation <- genic_permutation[!str_detect(genic_permutation$gene, ","), ]
	genic_permutation$group = group
	genic_permutation_all = rbind(genic_permutation_all, genic_permutation)
	
	length(levels(as.factor(genic_permutation$gene)))
	
	permutation_num <- data.frame(table(genic_permutation$gene, genic_permutation$perm.id))
	colnames(permutation_num) <- c("gene", "perm.id", "permutation_number")
	
	rdsfilename=paste0(workdir,group,'_snv_permutation.rds')
	saveRDS(permutation_num, rdsfilename)
}	
rdsfilename=paste0(workdir,'genic_permutation.rds')
saveRDS(genic_permutation_all, rdsfilename)

permutation_tau_num=readRDS("/n/data1/bwh/pathology/miller/lab/AD_Tau/results_hg19/Tau_snv_permutation.rds")
permutation_notau_num=readRDS("/n/data1/bwh/pathology/miller/lab/AD_Tau/results_hg19/noTau_snv_permutation.rds")

snv_mutation_perm <- c()
for (group in c('Tau','noTau')){
	if (group=='Tau'){
		permutation_num= permutation_tau_num
	}else if (group=='noTau'){
		permutation_num= permutation_notau_num
	}
		snv_mutation_num_sub = snv_mutation_num[snv_mutation_num$group==group,]
	snv_mutation_perm_sub = right_join(snv_mutation_num_sub, permutation_num, by = c("gene"))
	snv_mutation_perm_sub$clinic = group
	snv_mutation_perm = rbind(snv_mutation_perm,snv_mutation_perm_sub)
}


snv_mutation_perm$mut_number[is.na(snv_mutation_perm$mut_number)] <- 0
snv_mutation_perm$permutation_number[is.na(snv_mutation_perm$permutation_number)] <- 0
snv_mutation_perm$group=snv_mutation_perm$clinic
rdsfilename=paste0(workdir,'snv_mutation_perm.rds')
saveRDS(snv_mutation_perm, rdsfilename)

#########################################################################################
group_num = 10
gene_expression_dir='/n/data1/bwh/pathology/miller/lab/AD_Tau/results_hg19/'
cte_neuron_data <- readRDS("/n/data1/bch/genetics/lee/shulin/CTE/results/CTE_scRNAseq/inhouse_data_neurons.rds")
cte_expr_level <- data.frame(Seurat::AverageExpression(cte_neuron_data, group.by = "condition", slot = "data")$RNA)
cte_expr_level["gene"] <- row.names(cte_expr_level)
colnames(cte_expr_level) <- c("average_expr_level_control", "average_expr_level_CTE", "gene")
cte_expr_level["decile"] <- ntile(cte_expr_level$average_expr_level_control, n = group_num)
cte_expr_level$decile <- as.factor(cte_expr_level$decile)
rdsfilename=paste0(gene_expression_dir,'cte_expr_level.rds')
saveRDS(cte_expr_level, rdsfilename)

############## expression by Tau status
group_num = 10
neuron <-  readRDS('/n/data1/bwh/pathology/miller/lab/AD_Tau/results_hg19/single-soma_Excitatory.rds')
gene_table = read.csv('/n/data1/bwh/pathology/miller/lab/ref/gene_names_table.csv', header=T)
expr_level <- data.frame(Seurat::AverageExpression(neuron,  slot = "data")$RNA)
expr_level["gene"] <- row.names(expr_level)
colnames(expr_level) <- c("average_expr_level", "gene")
expr_level = inner_join(expr_level,gene_table, by=c("gene" = "gene"))
expr_level["decile"] <- ntile(expr_level$average_expr_level, n = group_num)
expr_level$decile <- as.factor(expr_level$decile)

Tau_expr_level =   data.frame(Seurat::AverageExpression(neuron,group.by = "SORT", slot = "data")$RNA)
Tau_expr_level["gene"] = row.names(Tau_expr_level)
Tau_expr_level["AT8_decile"] = ntile(Tau_expr_level$AT8, n=group_num)
Tau_expr_level$AT8_decile <- as.factor(Tau_expr_level$AT8_decile)

Tau_expr_level["MAP2_decile"] = ntile(Tau_expr_level$MAP2, n=group_num)
Tau_expr_level$MAP2_decile <- as.factor(Tau_expr_level$MAP2_decile)


All_expr_level = inner_join(Tau_expr_level,expr_level, by=c("gene" = "gene"))
rdsfilename=paste0(gene_expression_dir,'single-soma_expr_level.rds')
saveRDS(All_expr_level, rdsfilename)

#################################gtex expression####################
group_num = 10
gtex_neuron_data <- read.table("/home/boj924/AD_Tau_PTA/metafiles/gene_tpm_2017-06-05_v8_brain_frontal_cortex_ba9.gct", header=T)
average_expression = rowMeans(select(gtex_neuron_data,contains("GTEX")))

gtex_expr_level <- gtex_neuron_data$Description
gtex_expr_level <- as.data.frame(gtex_expr_level)
gtex_expr_level["average_expr_level"] <- average_expression
gtex_expr_level["decile"] <- ntile(gtex_expr_level$average_expr_level, n = group_num)
gtex_expr_level$decile <- as.factor(gtex_expr_level$decile)
gtex_expr_level$gene <- gtex_neuron_data$Description
gtex_expr_level = gtex_expr_level[!duplicated(gtex_expr_level[ , c("gene")]),]

rdsfilename=paste0(gene_expression_dir,'gtex_expr_level.rds')
saveRDS(gtex_expr_level, rdsfilename)


#########################################################################################
expr_level=readRDS(paste0(workdir, 'gtex_expr_level.rds'))
expr_level=readRDS(paste0(workdir, 'cte_expr_level.rds'))

expr_level_mutation <- inner_join(expr_level, snv_mutation_perm, by = c("gene_symbol" = "gene"))
expr_level_mutation$group = expr_level_mutation$clinic

average_mut_num <- expr_level_mutation %>% group_by(decile, group, perm.id) %>% summarise("mutation_number" = sum(mut_number),
                                                                                   "permutation_number" = sum(permutation_number))
average_mut_num$enrichment_ratio <- average_mut_num$mutation_number/average_mut_num$permutation_number
rdsfilename=paste0(workdir, 'snv_average_mut_num.rds')
saveRDS(average_mut_num, rdsfilename)

average_mut_num_summary <- average_mut_num %>% group_by(decile, group) %>% summarise("average_enrichment" = mean(enrichment_ratio),
                                                             "sd_enrichment" = sd(enrichment_ratio),
                                                             "mutation_num" = sum(mutation_number),
                                                             "permutation_num" = sum(permutation_number))

ggplot(average_mut_num_summary, aes(x = decile, y = average_enrichment, group = group, color = group)) +
  geom_line(size=2) +
  geom_point(size=2, colour="black") +
  scale_color_manual(values = cols) +
  geom_errorbar(aes(ymin = average_enrichment-sd_enrichment, ymax = average_enrichment+sd_enrichment), width = 0.2) +
  theme_classic() +
  ylim(c(0.5, 1.5)) +
  labs(x = "Gene expression levels",
       y = "Mutation enrichment ratio \n (observed/expected)",
       title = "Total sSNVs")
ggsave(paste0(workdir,"snv_enrichment_expr.pdf"), height = 4, width = 6)

################################################################################################
expr_level=readRDS(paste0(workdir, 'single-soma_expr_level.rds'))

Tau_snv_mutation_perm = snv_mutation_perm[snv_mutation_perm$group=='Tau',]
Tau_expr_level_mutation <- inner_join(expr_level,Tau_snv_mutation_perm, by = c("gene" = "gene"))
Tau_average_mut_num <- Tau_expr_level_mutation  %>% group_by(AT8_decile, group, perm.id) %>% summarise("mutation_number" = sum(mut_number),
                                                                                   "permutation_number" = sum(permutation_number))
Tau_average_mut_num$enrichment_ratio <- Tau_average_mut_num$mutation_number/Tau_average_mut_num$permutation_number
Tau_average_mut_num_summary <- Tau_average_mut_num %>% group_by(AT8_decile,group) %>% summarise("average_enrichment" = mean(enrichment_ratio),
                                                             "sd_enrichment" = sd(enrichment_ratio),
                                                             "mutation_num" = sum(mutation_number),
                                                             "permutation_num" = sum(permutation_number))
Tau_average_mut_num$decile = Tau_average_mut_num$AT8_decile
Tau_average_mut_num_summary$decile = Tau_average_mut_num_summary$AT8_decile

noTau_snv_mutation_perm = snv_mutation_perm[snv_mutation_perm$group=='noTau',]
noTau_expr_level_mutation <- inner_join(expr_level,noTau_snv_mutation_perm, by = c("gene_symbol" = "gene"))
noTau_average_mut_num <- noTau_expr_level_mutation  %>% group_by(MAP2_decile, group, perm.id) %>% summarise("mutation_number" = sum(mut_number),
                                                                                   "permutation_number" = sum(permutation_number))
noTau_average_mut_num$enrichment_ratio <- noTau_average_mut_num$mutation_number/noTau_average_mut_num$permutation_number
noTau_average_mut_num_summary <- noTau_average_mut_num %>% group_by(MAP2_decile,group) %>% summarise("average_enrichment" = mean(enrichment_ratio),
                                                             "sd_enrichment" = sd(enrichment_ratio),
                                                             "mutation_num" = sum(mutation_number),
                                                             "permutation_num" = sum(permutation_number))
noTau_average_mut_num$decile = noTau_average_mut_num$MAP2_decile
noTau_average_mut_num_summary$decile = noTau_average_mut_num_summary$MAP2_decile
All_average_mut_num_summary = rbind(Tau_average_mut_num_summary,noTau_average_mut_num_summary)                                                    
All_average_mut_num = rbind(Tau_average_mut_num,noTau_average_mut_num)
rdsfilename=paste0(workdir, 'All_average_mut_num.rds')
saveRDS(All_average_mut_num, rdsfilename)


ggplot(All_average_mut_num_summary, aes(x = decile, y = average_enrichment, group = group, color = group)) +
  geom_line(size=2) +
  geom_point(size=2, colour="black") +
  scale_color_manual(values = cols) +
  geom_errorbar(aes(ymin = average_enrichment-sd_enrichment, ymax = average_enrichment+sd_enrichment), width = 0.2) +
  theme_classic() +
  ylim(c(0, 2)) +
  labs(x = "Gene expression levels",
       y = "Mutation enrichment ratio \n (observed/expected)",
       title = "Total sSNVs")
ggsave(paste0(workdir,"Tau_snv_enrichment_expr.pdf"), height = 4, width = 6)

############################################################################################
ref_genome="BSgenome.Hsapiens.UCSC.hg19"
library(ref_genome, character.only = T)
chr_orders <- c(paste0("chr", 1:22), "chrX", "chrY", "chrM")
library(MutationalPatterns)
library(GenomicRanges)

average_mut_num = readRDS(paste0(workdir, 'All_average_mut_num.rds'))
average_mut_num["sigA_mutation"] <- NA
average_mut_num["sigC_mutation"] <- NA

for (group in c('Tau','noTau')){
	decile_grange_list <- GRangesList()
	for (i in 1:group_num){
		if (group == "Tau"){
			gene_name <- unique(Tau_expr_level_mutation[Tau_expr_level_mutation$AT8_decile == i, "gene_symbol"])
			genic_mutation_decile <- genic_snv_mutation_meta[(genic_snv_mutation_meta$gene %in% gene_name) & genic_snv_mutation_meta$group == "Tau",]
		} else if (group == "noTau" ) {
			gene_name <- unique(noTau_expr_level_mutation[noTau_expr_level_mutation$MAP2_decile == i, "gene_symbol"])
			genic_mutation_decile <- genic_snv_mutation_meta[(genic_snv_mutation_meta$gene %in% gene_name) & genic_snv_mutation_meta$group == "noTau",]
		}

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
	if (group == "Tau"){
		average_mut_num$sigA_mutation[average_mut_num$group == "Tau"] <- aging_signature_contribution[1, unlist(lapply(1:group_num, rep, 1000))]
		average_mut_num$sigC_mutation[average_mut_num$group == "Tau"] <- aging_signature_contribution[2, unlist(lapply(1:group_num, rep, 1000))]
	} else if(group == "noTau"){
		average_mut_num$sigA_mutation[average_mut_num$group == "noTau"] <- aging_signature_contribution[1, unlist(lapply(1:group_num, rep, 1000))]
		average_mut_num$sigC_mutation[average_mut_num$group == "noTau"] <- aging_signature_contribution[2, unlist(lapply(1:group_num, rep, 1000))]
	}
}


#############################################################################
genic_permutation = readRDS(paste0(workdir,'genic_permutation.rds'))
aging_signature_contribution_summary <- c()
for (group in c('Tau','noTau')){ 
	for (perm_round in 1:1000){
	decile_grange_list <- GRangesList()
	for (i in 1:group_num){
		if (group == "Tau"){
			gene_name <- unique(Tau_expr_level_mutation[Tau_expr_level_mutation$AT8_decile == i, "gene_symbol"])
			genic_permutation_decile <- genic_permutation[(genic_permutation$gene %in% gene_name) & (genic_permutation$group == "Tau") &
									(genic_permutation$perm.id == perm_round), ]
		} else if (group == "noTau" ) {
			gene_name <- unique(noTau_expr_level_mutation[noTau_expr_level_mutation$MAP2_decile == i, "gene_symbol"])
			genic_permutation_decile <- genic_permutation[(genic_permutation$gene %in% gene_name) & (genic_permutation$group == "noTau") & 
									(genic_permutation$perm.id == perm_round),]
		}

		genic_permutation_decile <- genic_permutation_decile[order(genic_permutation_decile$start),]
		genic_permutation_decile <- genic_permutation_decile[order(genic_permutation_decile$chr),]
		decile_grange_list[[paste0("decile_", i)]] <- GRanges(seqnames = genic_permutation_decile$chr,
					ranges = IRanges(start = genic_permutation_decile$start,
							end = genic_permutation_decile$end),
							ref = genic_permutation_decile$ref,
							alt = genic_permutation_decile$alt)
		print(i)
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
	row.names(aging_signature_contribution) <- 1:group_num
	colnames(aging_signature_contribution) <- c("sigA_permutation", "sigC_permutation")
	aging_signature_contribution[,"perm.id"] <- perm_round	
	aging_signature_contribution[,"decile"] <- 1:group_num
	aging_signature_contribution[,"group"] <- group

	aging_signature_contribution_summary <- rbind(aging_signature_contribution_summary, aging_signature_contribution)
	if(perm_round %% 100 == 0){print(perm_round)}
}}

rdsfilename=paste0(workdir,'snv_permutation_sig.rds')
saveRDS(aging_signature_contribution_summary, rdsfilename)
aging_signature_contribution_summary = readRDS(paste0(workdir,'snv_permutation_sig.rds'))


###########################################################################
rdsfilename=paste0(workdir,'snv_permutation_sig.rds')
saveRDS(average_mut_num, rdsfilename)
average_mut_num <- readRDS(rdsfilename)
average_mut_num$AT8_decile_rescale = ceiling(as.numeric(average_mut_num$AT8_decile)/2)
average_mut_num$AT8_decile_rescale = as.factor(average_mut_num$AT8_decile_rescale)

average_mut_num$MAP2_decile_rescale = ceiling(as.numeric(average_mut_num$MAP2_decile)/2)
average_mut_num$MAP2_decile_rescale = as.factor(average_mut_num$MAP2_decile_rescale)

select = average_mut_num$group=='Tau'
average_mut_num[select,'decile'] = average_mut_num[select, 'AT8_decile_rescale']

select = average_mut_num$group=='noTau'
average_mut_num[select,'decile'] = average_mut_num[select, 'MAP2_decile_rescale']

average_mut_num$decile <- as.factor(average_mut_num$decile)
average_mut_num$perm.id <- as.integer(as.character(average_mut_num$perm.id))

average_mut_num_rescale <- average_mut_num %>% group_by(decile, group, perm.id) %>% summarise( 
							mutation_number = sum(mutation_number),
							sigA_mutation = sum(sigA_mutation),
							sigC_mutation = sum(sigC_mutation),
							sigA_permutation = sum(sigA_permutation),
							sigC_permutation = sum(sigC_permutation),
							permutation_number = sum(permutation_number))


average_mut_num_rescale["sigA_enrichment"] <-
  (average_mut_num_rescale$sigA_mutation * average_mut_num_rescale$mutation_number) / (average_mut_num_rescale$sigA_permutation * average_mut_num_rescale$permutation_number)
average_mut_num_rescale[is.infinite(average_mut_num_rescale$sigA_enrichment), "sigA_enrichment"] <- 1
average_mut_num_rescale["sigC_enrichment"] <-
  (average_mut_num_rescale$sigC_mutation * average_mut_num_rescale$mutation_number) / (average_mut_num_rescale$sigC_permutation * average_mut_num_rescale$permutation_number)
average_mut_num_rescale[is.infinite(average_mut_num_rescale$sigC_enrichment), "sigC_enrichment"] <- 1



##############################################################
average_mut_num_summary_sigA <- average_mut_num_rescale %>% group_by(decile, group) %>% summarise("average_enrichment" = mean(sigA_enrichment),
                                                                                        "sd_enrichment" = sd(sigA_enrichment))
pvalue_clinical_sigA <- average_mut_num %>% group_by(group) %>%
  summarise("pvalue" = calculate_pvalue(decile, sigA_enrichment))

ggplot(average_mut_num_summary_sigA, aes(x = decile, y = average_enrichment, group = group, color = group)) +
  geom_line(size=2) +
  geom_point(size=2, colour="black") +
  scale_color_manual(values = cols) +
  geom_errorbar(aes(ymin = average_enrichment-sd_enrichment, ymax = average_enrichment+sd_enrichment), width = 0.2) +
  ylim(c(0, 2)) +
  theme_classic() +
  labs(x = "Gene expression levels",
       y = "Mutation enrichment ratio \n (observed/expected)",
       title = "Signature A")
ggsave(paste0(workdir,"signatureA_enrich_analysis.pdf"), height = 4, width = 6)

############################################################
average_mut_num_summary_sigC <- average_mut_num_rescale %>% group_by(decile, group) %>% summarise("average_enrichment" = mean(sigC_enrichment),
                                                                                        "sd_enrichment" = sd(sigC_enrichment))
pvalue_clinical_sigC <- average_mut_num %>% group_by(group) %>%
  summarise("pvalue" = calculate_pvalue(decile, sigC_enrichment))

ggplot(average_mut_num_summary_sigC, aes(x = decile, y = average_enrichment, group = group, color = group)) +
  geom_line(size=2) +
  geom_point(size=2, colour="black") +
  scale_color_manual(values = cols) +
  geom_errorbar(aes(ymin = average_enrichment-sd_enrichment, ymax = average_enrichment+sd_enrichment), width = 0.2) +
  theme_classic() +
  labs(x = "Gene expression levels",
       y = "Mutation enrichment ratio \n (observed/expected)",
       title = "Signature C")
ggsave(paste0(workdir,"signatureC_enrich_analysis.pdf"), height = 4, width = 6)




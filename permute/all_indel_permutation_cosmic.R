library(ggplot2)
library(stringr)
library(dplyr)
library(reshape2)
library(tidyverse)
library(tidyr)

cols <- c("Tau" = "#D62728", "noTau" = "#FFC0CB")
workdir='/n/data1/bwh/pathology/miller/lab/AD_Tau_PTA_results/'
metafile = '~/AD_Tau_PTA/metafiles/all_meta.tab'
indel_mutation_file='/n/data1/bwh/pathology/miller/lab/AD_Tau_PTA_results/indel_list.annovar.variant_function'
indel_permutation_file='/n/data1/bwh/pathology/miller/lab/AD_Tau_PTA_results/all_perms_by_cell_indel.combined.variant_function'
group_num = 10

# read metafile
metadata <- read.table(metafile,sep=",", header = T,  colClasses = c("character", "character", "character"))

library(tidyr)
# read in indel data
indel_mutation <- read.table(indel_mutation_file)
indel_mutation <- indel_mutation[,c(1:5, 11,12,15)]
colnames(indel_mutation) <- c("region", "gene", "chr", "start", "end", "ref", "alt", "sample_infor")
indel_mutation = indel_mutation %>% separate(sample_infor, sep=',', c('donor', 'single_cell_ID')) %>% separate(donor, sep='_', c('donor', 'batch'))

genic_indel_mutation <- indel_mutation[indel_mutation$region %in% c("exonic", "exonic;splicing",
                                                  "intronic", "splicing","UTR3", "UTR5", "UTR5;UTR3"),]

genic_indel_mutation$gene <- str_remove(genic_indel_mutation$gene, "\\(.*\\)$") # remove mutations corresponding to more than one genes
genic_indel_mutation <- genic_indel_mutation[!str_detect(genic_indel_mutation$gene, ","), ]
head(genic_indel_mutation)

library(dplyr)
genic_indel_mutation_meta = dplyr::left_join(genic_indel_mutation, metadata, by=c("single_cell_ID" = "sample"), suffix = c(".x", ""),)

indel_mutation_num <- data.frame(table(genic_indel_mutation_meta$gene))
colnames(indel_mutation_num) <- c("gene", "mut_number")

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
#
#
#average_mut_num_summary <- average_mut_num %>% group_by(decile) %>% summarise("average_enrichment" = mean(enrichment_ratio),
#                                                             "sd_enrichment" = sd(enrichment_ratio),
#                                                             "mutation_num" = sum(mutation_number),
#                                                             "permutation_num" = sum(permutation_number))
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

###############################################################################################
ref_genome="BSgenome.Hsapiens.UCSC.hg19"
library(ref_genome, character.only = T)
library(MutationalPatterns)
library(GenomicRanges)
chr_orders <- c(paste0("chr", 1:22), "chrX", "chrY", "chrM")

average_mut_num = readRDS(paste0(workdir, 'indel_average_mut_num.rds'))
average_mut_num["ID83A_mutation"] <- NA

decile_grange_list <- GRangesList()
for (i in 1:group_num){
        gene_name <- unique(expr_level_mutation[expr_level_mutation$decile == i, "gene"])
        genic_mutation_decile <- genic_indel_mutation[(genic_indel_mutation$gene %in% gene_name),]

        genic_mutation_decile <- genic_mutation_decile[order(genic_mutation_decile$start),]
        genic_mutation_decile <- genic_mutation_decile[order(genic_mutation_decile$chr),]

        decile_grange_list[[paste0('decile_',i)]] <- GRanges(seqnames = paste0("chr", genic_mutation_decile$chr),
                                ranges = IRanges(start = genic_mutation_decile$start,
                                                end = genic_mutation_decile$end),
                                                ref = genic_mutation_decile$ref,
                                                alt = genic_mutation_decile$alt)
}
chr_length <- seqlengths(Hsapiens)
seqlengths(decile_grange_list) <- chr_length[names(seqlengths(decile_grange_list))]
seqlevels(decile_grange_list) <- seqlevels(decile_grange_list)[order(factor(seqlevels(decile_grange_list), levels = chr_orders))]
genome(decile_grange_list) = 'hg19'
indel_grl <- get_mut_type(decile_grange_list, type = "indel")
indel_grl <- get_indel_context(indel_grl, ref_genome =ref_genome)
mut_mat = count_indel_contexts(indel_grl)
aging_signature_all <- readRDS("/home/boj924/AD_Tau_PTA/results/indel_cosmic_sig.rds")
mut_mat = mut_mat[rownames(aging_signature_all),]
aging_signature_fit <- fit_to_signatures(mut_mat, as.matrix(aging_signature_all))
aging_signature_contribution <- apply(aging_signature_fit$contribution, 2, function(x){x/sum(x)})
colnames(aging_signature_contribution) <- gsub("decile_", "", colnames(aging_signature_contribution))
average_mut_num$ID83A_mutation <- aging_signature_contribution[6, unlist(lapply(colnames(aging_signature_contribution), rep, 1000))]

rdsfilename=paste0(workdir,'indel_average_mut_num.rds')
saveRDS(average_mut_num, rdsfilename)

###############################################################################################
average_mut_num = readRDS(paste0(workdir, 'indel_average_mut_num.rds'))
average_mut_num["ID83A_mutation"] <- NA

decile_grange_list <- GRangesList()
for (i in 1:group_num){
        gene_name <- unique(expr_level_mutation[expr_level_mutation$decile == i, "gene"])
        genic_mutation_decile <- genic_indel_mutation[(genic_indel_mutation$gene %in% gene_name),]

        genic_mutation_decile <- genic_mutation_decile[order(genic_mutation_decile$start),]
        genic_mutation_decile <- genic_mutation_decile[order(genic_mutation_decile$chr),]

        decile_grange_list[[paste0('decile_',i)]] <- GRanges(seqnames = paste0("chr", genic_mutation_decile$chr),
                                ranges = IRanges(start = genic_mutation_decile$start,
                                                end = genic_mutation_decile$end),
                                                ref = genic_mutation_decile$ref,
                                                alt = genic_mutation_decile$alt)
}
chr_length <- seqlengths(Hsapiens)
seqlengths(decile_grange_list) <- chr_length[names(seqlengths(decile_grange_list))]
seqlevels(decile_grange_list) <- seqlevels(decile_grange_list)[order(factor(seqlevels(decile_grange_list), levels = chr_orders))]
genome(decile_grange_list) = 'hg19'
indel_grl <- get_mut_type(decile_grange_list, type = "indel")
indel_grl <- get_indel_context(indel_grl, ref_genome =ref_genome)
mut_mat = count_indel_contexts(indel_grl)
aging_signature_all <- readRDS("/home/boj924/AD_Tau_PTA/results/indel_cosmic_sig.rds")
mut_mat = mut_mat[rownames(aging_signature_all),]
aging_signature_fit <- fit_to_signatures(mut_mat, as.matrix(aging_signature_all))
aging_signature_contribution <- apply(aging_signature_fit$contribution, 2, function(x){x/sum(x)})
colnames(aging_signature_contribution) <- gsub("decile_", "", colnames(aging_signature_contribution))
average_mut_num$ID83A_mutation <- aging_signature_contribution[6, unlist(lapply(colnames(aging_signature_contribution), rep, 1000))]

rdsfilename=paste0(workdir,'indel_average_mut_num.rds')
saveRDS(average_mut_num, rdsfilename)


###############################################################################################
genic_permutation = readRDS(paste0(workdir,'genic_indel_mutation.rds'))
aging_signature_contribution_summary <- c()
for (perm_round in 1:1000){
       decile_mut_list <- data.frame()
       for (i in 1:group_num){
               gene_name <- unique(expr_level_mutation[expr_level_mutation$decile == i, "gene"])
               genic_permutation_decile <- genic_permutation[(genic_permutation$gene %in% gene_name)  &
                                                               (genic_permutation$perm.id == perm_round), ]
               genic_permutation_decile <- genic_permutation_decile[order(genic_permutation_decile$start),]
               genic_permutation_decile <- genic_permutation_decile[order(genic_permutation_decile$chr),]
	       decile_mut_list_temp = as.data.frame(genic_permutation_decile$muttype, col.names=c('muttype'))
	       decile_mut_list_temp[,"decile"] = i
               decile_mut_list = rbind(decile_mut_list, decile_mut_list_temp)
       }
       mut_mat = as.data.frame(table(decile_mut_list$"genic_permutation_decile$muttype", decile_mut_list$decile))
       mut_mat_reform <- reshape(mut_mat, idvar='Var1', timevar='Var2',direction = "wide")
       colnames(mut_mat_reform) <- gsub("Freq.", "", colnames(mut_mat_reform))
       meta_file=read.delim('/home/boj924/AD_Tau_PTA/custome_sig/indel83_lookup.tab', sep=',')
       mut_mat_meta = as.data.frame(merge(mut_mat_reform, meta_file, by.x='Var1', by.y='indeltype1', all.y=T))  %>% replace(is.na(.), 0)
       rownames(mut_mat_meta) = mut_mat_meta$indeltype2
       mut_mat_final <- subset(mut_mat_meta, select = -c(Var1, indeltype2))
       aging_signature_all <- readRDS("/home/boj924/AD_Tau_PTA/results/indel_cosmic_sig.rds")
       mut_mat_final = mut_mat_final[rownames(aging_signature_all),]

       aging_signature_fit <- fit_to_signatures(mut_mat_final, as.matrix(aging_signature_all))
       aging_signature_contribution <- apply(aging_signature_fit$contribution, 2, function(x){x/sum(x)})

       aging_signature_contribution <- data.frame(t(aging_signature_contribution))
       colnames(aging_signature_contribution) <- c("ID83A_permutation")
       aging_signature_contribution[,"perm.id"] <- perm_round
       aging_signature_contribution[,"decile"] <- rownames(aging_signature_contribution)

       aging_signature_contribution_summary <- rbind(aging_signature_contribution_summary, aging_signature_contribution)
       if(perm_round %% 100 == 0){print(perm_round)}
}


rdsfilename=paste0(workdir,'indel_permutation_sig.rds')
saveRDS(aging_signature_contribution_summary, rdsfilename)


workdir='/n/data1/bwh/pathology/miller/lab/AD_Tau_PTA_results/'
mut_num = readRDS(paste0(workdir, 'indel_average_mut_num.rds'))

average_mut_num_summary <- mut_num %>% group_by(decile) %>% summarise("average_enrichment" = mean(enrichment_ratio),
                                                             "sd_enrichment" = sd(enrichment_ratio),
                                                             "mutation_num" = sum(mutation_number),
                                                             "permutation_num" = sum(permutation_number))
average_mut_num_summary$decile = as.numeric(average_mut_num_summary$decile)
average_mut_num_summary$group = 'AD'
model = lm(average_enrichment ~ decile, data = average_mut_num_summary)
summary(model)

p <- ggplot(average_mut_num_summary, aes(x = decile, y = average_enrichment, group = group, color = group)) +
#  geom_line(size=1,  linetype = 2) +
  geom_smooth(method=lm, se=F, size=1) +
  geom_point(size=2, colour="black") +
  scale_color_manual(values = cols) +
  geom_errorbar(aes(ymin = average_enrichment-sd_enrichment, ymax = average_enrichment+sd_enrichment), width = 0.2) +
  theme_classic() +
  ylim(c(0, 2)) +
  labs(x = "Gene expression levels",
       y = "Mutation enrichment ratio \n (observed/expected)",
       title = "Total sIndels")
ggsave(paste0(workdir,"indel_enrichment_sc.pdf"), plot=p, height = 4, width = 6)

################################################################################################
AD_aging_signature_contribution_summary <- readRDS(paste0(workdir,'indel_permutation_sig.rds'))
AD_aging_signature_contribution_summary$decile = as.factor(AD_aging_signature_contribution_summary$decile)
AD_aging_signature_contribution_summary$perm.id = as.integer(AD_aging_signature_contribution_summary$perm.id)
average_mut_num$decile <- as.factor(average_mut_num$decile)
average_mut_num$perm.id <- as.integer(as.character(average_mut_num$perm.id))

colnames(AD_aging_signature_contribution_summary)=c('ID3_permutation','ID4_permutation','ID5_permutation','ID8_permutation','ID11_permutation','ID83A_permutation','perm.id','decile')
AD_average_mut_num_meta <- left_join(average_mut_num, AD_aging_signature_contribution_summary, by = c("decile", "perm.id"))
AD_average_mut_num_rescale <- AD_average_mut_num_meta %>% group_by(decile, perm.id) %>% summarise( 
							mutation_number = sum(mutation_number),
							ID83A_mutation = sum(ID83A_mutation),
							ID83A_permutation = sum(ID83A_permutation),
							permutation_number = sum(permutation_number))

AD_average_mut_num_rescale["ID83A_enrichment"] <-
  (AD_average_mut_num_rescale$ID83A_mutation * AD_average_mut_num_rescale$mutation_number) / (AD_average_mut_num_rescale$ID83A_permutation * AD_average_mut_num_rescale$permutation_number)
AD_average_mut_num_rescale[is.infinite(AD_average_mut_num_rescale$ID83A_enrichment), "ID83A_enrichment"] <- 1

average_mut_num_summary_ID83A <- AD_average_mut_num_rescale %>% group_by(decile) %>% summarise("average_enrichment" = mean(ID83A_enrichment),
                                                                                        "sd_enrichment" = sd(ID83A_enrichment))



p <- ggplot(average_mut_num_summary_ID83A, aes(x = decile, y = average_enrichment)) +
  geom_smooth(method=lm, se=T, size=1, colour="red") +
#  geom_line(size=2, color='red') +
  geom_point(size=2, colour="red") +
  scale_color_manual(values = cols) +
  geom_errorbar(aes(ymin = average_enrichment-sd_enrichment, ymax = average_enrichment+sd_enrichment), width = 0.2) +
  theme_classic() +
  ylim(c(0.5,1.5))+
  labs(x = "Gene expression levels",
       y = "Mutation enrichment ratio \n (observed/expected)",
       title = "ID83A")
ggsave(paste0(workdir,"ID83A_enrichment_sc.pdf"), plot=p, height = 4, width = 5)
average_mut_num_summary_ID83A$decile = as.numeric(average_mut_num_summary_ID83A$decile)
model = lm(average_enrichment ~ decile,  data = average_mut_num_summary_ID83A)
summary(model)




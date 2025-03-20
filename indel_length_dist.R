library(ggplot2)
library(stringr)
library(dplyr)
library(reshape2)
library(tidyverse)
library(tidyr)


workdir='/n/data1/bwh/pathology/miller/lab/AD_Tau_PTA/results_hg19/'
metafile = '~/AD_Tau_PTA/metafiles/all_meta.tab'
indel_mutation_file='/n/data1/bwh/pathology/miller/lab/AD_Tau/results_hg19/annovar/indel_list.annovar.variant_function'
AD_indel_mutation_file='/n/data1/bwh/pathology/miller/lab/AD_hg19/annovar/indel_list.combined.variant_function'
ctrl_indel_mutation_file='/n/data1/bwh/pathology/miller/lab/ctrl_hg19/annovar/indel_list.annovar.variant_function'

# read metafile
metadata <- read.table(metafile,sep=",", header = T,  colClasses = c("character", "character", "character"))

# read in indel data
indel_mutation <- read.table(indel_mutation_file)
colnames(indel_mutation) <- c("region", "gene", "chr", "start", "end", "ref", "alt",
			"len","pos","id","ref_nt","alt_nt","Score","QC", "sample_infor")
indel_mutation = indel_mutation %>% separate(sample_infor, sep=',', c('donor', 'single_cell_ID'))
indel_mutation_meta = dplyr::inner_join(indel_mutation, metadata, by=c("single_cell_ID" = "sample"), suffix = c(".x", ""),)
indel_mutation_meta$length = nchar(indel_mutation_meta$alt_nt) - nchar(indel_mutation_meta$ref_nt)
indel_mutation_meta$length = ifelse( indel_mutation_meta$length >=5, "5", ifelse(indel_mutation_meta$length <= -5, "-5", indel_mutation_meta$length))


AD_indel_mutation <- read.table(AD_indel_mutation_file)
colnames(AD_indel_mutation) <- c("region", "gene", "chr", "start", "end", "ref", "alt",
			"len","pos","id","ref_nt","alt_nt","Score","QC", "sample_infor")
AD_indel_mutation = AD_indel_mutation %>% separate(sample_infor, sep=',', c('donor', 'single_cell_ID'))
AD_indel_mutation_meta = dplyr::inner_join(AD_indel_mutation, metadata, by=c("single_cell_ID" = "sample"), suffix = c(".x", ""),)
AD_indel_mutation_meta$length = nchar(AD_indel_mutation_meta$alt_nt) - nchar(AD_indel_mutation_meta$ref_nt)
AD_indel_mutation_meta$length = ifelse( AD_indel_mutation_meta$length >=5, "5", ifelse(AD_indel_mutation_meta$length <= -5, "-5", AD_indel_mutation_meta$length))

ctrl_indel_mutation <- read.table(ctrl_indel_mutation_file)
colnames(ctrl_indel_mutation) <- c("region", "gene", "chr", "start", "end", "ref", "alt",
			"len","pos","id","ref_nt","alt_nt","Score","QC", "sample_infor")
ctrl_indel_mutation = ctrl_indel_mutation %>% separate(sample_infor, sep=',', c('donor', 'single_cell_ID'))
ctrl_indel_mutation_meta = dplyr::inner_join(ctrl_indel_mutation, metadata, by=c("single_cell_ID" = "sample"), suffix = c(".x", ""),)
ctrl_indel_mutation_meta$length = nchar(ctrl_indel_mutation_meta$alt_nt) - nchar(ctrl_indel_mutation_meta$ref_nt)
ctrl_indel_mutation_meta$length = ifelse( ctrl_indel_mutation_meta$length >=5, "5", ifelse(ctrl_indel_mutation_meta$length <= -5, "-5", ctrl_indel_mutation_meta$length))


all_indel_mutation = rbind(indel_mutation_meta,AD_indel_mutation_meta)
all_indel_mutation = rbind(all_indel_mutation,ctrl_indel_mutation_meta)
all_indel_mutation$length = as.numeric(all_indel_mutation$length)
all_indel_mutation$group = factor(all_indel_mutation$group, levels=c('ctrl','AD','Tau','noTau'))

indel_table = all_indel_mutation %>% count(group, length) %>% group_by(group) %>%  mutate(percent = n / sum(n) * 100)

p <- indel_table %>% 
ggplot(aes(x = length,  y=percent, fill=group)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  xlim(-6, 6) +
  scale_fill_manual(values = c("Tau"="#e377c2", "noTau"="#d62728", 'AD'="#F57F20", 'ctrl'="#2278B5")) +  
  labs(x = "Indel size", 
       y = "Count",
       title = "sIndels size distribution") +
  facet_wrap(~group, ncol=1)

ggsave(paste0("/home/boj924/AD_Tau_PTA/figures/indel_length_distribution.pdf"), height = 5, width = 4)




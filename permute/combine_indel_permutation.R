library(ggplot2)
library(stringr)
library(dplyr)
library(reshape2)
library(tidyverse)
library(tidyr)


res_dir = '/home/boj924/AD_Tau_PTA/results/'
cols <- c("Tau" = "#D62A28", "noTau" = "#D87AB2", "AD" = "#F58020", "ctrl" = "#2579B6")
workdir_ctrl='/n/data1/bwh/pathology/miller/lab/ctrl_hg19/'
workdir_AD='/n/data1/bwh/pathology/miller/lab/AD_hg19/'
workdir_tau='/n/data1/bwh/pathology/miller/lab/AD_Tau/results_hg19/'
group_num = 10

expr_level=readRDS('/n/data1/bwh/pathology/miller/lab/AD_Tau/results_hg19/cte_expr_level.rds')
ctrl_mut_num = readRDS(paste0(workdir_ctrl, 'indel_average_mut_num.rds'))
AD_mut_num = readRDS(paste0(workdir_AD, 'indel_average_mut_num.rds'))
tau_mut_num = readRDS(paste0(workdir_tau, 'indel_average_mut_num.rds'))

ctrl_mut_num_summary <- ctrl_mut_num %>% group_by(decile) %>% summarise("average_enrichment" = mean(enrichment_ratio),
                                                             "sd_enrichment" = sd(enrichment_ratio),
                                                             "mutation_num" = sum(mutation_number),
                                                             "permutation_num" = sum(permutation_number))
AD_mut_num_summary <- AD_mut_num %>% group_by(decile) %>% summarise("average_enrichment" = mean(enrichment_ratio),
                                                             "sd_enrichment" = sd(enrichment_ratio),
                                                             "mutation_num" = sum(mutation_number),
                                                             "permutation_num" = sum(permutation_number))
tau_mut_num_summary <- tau_mut_num %>% group_by(decile, group) %>% summarise("average_enrichment" = mean(enrichment_ratio),
                                                             "sd_enrichment" = sd(enrichment_ratio),
                                                             "mutation_num" = sum(mutation_number),
                                                             "permutation_num" = sum(permutation_number))
ctrl_mut_num_summary$group = 'ctrl'
AD_mut_num_summary$group = 'AD'

average_mut_num_summary = rbind(ctrl_mut_num_summary, AD_mut_num_summary)
average_mut_num_summary = rbind(average_mut_num_summary, tau_mut_num_summary)

average_mut_num_summary$decile = as.numeric(average_mut_num_summary$decile)
model = lm(average_enrichment ~ decile * group, data = average_mut_num_summary)
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
ggsave(paste0(res_dir,"indel_enrichment_sc.pdf"), plot=p, height = 4, width = 6)

################################################################################################
tau_aging_signature_contribution_summary <- readRDS(paste0(workdir_tau,'indel_permutation_sig.rds'))
tau_aging_signature_contribution_summary$decile = as.factor(tau_aging_signature_contribution_summary$decile)
tau_aging_signature_contribution_summary$perm.id = as.integer(tau_aging_signature_contribution_summary$perm.id)
tau_mut_num$decile <- as.factor(tau_mut_num$decile)
tau_mut_num$perm.id <- as.integer(as.character(tau_mut_num$perm.id))

tau_average_mut_num_meta <- left_join(tau_mut_num, tau_aging_signature_contribution_summary, by = c("decile", "perm.id", "group"))
tau_average_mut_num_rescale <- tau_average_mut_num_meta %>% group_by(decile, group, perm.id) %>% summarise( 
							mutation_number = sum(mutation_number),
							ID4_mutation = sum(ID4_mutation),
							ID22_mutation = sum(ID22_mutation),
							ID4_permutation = sum(ID4_permutation),
							ID22_permutation = sum(ID22_permutation),
							permutation_number = sum(permutation_number))


tau_average_mut_num_rescale["ID4_enrichment"] <-
  (tau_average_mut_num_rescale$ID4_mutation * tau_average_mut_num_rescale$mutation_number) / (tau_average_mut_num_rescale$ID4_permutation * tau_average_mut_num_rescale$permutation_number)
tau_average_mut_num_rescale[is.infinite(tau_average_mut_num_rescale$ID4_enrichment), "ID4_enrichment"] <- 1
tau_average_mut_num_rescale["ID22_enrichment"] <-
  (tau_average_mut_num_rescale$ID22_mutation * tau_average_mut_num_rescale$mutation_number) / (tau_average_mut_num_rescale$ID22_permutation * tau_average_mut_num_rescale$permutation_number)
tau_average_mut_num_rescale[is.infinite(tau_average_mut_num_rescale$ID22_enrichment), "ID22_enrichment"] <- 1

tau_average_mut_num_summary_ID4 <- tau_average_mut_num_rescale %>% group_by(decile, group) %>% summarise("average_enrichment" = mean(ID4_enrichment),
                                                                                        "sd_enrichment" = sd(ID4_enrichment))
tau_average_mut_num_summary_ID22 <- tau_average_mut_num_rescale %>% group_by(decile, group) %>% summarise("average_enrichment" = mean(ID22_enrichment),
                                                                                        "sd_enrichment" = sd(ID22_enrichment))


############################################################
AD_aging_signature_contribution_summary <- readRDS(paste0(workdir_AD,'indel_permutation_sig.rds'))
AD_aging_signature_contribution_summary$decile = as.factor(AD_aging_signature_contribution_summary$decile)
AD_aging_signature_contribution_summary$perm.id = as.integer(AD_aging_signature_contribution_summary$perm.id)
AD_mut_num$decile <- as.factor(AD_mut_num$decile)
AD_mut_num$perm.id <- as.integer(as.character(AD_mut_num$perm.id))

AD_average_mut_num_meta <- left_join(AD_mut_num, AD_aging_signature_contribution_summary, by = c("decile", "perm.id"))
AD_average_mut_num_rescale <- AD_average_mut_num_meta %>% group_by(decile, perm.id) %>% summarise( 
							mutation_number = sum(mutation_number),
							ID4_mutation = sum(ID4_mutation),
							ID22_mutation = sum(ID22_mutation),
							ID4_permutation = sum(ID4_permutation),
							ID22_permutation = sum(ID22_permutation),
							permutation_number = sum(permutation_number))

AD_average_mut_num_rescale["ID4_enrichment"] <-
  (AD_average_mut_num_rescale$ID4_mutation * AD_average_mut_num_rescale$mutation_number) / (AD_average_mut_num_rescale$ID4_permutation * AD_average_mut_num_rescale$permutation_number)
AD_average_mut_num_rescale[is.infinite(AD_average_mut_num_rescale$ID4_enrichment), "ID4_enrichment"] <- 1
AD_average_mut_num_rescale["ID22_enrichment"] <-
  (AD_average_mut_num_rescale$ID22_mutation * AD_average_mut_num_rescale$mutation_number) / (AD_average_mut_num_rescale$ID22_permutation * AD_average_mut_num_rescale$permutation_number)
AD_average_mut_num_rescale[is.infinite(AD_average_mut_num_rescale$ID22_enrichment), "ID22_enrichment"] <- 1

AD_average_mut_num_summary_ID4 <- AD_average_mut_num_rescale %>% group_by(decile) %>% summarise("average_enrichment" = mean(ID4_enrichment),
                                                                                        "sd_enrichment" = sd(ID4_enrichment))

AD_average_mut_num_summary_ID22 <- AD_average_mut_num_rescale %>% group_by(decile) %>% summarise("average_enrichment" = mean(ID22_enrichment),
                                                                                        "sd_enrichment" = sd(ID22_enrichment))
AD_average_mut_num_summary_ID4$group='AD'
AD_average_mut_num_summary_ID22$group='AD'


average_mut_num_summary_ID4 = rbind(AD_average_mut_num_summary_ID4, tau_average_mut_num_summary_ID4)
average_mut_num_summary_ID22 = rbind(AD_average_mut_num_summary_ID22, tau_average_mut_num_summary_ID22)


p <- ggplot(average_mut_num_summary_ID4, aes(x = decile, y = average_enrichment, group = group, color = group)) +
  geom_line(size=2) +
  geom_point(size=2, colour="black") +
  scale_color_manual(values = cols) +
  geom_errorbar(aes(ymin = average_enrichment-sd_enrichment, ymax = average_enrichment+sd_enrichment), width = 0.2) +
  theme_classic() +
  ylim(c(0, 2)) +
  labs(x = "Gene expression levels",
       y = "Mutation enrichment ratio \n (observed/expected)",
       title = "ID4")
ggsave(paste0(res_dir,"ID4_enrichment_sc.pdf"), plot=p, height = 4, width = 6)
average_mut_num_summary_ID4$decile = as.numeric(average_mut_num_summary_ID4$decile)
model = lm(average_enrichment ~ decile,  data = average_mut_num_summary_ID4)
summary(model)


p <- ggplot(average_mut_num_summary_ID22, aes(x = decile, y = average_enrichment, group = group, color = group)) +
  geom_line(size=2) +
  geom_point(size=2, colour="black") +
  scale_color_manual(values = cols) +
  geom_errorbar(aes(ymin = average_enrichment-sd_enrichment, ymax = average_enrichment+sd_enrichment), width = 0.2) +
  theme_classic() +
  ylim(c(0, 2)) +
  labs(x = "Gene expression levels",
       y = "Mutation enrichment ratio \n (observed/expected)",
       title = "ID22")
ggsave(paste0(res_dir,"ID22_enrichment_sc.pdf"), plot=p, height = 4, width = 6)
average_mut_num_summary_ID22$decile = as.numeric(average_mut_num_summary_ID22$decile)
model = lm(average_enrichment ~ decile, data = average_mut_num_summary_ID22)
summary(model)




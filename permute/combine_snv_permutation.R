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
ctrl_mut_num = readRDS(paste0(workdir_ctrl, 'snv_average_mut_num.rds'))
AD_mut_num = readRDS(paste0(workdir_AD, 'snv_average_mut_num.rds'))
tau_mut_num = readRDS(paste0(workdir_tau, 'snv_average_mut_num.rds'))

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
model = lm(average_enrichment ~ decile + group, data = average_mut_num_summary)
summary(model)


p <- ggplot(average_mut_num_summary, aes(x = decile, y = average_enrichment, group = group, color = group)) +
  geom_line(size=2) +
  geom_point(size=2, colour="black") +
  scale_color_manual(values = cols) +
  geom_errorbar(aes(ymin = average_enrichment-sd_enrichment, ymax = average_enrichment+sd_enrichment), width = 0.2) +
  theme_classic() +
  ylim(c(0.5, 1.5)) +
  labs(x = "Gene expression levels",
       y = "Mutation enrichment ratio \n (observed/expected)",
       title = "Total sSNVs")
ggsave(paste0(res_dir,"snv_enrichment_sc.pdf"), plot=p, height = 4, width = 6)

################################################################################################
tau_aging_signature_contribution_summary <- readRDS(paste0(workdir_tau,'snv_permutation_sig.rds'))
tau_aging_signature_contribution_summary$decile = as.factor(tau_aging_signature_contribution_summary$decile)
tau_aging_signature_contribution_summary$perm.id = as.integer(tau_aging_signature_contribution_summary$perm.id)
tau_mut_num$decile <- as.factor(tau_mut_num$decile)
tau_mut_num$perm.id <- as.integer(as.character(tau_mut_num$perm.id))

tau_average_mut_num_meta <- left_join(tau_mut_num, tau_aging_signature_contribution_summary, by = c("decile", "perm.id", "group"))
tau_average_mut_num_rescale <- tau_average_mut_num_meta %>% group_by(decile, group, perm.id) %>% summarise( 
							mutation_number = sum(mutation_number),
							sigA_mutation = sum(sigA_mutation),
							sigC_mutation = sum(sigC_mutation),
							sigA_permutation = sum(sigA_permutation),
							sigC_permutation = sum(sigC_permutation),
							permutation_number = sum(permutation_number))


tau_average_mut_num_rescale["sigA_enrichment"] <-
  (tau_average_mut_num_rescale$sigA_mutation * tau_average_mut_num_rescale$mutation_number) / (tau_average_mut_num_rescale$sigA_permutation * tau_average_mut_num_rescale$permutation_number)
tau_average_mut_num_rescale[is.infinite(tau_average_mut_num_rescale$sigA_enrichment), "sigA_enrichment"] <- 1
tau_average_mut_num_rescale["sigC_enrichment"] <-
  (tau_average_mut_num_rescale$sigC_mutation * tau_average_mut_num_rescale$mutation_number) / (tau_average_mut_num_rescale$sigC_permutation * tau_average_mut_num_rescale$permutation_number)
tau_average_mut_num_rescale[is.infinite(tau_average_mut_num_rescale$sigC_enrichment), "sigC_enrichment"] <- 1

tau_average_mut_num_summary_sigA <- tau_average_mut_num_rescale %>% group_by(decile, group) %>% summarise("average_enrichment" = mean(sigA_enrichment),
                                                                                        "sd_enrichment" = sd(sigA_enrichment))
tau_average_mut_num_summary_sigC <- tau_average_mut_num_rescale %>% group_by(decile, group) %>% summarise("average_enrichment" = mean(sigC_enrichment),
                                                                                        "sd_enrichment" = sd(sigC_enrichment))


############################################################
AD_aging_signature_contribution_summary <- readRDS(paste0(workdir_AD,'AD_snv_permutation_sig.rds'))
AD_aging_signature_contribution_summary$decile = as.factor(AD_aging_signature_contribution_summary$decile)
AD_aging_signature_contribution_summary$perm.id = as.integer(AD_aging_signature_contribution_summary$perm.id)
AD_mut_num$decile <- as.factor(AD_mut_num$decile)
AD_mut_num$perm.id <- as.integer(as.character(AD_mut_num$perm.id))

AD_average_mut_num_meta <- left_join(AD_mut_num, AD_aging_signature_contribution_summary, by = c("decile", "perm.id"))
AD_average_mut_num_rescale <- AD_average_mut_num_meta %>% group_by(decile, perm.id) %>% summarise( 
							mutation_number = sum(mutation_number),
							sigA_mutation = sum(sigA_mutation),
							sigC_mutation = sum(sigC_mutation),
							sigA_permutation = sum(sigA_permutation),
							sigC_permutation = sum(sigC_permutation),
							permutation_number = sum(permutation_number))

AD_average_mut_num_rescale["sigA_enrichment"] <-
  (AD_average_mut_num_rescale$sigA_mutation * AD_average_mut_num_rescale$mutation_number) / (AD_average_mut_num_rescale$sigA_permutation * AD_average_mut_num_rescale$permutation_number)
AD_average_mut_num_rescale[is.infinite(AD_average_mut_num_rescale$sigA_enrichment), "sigA_enrichment"] <- 1
AD_average_mut_num_rescale["sigC_enrichment"] <-
  (AD_average_mut_num_rescale$sigC_mutation * AD_average_mut_num_rescale$mutation_number) / (AD_average_mut_num_rescale$sigC_permutation * AD_average_mut_num_rescale$permutation_number)
AD_average_mut_num_rescale[is.infinite(AD_average_mut_num_rescale$sigC_enrichment), "sigC_enrichment"] <- 1

AD_average_mut_num_summary_sigA <- AD_average_mut_num_rescale %>% group_by(decile) %>% summarise("average_enrichment" = mean(sigA_enrichment),
                                                                                        "sd_enrichment" = sd(sigA_enrichment))

AD_average_mut_num_summary_sigC <- AD_average_mut_num_rescale %>% group_by(decile) %>% summarise("average_enrichment" = mean(sigC_enrichment),
                                                                                        "sd_enrichment" = sd(sigC_enrichment))
AD_average_mut_num_summary_sigA$group='AD'
AD_average_mut_num_summary_sigC$group='AD'


############################################################
ctrl_aging_signature_contribution_summary <- readRDS(paste0(workdir_ctrl,'AD_snv_permutation_sig.rds'))
ctrl_aging_signature_contribution_summary$decile = as.factor(ctrl_aging_signature_contribution_summary$decile)
ctrl_aging_signature_contribution_summary$perm.id = as.integer(ctrl_aging_signature_contribution_summary$perm.id)
ctrl_mut_num$decile <- as.factor(ctrl_mut_num$decile)
ctrl_mut_num$perm.id <- as.integer(as.character(ctrl_mut_num$perm.id))

ctrl_average_mut_num_meta <- left_join(ctrl_mut_num, ctrl_aging_signature_contribution_summary, by = c("decile", "perm.id"))
ctrl_average_mut_num_rescale <- ctrl_average_mut_num_meta %>% group_by(decile,  perm.id) %>% summarise( 
							mutation_number = sum(mutation_number),
							sigA_mutation = sum(sigA_mutation),
							sigC_mutation = sum(sigC_mutation),
							sigA_permutation = sum(sigA_permutation),
							sigC_permutation = sum(sigC_permutation),
							permutation_number = sum(permutation_number))


ctrl_average_mut_num_rescale["sigA_enrichment"] <-
  (ctrl_average_mut_num_rescale$sigA_mutation * ctrl_average_mut_num_rescale$mutation_number) / (ctrl_average_mut_num_rescale$sigA_permutation * ctrl_average_mut_num_rescale$permutation_number)
ctrl_average_mut_num_rescale[is.infinite(ctrl_average_mut_num_rescale$sigA_enrichment), "sigA_enrichment"] <- 1
ctrl_average_mut_num_rescale["sigC_enrichment"] <-
  (ctrl_average_mut_num_rescale$sigC_mutation * ctrl_average_mut_num_rescale$mutation_number) / (ctrl_average_mut_num_rescale$sigC_permutation * ctrl_average_mut_num_rescale$permutation_number)
ctrl_average_mut_num_rescale[is.infinite(ctrl_average_mut_num_rescale$sigC_enrichment), "sigC_enrichment"] <- 1

ctrl_average_mut_num_summary_sigA <- ctrl_average_mut_num_rescale %>% group_by(decile) %>% summarise("average_enrichment" = mean(sigA_enrichment),
                                                                                        "sd_enrichment" = sd(sigA_enrichment))
ctrl_average_mut_num_summary_sigC <- ctrl_average_mut_num_rescale %>% group_by(decile) %>% summarise("average_enrichment" = mean(sigC_enrichment),
                                                                                        "sd_enrichment" = sd(sigC_enrichment))

ctrl_average_mut_num_summary_sigA$group='ctrl'
ctrl_average_mut_num_summary_sigC$group='ctrl'


average_mut_num_summary_sigA = rbind(ctrl_average_mut_num_summary_sigA, AD_average_mut_num_summary_sigA)
average_mut_num_summary_sigA = rbind(average_mut_num_summary_sigA, tau_average_mut_num_summary_sigA)

average_mut_num_summary_sigC = rbind(ctrl_average_mut_num_summary_sigC, AD_average_mut_num_summary_sigC)
average_mut_num_summary_sigC = rbind(average_mut_num_summary_sigC, tau_average_mut_num_summary_sigC)

average_mut_num_summary_sigA$decile = as.numeric(average_mut_num_summary_sigA$decile)
model = lm(average_enrichment ~ decile + group, data = average_mut_num_summary_sigA)
summary(model)

average_mut_num_summary_sigC$decile = as.numeric(average_mut_num_summary_sigC$decile)
model = lm(average_enrichment ~ decile + group, data = average_mut_num_summary_sigC)
summary(model)

p <- ggplot(average_mut_num_summary_sigA, aes(x = decile, y = average_enrichment, group = group, color = group)) +
  geom_line(size=2) +
  geom_point(size=2, colour="black") +
  scale_color_manual(values = cols) +
  geom_errorbar(aes(ymin = average_enrichment-sd_enrichment, ymax = average_enrichment+sd_enrichment), width = 0.2) +
  theme_classic() +
  ylim(c(0.5, 1.5)) +
  labs(x = "Gene expression levels",
       y = "Mutation enrichment ratio \n (observed/expected)",
       title = "Signature A")
ggsave(paste0(res_dir,"sigA_enrichment_sc.pdf"), plot=p, height = 4, width = 6)

p <- ggplot(average_mut_num_summary_sigC, aes(x = decile, y = average_enrichment, group = group, color = group)) +
  geom_line(size=2) +
  geom_point(size=2, colour="black") +
  scale_color_manual(values = cols) +
  geom_errorbar(aes(ymin = average_enrichment-sd_enrichment, ymax = average_enrichment+sd_enrichment), width = 0.2) +
  theme_classic() +
  ylim(c(-5.5, 7.5)) +
  scale_y_continuous(breaks=seq(-5,7,1))
  labs(x = "Gene expression levels",
       y = "Mutation enrichment ratio \n (observed/expected)",
       title = "Signature C")
ggsave(paste0(res_dir,"sigC_enrichment_sc.pdf"), plot=p, height = 4, width = 6)




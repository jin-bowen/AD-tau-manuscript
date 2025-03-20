library(MutationalPatterns)
library(ggplot2)
library(reshape2)
library(BSgenome)
library(dplyr)
library(tidyr)

source('/home/boj924/AD_Tau_PTA/analysis_scripts/specturm.R')

#### Fit SCAN2 calls to META-CS ds/ssIndel signatures ####
# Get original ds/ss signatures
nmf_res.ds <- readRDS("/home/boj924/Tn5-duplex-calling/results/tlpk_indel_ds_sig_ctrl.rds")
nmf_res.ss <- readRDS("/home/boj924/Tn5-duplex-calling/results/tlpk_indel_ss_sig_ctrl.rds")
metacs_signatures <- as.data.frame(cbind(nmf_res.ds, nmf_res.ss))
names(metacs_signatures) <- c("Signature_dsIndel_filtered", "Signature_ssIndel_filtered")
#normalize by signature (96 fractions for each signature sum up to 1)
metacs_signatures <- t(t(metacs_signatures)/colSums(metacs_signatures))

# NMF on ds/ss signatures
r<-2
nmf_res <- extract_signatures(metacs_signatures, rank = r, single_core = T)
colnames(nmf_res$signatures) <- c("Signature_METACS_S1","Signature_METACS_S2")
rownames(nmf_res$contribution) <- c("Signature_METACS_S1","Signature_METACS_S2")
saveRDS(nmf_res, "/home/boj924/Tn5-duplex-calling/results/nmf_indel_ctrl.S1S2.rds")

S1S2 <- nmf_res$signatures
S1S2_contribution <- nmf_res$contribution
#normalize by signature (96 fractions for each signature sum up to 1)
S1S2 <- t(t(S1S2)/colSums(S1S2))
# normalize S1/S2 contribution to ds/ss
S1S2_contribution_norm <- S1S2_contribution/rowSums(S1S2_contribution)
S1S2_contribution_norm <- as.data.frame((S1S2_contribution_norm))
S1S2_contribution_norm$Signature <- rownames(S1S2_contribution_norm)

# Read input matrix and fit to S1S2
mut_mat <- readRDS("/home/boj924/AD_Tau_PTA/results/tauADCtrl_indelcounts.rds")
fit_res <- fit_to_signatures(mut_mat, S1S2)
df <- as.data.frame(fit_res$contribution)
df$Signature <- rownames(df)
df <- melt(df, id.vars="Signature", variable.name="Cell_ID", value.name="Contribution")

# Reconstruct ds/ss contribution to PTA using S1/S2 contribution to ds/ss
df <- merge(df, S1S2_contribution_norm, by="Signature")
df <- df %>% group_by(Cell_ID, Signature) %>%
  mutate(Contribution_dsIndel = Contribution * Signature_dsIndel_filtered,
         Contribution_ssIndel = Contribution * Signature_ssIndel_filtered)
df1 <- df %>% group_by(Cell_ID) %>%
  summarize(Contribution_dsIndel_sum = sum(Contribution_dsIndel),
            Contribution_ssIndel_sum = sum(Contribution_ssIndel))
df1$Cell_ID <- as.character(df1$Cell_ID)

mut_mat_re <- c()
for (cell in colnames(mut_mat)){
  reconstruct <- metacs_signatures[,1]*df1$Contribution_dsIndel_sum[df1$Cell_ID==cell] + metacs_signatures[,2]*df1$Contribution_ssIndel_sum[df1$Cell_ID==cell]
  mut_mat_re <- cbind(mut_mat_re, reconstruct)
}
colnames(mut_mat_re) <- colnames(mut_mat)
names(df1)[2:3] <- c("Signature_dsIndel_re", "Signature_ssIndel_re")
df1_sum = rowSums(df1[,c("Signature_dsIndel_re", "Signature_ssIndel_re")])
df1_frac = df1[,c("Signature_dsIndel_re", "Signature_ssIndel_re")] / df1_sum
df1_frac$Cell_ID = df1$Cell_ID

meta <- read.csv('/home/boj924/AD_Tau_PTA/metafiles/all_meta.tab', header=TRUE)
meta_df = merge(df1_frac, meta, by.x="Cell_ID", by.y="sample")

clinic = read.csv('/home/boj924/AD_Tau_PTA/metafiles/all_clinic', sep="\t", header=TRUE)
meta_clinic_df = merge(meta_df, clinic, by.x="donor", by.y="Case_ID")
meta_clinic_df = meta_clinic_df[order(meta_clinic_df$group, meta_clinic_df$Age), ]

burden =as.data.frame(read.csv("/home/boj924/AD_Tau_PTA/results/tau_AD_ctrl_res.tab", header=TRUE))
snv_burden = burden[burden['mutation_type']=='indel', ]
meta_clinic_burden = merge(meta_clinic_df, snv_burden, by="Cell_ID")
meta_clinic_burden$Signature_dsIndel_re_burden = meta_clinic_burden$Signature_dsIndel_re * meta_clinic_burden$burden
meta_clinic_burden$Signature_ssIndel_re_burden = meta_clinic_burden$Signature_ssIndel_re * meta_clinic_burden$burden
saveRDS(meta_clinic_burden, "/home/boj924/Tn5-duplex-calling/results/PTA_fit_indel-ctrl-ss.rds")

for (grp in c('ctrl')){
meta_clinic_sub = meta_clinic_burden[meta_clinic_burden$group==grp,]
meta_clinic_sub = meta_clinic_sub[meta_clinic_sub['Cell_ID']!='1995P_201001E3',]
meta_clinic_sub = meta_clinic_sub %>% drop_na()

meta_clinic_sub_melt <- melt(meta_clinic_sub[c('Cell_ID','Signature_dsIndel_re_burden','Signature_ssIndel_re_burden')], 
			id.vars=c("Cell_ID"), measure.vars = c('Signature_dsIndel_re_burden','Signature_ssIndel_re_burden'))
meta_clinic_sub_melt$Cell_ID <- factor(meta_clinic_sub_melt$Cell_ID, levels = meta_clinic_sub_melt$Cell_ID[order(meta_clinic_sub$Age)])


meta_clinic_sub_percent <- melt(meta_clinic_sub[c('Cell_ID','Signature_dsIndel_re','Signature_ssIndel_re')], 
			id.vars=c("Cell_ID"), measure.vars = c('Signature_dsIndel_re','Signature_ssIndel_re'))
meta_clinic_sub_percent$Cell_ID <- factor(meta_clinic_sub_percent$Cell_ID, levels = meta_clinic_sub_percent$Cell_ID[order(meta_clinic_sub$Age)])


p=ggplot(meta_clinic_sub_melt, aes(fill=variable, y=value, x=Cell_ID)) + 
    geom_bar(position="stack", stat="identity") +
    ylim(0, 8000) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(paste0("/home/boj924/AD_Tau_PTA/results/",grp,"_ds_ss_indel-ctrl",".pdf"), plot = p, width = 8, height = 4, dpi = 300)

p=ggplot(meta_clinic_sub_percent, aes(fill=variable, y=value, x=Cell_ID)) + 
    geom_bar(position="stack", stat="identity") +
    ylim(0, 1) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(paste0("/home/boj924/AD_Tau_PTA/results/",grp,"_ds_ss_indel-ctrl_percent",".pdf"), plot = p, width = 8, height = 4, dpi = 300)
}

p=plot_indel_profile(S1S2)
ggsave(paste0("/home/boj924/AD_Tau_PTA/results/","ds_ss_indel-ctrl_spectrum",".pdf"), plot = p, width = 8, height = 6, dpi = 300)


library(MutationalPatterns)
library(ggplot2)
library(reshape2)
library(BSgenome)
library(dplyr)
library(tidyr)
library(forcats)

source('/home/boj924/AD_Tau_PTA/analysis_scripts/specturm.R')
source('/home/boj924/Tn5-duplex-calling/qc_script/plot_spectrum_mod.R')

#### Fit SCAN2 calls to META-CS ds/ssSNV signatures ####
# Get original ds/ss signatures
nmf_res.ds <- readRDS("/home/boj924/Tn5-duplex-calling/results/tlpk_AD_ds_sig.rds")
nmf_res.ss <- readRDS("/home/boj924/Tn5-duplex-calling/results/tlpk_AD_ss_sig.rds")
metacs_signatures <- as.data.frame(cbind(nmf_res.ds, nmf_res.ss))
names(metacs_signatures) <- c("Signature_dsSNV_filtered", "Signature_ssSNV_filtered")
#normalize by signature (96 fractions for each signature sum up to 1)
metacs_signatures <- t(t(metacs_signatures)/colSums(metacs_signatures))

# NMF on ds/ss signatures
r<-2
nmf_res <- extract_signatures(metacs_signatures, rank = r, single_core = T)
colnames(nmf_res$signatures) <- c("Signature_METACS_S1","Signature_METACS_S2")
rownames(nmf_res$contribution) <- c("Signature_METACS_S1","Signature_METACS_S2")
saveRDS(nmf_res, "/home/boj924/Tn5-duplex-calling/results/nmf_res.S1S2.rds")

nmf_res=readRDS("/home/boj924/Tn5-duplex-calling/results/nmf_res.S1S2.rds")
S1S2 <- nmf_res$signatures
S1S2_contribution <- nmf_res$contribution
#normalize by signature (96 fractions for each signature sum up to 1)
S1S2 <- t(t(S1S2)/colSums(S1S2))
# normalize S1/S2 contribution to ds/ss
S1S2_contribution_norm <- S1S2_contribution/rowSums(S1S2_contribution)
S1S2_contribution_norm <- as.data.frame((S1S2_contribution_norm))
S1S2_contribution_norm$Signature <- rownames(S1S2_contribution_norm)





# Read input matrix and fit to S1S2
mut_mat <- readRDS("/home/boj924/AD_Tau_PTA/results/tauADctrl_mut.rds")
mut_mat = mut_mat[, !(colnames(mut_mat) %in% c('1995P_201001E3'))]
fit_res <- fit_to_signatures(mut_mat, S1S2)
df <- as.data.frame(fit_res$contribution)
df$Signature <- rownames(df)
df <- melt(df, id.vars="Signature", variable.name="Cell_ID", value.name="Contribution")

# Reconstruct ds/ss contribution to PTA using S1/S2 contribution to ds/ss
df <- merge(df, S1S2_contribution_norm, by="Signature")
df <- df %>% group_by(Cell_ID, Signature) %>%
  mutate(Contribution_dsSNV = Contribution * Signature_dsSNV_filtered,
         Contribution_ssSNV = Contribution * Signature_ssSNV_filtered)
df1 <- df %>% group_by(Cell_ID) %>%
  summarize(Contribution_dsSNV_sum = sum(Contribution_dsSNV),
            Contribution_ssSNV_sum = sum(Contribution_ssSNV))
df1$Cell_ID <- as.character(df1$Cell_ID)

mut_mat_re <- c()
for (cell in colnames(mut_mat)){
  reconstruct <- metacs_signatures[,1]*df1$Contribution_dsSNV_sum[df1$Cell_ID==cell] + metacs_signatures[,2]*df1$Contribution_ssSNV_sum[df1$Cell_ID==cell]
  mut_mat_re <- cbind(mut_mat_re, reconstruct)
}
colnames(mut_mat_re) <- colnames(mut_mat)
names(df1)[2:3] <- c("Signature_dsSNV_re", "Signature_ssSNV_re")
df1_sum = rowSums(df1[,c("Signature_dsSNV_re", "Signature_ssSNV_re")])
df1_frac = df1[,c("Signature_dsSNV_re", "Signature_ssSNV_re")] / df1_sum
df1_frac$Cell_ID = df1$Cell_ID

meta <- read.csv('/home/boj924/AD_Tau_PTA/metafiles/all_meta.tab', header=TRUE)
meta_df = merge(df1_frac, meta, by.x="Cell_ID", by.y="sample")

clinic = read.csv('/home/boj924/AD_Tau_PTA/metafiles/all_clinic', sep="\t", header=TRUE)
meta_clinic_df = merge(meta_df, clinic, by.x="donor", by.y="Case_ID")
meta_clinic_df = meta_clinic_df[order(meta_clinic_df$group, meta_clinic_df$Age), ]

burden =as.data.frame(read.csv("/home/boj924/AD_Tau_PTA/results/tau_AD_ctrl_res.tab", header=TRUE))
snv_burden = burden[burden['mutation_type']=='snv', ]
meta_clinic_burden = merge(meta_clinic_df, snv_burden, by="Cell_ID")
meta_clinic_burden$Signature_dsSNV_re_burden = meta_clinic_burden$Signature_dsSNV_re * meta_clinic_burden$burden
meta_clinic_burden$Signature_ssSNV_re_burden = meta_clinic_burden$Signature_ssSNV_re * meta_clinic_burden$burden
saveRDS(meta_clinic_burden, "/home/boj924/Tn5-duplex-calling/results/PTA_fit-ss.rds")

meta_clinic_burden$class = if_else(meta_clinic_burden$group=='ctrl', 'ctrl', 'AD')
meta_clinic_burden$class = factor(meta_clinic_burden$class, level=c('ctrl','AD'))
meta_clinic_burden$group = factor(meta_clinic_burden$group, level=c('ctrl', 'AD','Tau','noTau'))

meta_clinic_burden = meta_clinic_burden[meta_clinic_burden['Cell_ID']!='1995P_201001E3',]
meta_clinic_burden = meta_clinic_burden %>% drop_na()

meta_clinic_melt <- melt(meta_clinic_burden[c('Cell_ID','class','donor','group','Age', 'Signature_dsSNV_re_burden','Signature_ssSNV_re_burden')],
                        id.vars=c("Cell_ID", "class","donor","group","Age"), measure.vars = c('Signature_dsSNV_re_burden','Signature_ssSNV_re_burden'))
meta_clinic_melt_sort =  meta_clinic_melt[order(meta_clinic_melt$class, meta_clinic_melt$Age, meta_clinic_melt$group), ]
meta_clinic_melt_sort$Cell_ID = fct_inorder(meta_clinic_melt_sort$Cell_ID)

meta_clinic_percent_melt <- melt(meta_clinic_burden[c('Cell_ID','class','donor','group','Age','Signature_dsSNV_re','Signature_ssSNV_re')],
                        id.vars=c("Cell_ID", "class","donor","group","Age"), measure.vars = c('Signature_dsSNV_re','Signature_ssSNV_re'))
meta_clinic_percent_melt_sort =  meta_clinic_percent_melt[order(meta_clinic_percent_melt$class, meta_clinic_percent_melt$Age, meta_clinic_percent_melt$group), ]
meta_clinic_percent_melt_sort$Cell_ID = fct_inorder(meta_clinic_percent_melt_sort$Cell_ID)

write.csv(meta_clinic_melt,'/home/boj924/AD_Tau_PTA/results/All_ss_abs.csv', quote=F)
write.csv(meta_clinic_percent_melt,'/home/boj924/AD_Tau_PTA/results/All_ss_percent.csv', quote=F)

p=ggplot(meta_clinic_melt_sort, aes(fill=variable, y=value, x=Cell_ID)) +
    geom_bar(position="stack", stat="identity") +
    ylim(0, 3000) +
    scale_fill_manual(values=c("#9999CC", "#66CC99"))+
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(paste0("/home/boj924/AD_Tau_PTA/results/all_ds_ss_snv-AD",".pdf"), plot = p, width = 16, height = 4, dpi = 300)

p=ggplot(meta_clinic_percent_melt_sort, aes(fill=variable, y=value, x=Cell_ID)) +
    geom_bar(position="stack", stat="identity") + 
    ylim(0, 1.1) +
    scale_fill_manual(values=c("#9999CC", "#66CC99"))+
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(paste0("/home/boj924/AD_Tau_PTA/results/all_ds_ss_snv-AD_percent",".pdf"), plot = p, width = 16, height = 4, dpi = 300)

p=plot_snv_profile(S1S2)
ggsave(paste0("/home/boj924/AD_Tau_PTA/results/","ds_ss_snv-AD_spectrum",".pdf"), plot = p, width = 8, height = 6, dpi = 300)














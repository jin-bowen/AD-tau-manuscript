library(MutationalPatterns)
library(ggplot2)
library(reshape2)
library(BSgenome)
library(dplyr)
library(tidyr)
library(forcats)

source('/home/boj924/AD_Tau_PTA/analysis_scripts/specturm.R')

#### Fit SCAN2 calls to META-CS ds/ssIndel signatures ####
# Get original ds/ss signatures
nmf_res.ds <- readRDS("/home/boj924/Tn5-duplex-calling/results/tlpk_indel_ds_sig_AD.rds")
nmf_res.ss <- readRDS("/home/boj924/Tn5-duplex-calling/results/tlpk_indel_ss_sig_AD.rds")
metacs_signatures <- as.data.frame(cbind(nmf_res.ds, nmf_res.ss))
names(metacs_signatures) <- c("Signature_dsIndel_filtered", "Signature_ssIndel_filtered")
#normalize by signature (96 fractions for each signature sum up to 1)
metacs_signatures <- t(t(metacs_signatures)/colSums(metacs_signatures))

# NMF on ds/ss signatures
r<-2
nmf_res <- extract_signatures(metacs_signatures, rank = r, single_core = T)
colnames(nmf_res$signatures) <- c("Signature_METACS_S1","Signature_METACS_S2")
rownames(nmf_res$contribution) <- c("Signature_METACS_S1","Signature_METACS_S2")
saveRDS(nmf_res, "/home/boj924/Tn5-duplex-calling/results/nmf_indel_AD.S1S2.rds")


nmf_res=readRDS("/home/boj924/Tn5-duplex-calling/results/nmf_indel_AD.S1S2.rds")
S1S2 <- nmf_res$signatures
S1S2_contribution <- nmf_res$contribution
#normalize by signature (96 fractions for each signature sum up to 1)
S1S2 <- t(t(S1S2)/colSums(S1S2))
# normalize S1/S2 contribution to ds/ss
S1S2_contribution_norm <- S1S2_contribution/rowSums(S1S2_contribution)
S1S2_contribution_norm <- as.data.frame((S1S2_contribution_norm))
S1S2_contribution_norm$Signature <- rownames(S1S2_contribution_norm)

signature_all <- readRDS("/home/boj924/AD_Tau_PTA/results/indel_musical_all.rds")
signature_fit <- fit_to_signatures(S1S2, as.matrix(signature_all))
signature_contribution <- apply(signature_fit$contribution, 2, function(x){x/sum(x)})

df <- melt(signature_contribution)
colnames(df) <- c("sig1", "sig2", "Value")

p=ggplot(df, aes(x = sig1, y = sig2, fill = Value)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "red") +
  theme_minimal() +
  labs(title = "Heatmap from Matrix",
       x = "X-axis",
       y = "Y-axis",
       fill = "Value")
ggsave(paste0("/home/boj924/AD_Tau_PTA/results/ss_indel_AD_musical",".pdf"), plot = p, width = 12, height = 4, dpi = 300)



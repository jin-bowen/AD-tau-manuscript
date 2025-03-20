args=commandArgs(trailingOnly=TRUE)

library(MutationalPatterns)
library(ggplot2)
library(reshape2)
library(BSgenome)


plot_profile <- function(mut_matrix, colors = NA) {
  library(dplyr)
  if (colors=='indel'){
  colors <- c(
    "#FDBE6F", "#FF8001", "#B0DD8B", "#36A12E", "#FDCAB5", "#FC8A6A",
    "#F14432", "#BC141A", "#D0E1F2", "#94C4DF", "#4A98C9", "#1764AB",
    "#E2E2EF", "#B6B6D8", "#8683BD", "#61409B"
  )} else {
  colors <- c(
  "#2EBAED", "#000000", "#DE1C14",
  "#D4D2D2", "#ADCC54", "#F0D0CE"
  )}

  # Get substitution and context from rownames and make long.
  tb <- mut_matrix %>%
    as.data.frame() %>%
    tibble::rownames_to_column("full_context") %>%
    dplyr::mutate(
      substitution = stringr::str_replace(full_context, "\\w\\[(.*)\\]\\w", "\\1"),
      context = stringr::str_replace(full_context, "\\[.*\\]", "\\.")
    ) %>%
    dplyr::select(-full_context) %>%
    tidyr::pivot_longer(c(-substitution, -context), names_to = "sample", values_to = "count") %>% 
    dplyr::mutate(sample = factor(sample, levels = unique(sample)))

    width <- 1
    spacing <- 0

  # Create figure
  plot <- ggplot(data = tb, aes(
    x = context,
    y = count,
    fill = substitution,
    width = width
  )) +
    geom_bar(stat = "identity", colour = "black", size = .2) +
    scale_fill_manual(values = colors) +
    facet_grid(sample ~ substitution) +
    ylab("Absolute difference(Tau+ - Tau-)") +
    guides(fill = FALSE) +
    theme_bw() +
    theme(
      axis.title.y = element_text(size = 12, vjust = 1),
      axis.text.y = element_text(size = 8),
      axis.title.x = element_text(size = 12),
      axis.text.x = element_text(size = 5, angle = 90, vjust = 0.5),
      strip.text.x = element_text(size = 9),
      strip.text.y = element_text(size = 9),
      strip.background = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.spacing.x = unit(spacing, "lines")
    )

  return(plot)
}

workdir = '/home/boj924/AD_Tau_PTA/results/'
dir1 = '/home/boj924/AD_Tau_PTA/results/control_neurons/snv/'
dir2 = '/home/boj924/AD_Tau_PTA/results/tau_neurons/snv/'
dir3 = '/home/boj924/AD_Tau_PTA/results/AD_neurons/snv/'

####compare different subset
vcf_files <- list.files(c(dir1,dir2,dir3),pattern = "*_snv.vcf", full.names=TRUE)
outname='tauADCtrl'

ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
library(ref_genome, character.only = TRUE)

getSampleName = function(x){ gsub('_snv.vcf','', basename(x))  }
sample_names = lapply(vcf_files,getSampleName )
sample_names = unlist(sample_names)

grl <- read_vcfs_as_granges(vcf_files, sample_names, ref_genome)
snv_grl <- get_mut_type(grl, type = "snv")
muts <- mutations_from_vcf(grl[[1]])
context <- mut_context(grl[[1]], ref_genome)
type_context <- type_context(grl[[1]], ref_genome)
type_occurrences <- mut_type_occurrences(grl, ref_genome)

mut_mat <- mut_matrix(vcf_list = grl, ref_genome = ref_genome)
mut_mat <- mut_mat + 0.0001
saveRDS(mut_mat, file = "tauADctrl_mut.rds")


mut_mat = readRDS('results/tauADctrl_mut.rds')
# combined sig
aging_signatures=as.matrix(read.csv("/n/data1/bwh/pathology/miller/lab/ref/aging.mutation_pattern.3sigs.tsv"))
aging_signatures_renamed=as.matrix(cbind(aging_signatures[,2],aging_signatures[,3],aging_signatures[,1]))
colnames(aging_signatures_renamed)=c("Signature.A","Signature.B","Signature.C")

artifact_signatures=read.csv("/n/data1/bwh/pathology/miller/lab/ref/Signature_scEF_Mia.txt",header=TRUE)
artifact_signatures=as.matrix(artifact_signatures[,-c(1:2)])
aging_signatures_all=cbind(aging_signatures_renamed,artifact_signatures)
aging_signatures_all=t(t(aging_signatures_all)/colSums(aging_signatures_all))
aging_signatures_all=aging_signatures_all[,-2]

if(sum(rownames(aging_signatures_all) != rownames(mut_mat)) >1){
        print("column name does not match")
}

fit_res_aging=fit_to_signatures(mut_mat,aging_signatures_all)

aging_signatures_contribute = fit_res_aging$contribution
aging_signatures_contribute = t(t(aging_signatures_contribute)/colSums(aging_signatures_contribute))
write.csv(t(aging_signatures_contribute), paste0(workdir,"/snv_AgingSigContribute.csv"),quote=F)


# combined sig
sig_file='/n/data1/bch/genetics/lee/lan/CTE/mut_sig/aging.mutation_pattern.3sigs.tsv'
aging_signatures=as.matrix(read.delim(sig_file))
chem_signatures <- read.table("/home/boj924/AD_Tau_PTA/custome_sig/Mutagen53_sub_signature.tsv",sep="\t",header=T,check.names=FALSE)
sigA = aging_signatures[,2] /sum(aging_signatures[,2])
all_sig = as.matrix(cbind(chem_signatures[2:54], 'SigA'=sigA))

strict_refit <- fit_to_signatures_strict(
  as.matrix(mut_mat),
  as.matrix(all_sig),
  max_delta = 0.004,
  method = "backward")

fit_res_strict <- strict_refit$fit_res
meta <- read.csv('/home/boj924/AD_Tau_PTA/metafiles/all_meta.tab',header=TRUE)
getGrpName = function(x){ meta[meta$sample==x,'group']  }
grpname = lapply(  colnames(fit_res_strict$contribution), getGrpName )
grpname =unlist(grpname)
mat = fit_res_strict$contribution
colnames(mat) = grpname


ctrl_res = rowSums( mat[,colnames(mat) %in% c('ctrl')] )
ctrl_res = ctrl_res/sum(ctrl_res)
Tau_res = rowSums( mat[,colnames(mat) %in% c('Tau')] )
Tau_res = Tau_res/sum(Tau_res)
noTau_res = rowSums( mat[,colnames(mat) %in% c('noTau')] )
noTau_res = noTau_res/sum(noTau_res)
AD_res = rowSums( mat[,colnames(mat) %in% c('AD')] )
AD_res = AD_res/sum(AD_res)
combine_res = t(rbind(ctrl_res, Tau_res, noTau_res,AD_res))
nonzero = combine_res[((combine_res[,"Tau_res"]>0.01) | 
			(combine_res[,"noTau_res"]>0.01) |
			(combine_res[,"AD_res"]>0.01) |
			(combine_res[,"ctrl_res"]>0.01)),]
strict_refit_sub <- fit_to_signatures_strict(
  mut_mat,
  all_sig[,row.names(nonzero)],
  max_delta = 0.004,
  method = "best_subset")
saveRDS(strict_refit_sub, file = "all_combineSig.rds")

cosmics=get_known_signatures()
strict_refit2 <- fit_to_signatures_strict(
  as.matrix(mut_mat),
  as.matrix(cosmics),
  max_delta = 0.004,
  method = "backward")

fit_res_strict2 <- strict_refit2$fit_res
meta <- read.csv('/home/boj924/AD_Tau_PTA/metafiles/all_meta.tab',header=TRUE)
getGrpName = function(x){ meta[meta$sample==x,'group']  }
grpname = lapply(  colnames(fit_res_strict2$contribution), getGrpName )
grpname =unlist(grpname)
mat = fit_res_strict2$contribution
colnames(mat) = grpname

ctrl_res = rowSums( mat[,colnames(mat) %in% c('ctrl')] )
ctrl_res = ctrl_res/sum(ctrl_res)
Tau_res = rowSums( mat[,colnames(mat) %in% c('Tau')] )
Tau_res = Tau_res/sum(Tau_res)
noTau_res = rowSums( mat[,colnames(mat) %in% c('noTau')] )
noTau_res = noTau_res/sum(noTau_res)
AD_res = rowSums( mat[,colnames(mat) %in% c('AD')] )
AD_res = AD_res/sum(AD_res)
combine_res = t(rbind(ctrl_res, Tau_res, noTau_res,AD_res))
nonzero = combine_res[((combine_res[,"Tau_res"]>0.01) | 
			(combine_res[,"noTau_res"]>0.01) |
			(combine_res[,"AD_res"]>0.01) |
			(combine_res[,"ctrl_res"]>0.01)),]
strict_refit_sub <- fit_to_signatures_strict(
  mut_mat,
  cosmics[,row.names(nonzero)],
  max_delta = 0.004,
  method = "best_subset")
saveRDS(strict_refit_sub, file = "all_cosmicSig.rds")

##diff
tau_diff <- read.csv("/home/boj924/AD_Tau_PTA/tau_diff.mat",header=T,check.names=FALSE,row.names='muttype')
plot_96_profile(tau_diff)


##de novo signature
#library("NMF")
#estimate <- nmf(mut_mat, rank = 2:10, method = "brunet", 
#                nrun = 200, seed = 123456, .opt = "v-p")
#fig=plot(estimate)
#ggsave("/home/boj924/AD_Tau_PTA/results/all_snv_denovo_estimate.png", plot = fig, width = 16, height = 8, dpi = 300)

#r=2
#nmf_res <- extract_signatures(mut_mat, rank = r, single_core = T)
#colnames(nmf_res$signatures) <- paste0("Signature N",1:r)
#rownames(nmf_res$contribution) <- paste0("Signature N",1:r)
#saveRDS(nmf_res,paste0("/home/boj924/AD_Tau_PTA/results/AD.nmf_res",r,".rds"))
#fit <- as.data.frame(nmf_res$contribution)
#fit$Signature <- rownames(fit)
#fit <- melt(fit, id.vars="Signature", variable.name="Cell_ID", value.name="Contribution")
#write.table(fit, paste0("/home/boj924/AD_Tau_PTA/results/AD.nmf_res",r,"_contribution.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
#
#res=readRDS('/home/boj924/AD_Tau_PTA/results/AD.nmf_res2.rds')
#meta <- read.csv('/home/boj924/AD_Tau_PTA/metafiles/all_meta.tab',header=TRUE)
#cos_sim = cos_sim_matrix(mut_mat, res$reconstructed)
#cosine_heatmap = plot_cosine_heatmap(round(cos_sim,2), plot_values = T, cluster_rows = F)
#
#cosmics=get_known_signatures()
#cos_sim_COSMIC = cos_sim_matrix(cosmics, res$signatures)
#
#strict_refit <- fit_to_signatures_strict(
#  mut_mat,
#  cosmics,
#  max_delta = 0.004,
#  method = "backward")
#fit_res_strict <- strict_refit$fit_res
#contribute = fit_res_strict$contribution
#pct_contribute = t(t(contribute)/colSums(contribute))
#sort_pct_contribute = sort(colMeans(t(pct_contribute)))
#sort_pct_contribute[sort_pct_contribute>0.01]
#
#
#submeta = meta[meta$sample %in% colnames(res$contribution),]
#sample_order = submeta[order(submeta$group),]$sample
#plot_contribution_bar = plot_contribution(res$contribution[,sample_order], coord_flip = T)
#ggsave("/home/boj924/AD_Tau_PTA/results/snv_denovo_contri2.png", plot = plot_contribution_bar, width = 8, height = 16, dpi = 300)
#ggsave("/home/boj924/AD_Tau_PTA/results/AD.nmf_res2.png", plot = cosine_heatmap_COSMIC, width = 8, height = 16, dpi = 300)
#
#
# transcription bias
# BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")
# BiocManager::install("AnnotationDbi")
# BiocManager::install("GenomicFeatures")
library("TxDb.Hsapiens.UCSC.hg19.knownGene")
Gene_hg19 <- GenomicFeatures::genes(TxDb.Hsapiens.UCSC.hg19.knownGene)
strand <- mut_strand(grl[[1]], Genes_hg19)
mut_mat_s  <- mut_matrix_stranded(grl, ref_genome, Gene_hg19)

meta=read.csv('/home/boj924/AD_Tau_PTA/metafiles/all_meta.tab')
getGrpName = function(x){ meta[meta$sample==x,'group']  }
grpname = lapply(  colnames(mut_mat_s), getGrpName )
grpname =unlist(grpname)

strand_counts <- strand_occurrences(mut_mat_s, by = grpname)
strand_bias <- strand_bias_test(strand_counts)
strand_bias_strict <- strand_bias_test(strand_counts,fdr_cutoffs=0.05) 
fig <- plot_strand_bias(strand_bias_strict, mode = "relative")
ggsave("/home/boj924/AD_Tau_PTA/results/transcription_bias_snv.png", plot = fig, width = 16, height = 8, dpi = 300)
#
##de novo signature decomposition
#mut_mat=readRDS("tauADctrl_mut.rds")
#signatures <- read.table("/home/boj924/AD_Tau_PTA/Mutagen53_COSMIC.tsv",sep="\t",header=T,check.names=FALSE)
#subset_sig = signatures[c('SBS1','SBS5','SBS16','SBS30','OTA (0.08 uM) + S9',
#			'SBS8', 'Potassium bromate (260 uM)','Potassium bromate (875 uM)')]
#subset_sig = as.matrix(subset_sig)
#strict_refit2 <- fit_to_signatures_strict(
#  mut_mat,
#  subset_sig,
#  max_delta = 0.004,
#  method = "best_subset")
#
#fit_res_strict2 <- strict_refit2$fit_res
#
#getGrpName = function(x){ meta[meta$sample==x,'group']  }
#grpname = lapply(  colnames(fit_res_strict2$contribution), getGrpName )
#mat = fit_res_strict2$contribution
#colnames(mat) = grpname
#ctrl_res = rowSums( mat[,colnames(mat) %in% c('ctrl')] )
#ctrl_res = ctrl_res/sum(ctrl_res)
#Tau_res = rowSums( mat[,colnames(mat) %in% c('Tau')] )
#Tau_res = Tau_res/sum(Tau_res)
#noTau_res = rowSums( mat[,colnames(mat) %in% c('noTau')] )
#noTau_res = noTau_res/sum(noTau_res)
#AD_res = rowSums( mat[,colnames(mat) %in% c('AD')] )
#AD_res = AD_res/sum(AD_res)
#combine_res = t(rbind(ctrl_res, Tau_res, noTau_res,AD_res))
#p = plot_contribution(combine_res,
#  coord_flip = TRUE,
#  mode = "relative")
#
#ggsave("/home/boj924/AD_Tau_PTA/results/test1.png", plot = p, width = 8, height = 4, dpi = 300)









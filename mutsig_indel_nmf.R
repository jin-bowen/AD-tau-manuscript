args=commandArgs(trailingOnly=TRUE)

library(MutationalPatterns)
library(ggplot2)
library(BSgenome)
library(reshape2)
library(dplyr)

dir =   '/home/boj924/AD_Tau_PTA/results/all_indel'
outname='tauADCtrl'
vcf_files <- list.files(c(dir),pattern = "*_indel.vcf", full.names=TRUE)

ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
library(ref_genome, character.only = TRUE)

getSampleName = function(x){ gsub('_indel.vcf','', basename(x))  }
sample_names = lapply(vcf_files,getSampleName )
sample_names = unlist(sample_names)
grl <- read_vcfs_as_granges(vcf_files, sample_names, ref_genome, type = "all")
indel_grl <- get_mut_type(grl, type = "indel")
indel_grl <- get_indel_context(indel_grl, ref_genome)
mut_mat <- count_indel_contexts(indel_grl)
mut_mat <- mut_mat + 0.0001
library("NMF")
estimate <- nmf(mut_mat, rank = 1:10, method = "brunet", 
                nrun = 10, seed = 123456, .opt = "v-p")
pdf(pointsize=5, file='indel_denovo_nmf.pdf')
plot(estimate)
dev.off()




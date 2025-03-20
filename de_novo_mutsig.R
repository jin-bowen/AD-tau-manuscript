args=commandArgs(trailingOnly=TRUE)

library(MutationalPatterns)
library(BSgenome)

dir = '/home/boj924/AD_Tau_PTA/results/AD_neurons/'
dir1 = '/home/boj924/AD_Tau_PTA/results/control_neurons/'
outname='ADnctrl'

ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
library(ref_genome, character.only = TRUE)

getSampleName = function(x){ gsub('_snv.vcf.gz','', basename(x))  }

vcf_files <- list.files(c(dir,dir1),pattern = "*_snv.vcf.gz", full.names=TRUE)

sample_names = lapply(vcf_files,getSampleName )
sample_names = unlist(sample_names)

grl <- read_vcfs_as_granges(vcf_files, sample_names, ref_genome)
snv_grl <- get_mut_type(grl, type = "snv")
muts <- mutations_from_vcf(grl[[1]])
context <- mut_context(grl[[1]], ref_genome)
type_context <- type_context(grl[[1]], ref_genome)
type_occurrences <- mut_type_occurrences(grl, ref_genome)

##general dist
p3 <- plot_spectrum(type_occurrences,indv_points = TRUE)

#signature analysis
mut_mat <- mut_matrix(vcf_list = grl, ref_genome = ref_genome)
mut_mat <- mut_mat + 0.0001
estimate <- nmf(mut_mat, rank = 2:10, method = "brunet",
                nrun = 200, seed = 123456, .opt = "v-p")
signatures = get_known_signatures()
for (r in 2:6){
    nmf_res <- extract_signatures(mut_mat, rank = r, single_core = T)
    cosmic_signatures <- get_known_signatures()
    fit_res <- fit_to_signatures(nmf_res$signatures, cosmic_signatures)
    write.table(t(fit_res$contribution), file=paste0(outname,'_denovo_rank',r,'_Sigcontribute.tab'),
		quote = FALSE, sep = ",")
    write.table(t(nmf_res$contribution), file=paste0(outname,'_denovo_rank',r,'_NMFcontribute.tab'),
		quote = FALSE, sep = ",")

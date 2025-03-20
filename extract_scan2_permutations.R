options(stringsAsFactors = FALSE)
library(stringr)

extract_perms_by_cell <- function(rda_path){
  load(rda_path)
  if (nrow(permdata$muts)==0) return(c()) #return empty df if no mut calls
  sample <- permdata$muts$sample[1]
  perms <- permdata$perms
  df.perms <- lapply(1:length(perms), function(i){
    x <- perms[[i]]
    x$sample <- sample
    x$perm.id <- i
    return(x)
  })
  df <- do.call("rbind", df.perms)
  return(df)
}


#for (base_dir in c("/n/data1/bwh/pathology/miller/lab/AD_Tau_PTA/results_hg19/Permute_Tau/", 
#base_dir="/n/data1/bwh/pathology/miller/lab/AD_Tau_PTA/results_hg19/Permute_noTau/"
#sample_to_subject_map <- '/home/boj924/AD_Tau_PTA/AD_tau_config/scan2_sample_to_subject.csv'

#base_dir="/n/scratch/users/b/boj924/AD_hg19/permute/"
#sample_to_subject_map <- '/home/boj924/AD_Tau_PTA/AD_case_config/scan2_sample_to_subject.csv'

base_dir="/n/data1/bwh/pathology/miller/lab/ctrl_hg19/permute/"
sample_to_subject_map <- '/home/boj924/AD_Tau_PTA/Control_config/scan2_sample_to_subject.csv'

perms_by_group_snv <- c()
perms_by_group_indel <- c()
print(base_dir)

meta <- read.table(sample_to_subject_map, header=T, sep=",")
case_ID_list <- unique(meta$subject)

for (case_ID in case_ID_list){
  cell_ID_list <- meta$sample[meta$subject==case_ID]
  for (cell_ID in cell_ID_list){
    rda_path_snv <- paste0(base_dir, "perms_by_sample/", cell_ID, "/snv_pass.rda")
    rda_path_indel <- paste0(base_dir, "perms_by_sample/", cell_ID, "/indel_pass.rda")

    if (!file.exists(rda_path_snv)) { next}

    perms_by_cell_snv <- extract_perms_by_cell(rda_path_snv)
    perms_by_cell_indel <- extract_perms_by_cell(rda_path_indel)
    perms_by_group_snv <- rbind(perms_by_group_snv, perms_by_cell_snv)
    perms_by_group_indel <- rbind(perms_by_group_indel, perms_by_cell_indel)
  }
}

write.table(perms_by_group_snv, paste0(base_dir, "/perms_by_cell_snv.combined.txt"), sep="\t", quote=F, row.names=F)
write.table(perms_by_group_indel, paste0(base_dir, "/perms_by_cell_indel.combined.txt"), sep="\t", quote=F, row.names=F)



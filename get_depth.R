args=commandArgs(trailingOnly=TRUE)

library(scan2)

rdat_file = args[1]
id= args[2]
load(rdat_file)

#depth_quantile <- quantile(results@gatk[results@gatk$chr %in% c(1:22) & 
#	(results@gatk$resampled.training.site==T) & 
#	(results@gatk$muttype=='snv'),]$dp, prob=0:4/4)
#
#write.csv(t(as.data.frame(depth_quantile)), paste0(id, '_snv_depth','.txt'), quote=F)


depth_quantile <- quantile(results@gatk[results@gatk$chr %in% c(1:22) & 
	(results@gatk$resampled.training.site==T) & 
	(results@gatk$muttype=='indel'),]$dp, prob=0:4/4)

write.csv(t(as.data.frame(depth_quantile)), paste0(id, '_indel_depth','.txt'), quote=F)



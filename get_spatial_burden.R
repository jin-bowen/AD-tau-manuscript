args=commandArgs(trailingOnly=TRUE)

library(scan2)

rdat_file = args[1]
id= args[2]
load(rdat_file)

s <- results@gatk[pass == TRUE & muttype == 'indel']
write.csv(s, paste0(id, '_indel','.tab'), row.names=FALSE)

#s_all <- results@gatk[ muttype == 'indel' & somatic.candidate == TRUE]
#write.csv(s_all, paste0(id, '_indel_all','.tab'), row.names=FALSE)

s <- results@gatk[pass == TRUE & muttype == 'snv']
write.csv(s, paste0(id, '_snv','.tab'), row.names=FALSE)

#s_all <- results@gatk[ muttype == 'snv' & somatic.candidate == TRUE]
#write.csv(s_all, paste0(id, '_snv_all','.tab'), row.names=FALSE)



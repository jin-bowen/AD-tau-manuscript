# Extract SCAN2 SNV/indel burden and calls

args = commandArgs(trailingOnly=TRUE)
wd <- args[1]
metafile <- args[2]

print(wd)

# setwd("/n/data1/bch/genetics/lee/lan/CTE/scan2")
setwd(wd)
options(stringsAsFactors = FALSE)

library(scan2)
library(stringr)

##REF: /n/data1/bch/genetics/lee/August/single-cell/extractScan2Result_new.R
extract_results_per_cell <- function(Case_ID, Cell_ID){
	if (!file.exists(paste0(Case_ID,"/call_mutations/",Cell_ID,"/scan2_object_update.rda"))) return()
	load(paste0(Case_ID,"/call_mutations/",Cell_ID,"/scan2_object_update.rda"))

	summary <- data.frame(
		mutation.type="snv", 
		num.pass=results@call.mutations$snv.pass,
		burden.genome=results@mutburden$snv$autosomal$rate.per.gb[2]*5.845,
		sim.negbin.burden=results@spatial.sensitivity$burden$snv$sim.negbin.estimates$equal.burden,
		poisson.burden=results@spatial.sensitivity$burden$snv$poisson.estimates$equal.burden,
		negbin.burden=results@spatial.sensitivity$burden$snv$negbin.estimates$equal.burden,
		mean.burden=results@spatial.sensitivity$burden$snv$mean.estimates$equal.burden,
		results@mutburden$snv$autosomal[2,])

	summary <- rbind(summary, data.frame(
		mutation.type="indel", 
		num.pass=results@call.mutations$indel.pass, 
		burden.genome=results@mutburden$indel$autosomal$rate.per.gb[2]*5.845,
		sim.negbin.burden=results@spatial.sensitivity$burden$indel$sim.negbin.estimates$equal.burden,
		poisson.burden=results@spatial.sensitivity$burden$indel$poisson.estimates$equal.burden,
		negbin.burden=results@spatial.sensitivity$burden$indel$negbin.estimates$equal.burden,
		mean.burden=results@spatial.sensitivity$burden$indel$mean.estimates$equal.burden,
		results@mutburden$indel$autosomal[2,]))
	write.table(summary, paste0(Case_ID,"/call_mutations/",Cell_ID,"/summary.txt"), quote=F,row.names=F,sep="\t")

        snvs <- results@gatk[pass == TRUE & muttype == 'snv']
        indels <- results@gatk[pass == TRUE & muttype == 'indel']
        if(nrow(snvs) > 0) {
                cell.idx <- which(colnames(snvs)==sprintf("X%s",Cell_ID) | colnames(snvs)==Cell_ID)
                df <- data.frame(snvs$chr,snvs$pos,snvs$refnt,snvs$altnt,snvs$scref,snvs$scalt,snvs[,..cell.idx],snvs$bref,snvs$balt,snvs$bulk.gt,snvs$dbsnp)
                colnames(df)=c("Chromosome","Position","REF_Base","ALT_Base","REF_Num","ALT_Num","Genotype","REF_Num_Bulk","ALT_Num_Bulk","Genotype_Bulk","dbSNP_ID")
        } else {
                df <- data.frame(matrix(ncol = 11, nrow = 0))
                colnames(df)=c("Chromosome","Position","REF_Base","ALT_Base","REF_Num","ALT_Num","Genotype","REF_Num_Bulk","ALT_Num_Bulk","Genotype_Bulk","dbSNP_ID")
        }
        write.table(df, paste0(Case_ID,"/call_mutations/",Cell_ID,"/snv_list.txt"), quote=F,row.names=F,sep="\t")

        if(nrow(indels) > 0) {
                cell.idx <- which(colnames(indels)==sprintf("X%s",Cell_ID) | colnames(indels)==Cell_ID)
                df <- data.frame(indels$chr,indels$pos,indels$refnt,indels$altnt,indels$scref,indels$scalt,indels[,..cell.idx],indels$bref,indels$balt,indels$bulk.gt,indels$dbsnp)
                colnames(df)=c("Chromosome","Position","REF_Base","ALT_Base","REF_Num","ALT_Num","Genotype","REF_Num_Bulk","ALT_Num_Bulk","Genotype_Bulk","dbSNP_ID")
        } else {
                df <- data.frame(matrix(ncol = 11, nrow = 0))
                colnames(df)=c("Chromosome","Position","REF_Base","ALT_Base","REF_Num","ALT_Num","Genotype","REF_Num_Bulk","ALT_Num_Bulk","Genotype_Bulk","dbSNP_ID")
        }
        write.table(df, paste0(Case_ID,"/call_mutations/",Cell_ID,"/indel_list.txt"), quote=F,row.names=F,sep="\t")

}

# extract results for each cell
print("extract results for each cell...")
# meta <- read.table("/n/data1/bch/genetics/lee/lan/CTE/data/scan2.panel_meta_with_bam.csv", header=T, sep=",")
meta <- read.table(metafile, header=T, sep=",",stringsAsFactors=F, colClasses=c('character','character'))
for (i in 1:nrow(meta)){
	Case_ID=meta[i, 'Case_ID']
	Cell_ID=meta[i, 'Cell_ID']
	extract_results_per_cell(Case_ID, Cell_ID)
}

# combine summary results
print("combine summary results...")
scan2 <- c()
for (i in 1:nrow(meta)){
	case=meta[i, 'Case_ID']
	cell=meta[i, 'Cell_ID']
        if (!file.exists(sprintf("%s/call_mutations/%s/summary.txt", case, cell))) next
        summary <- read.table(sprintf("%s/call_mutations/%s/summary.txt", case, cell), header=T, sep="\t")
        summary$Case_ID <- case
        summary$Cell_ID <- cell
        scan2 <- rbind(scan2, summary)   
}

# reorder columns
scan2 <- as.data.frame(scan2)
write.table(scan2,"summary.combined.txt",sep="\t",quote=F,row.names=F)


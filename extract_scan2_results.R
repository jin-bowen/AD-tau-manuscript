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
		burden.genome=results@mutburden$snv$rate.per.gb[2]*5.845,
		results@mutburden$snv[2,])
	summary <- rbind(summary, data.frame(
		mutation.type="indel", 
		num.pass=results@call.mutations$indel.pass, 
		burden.genome=results@mutburden$indel$rate.per.gb[2]*5.845,
		results@mutburden$indel[2,])
	)
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
scan2 <- scan2[,c(12,13,1:11)]
write.table(scan2,"summary.combined.txt",sep="\t",quote=F,row.names=F)

# combine mutation calls
print("combine mutation calls...")

all_snv <- c(); all_indel <- c()
for (i in 1:nrow(meta)){
    case=meta[i, 'Case_ID']
    cell=meta[i, 'Cell_ID']
    if (!file.exists(sprintf("%s/call_mutations/%s/snv_list.txt", case, cell))) next
    
    snv <- read.table(sprintf("%s/call_mutations/%s/snv_list.txt", case, cell), header=T, sep="\t")
    if (nrow(snv) > 0) {
	    snv$Case_ID <- case
	    snv$Cell_ID <- cell
	    all_snv <- rbind(all_snv, snv)    	
    }
    indel <- read.table(sprintf("%s/call_mutations/%s/indel_list.txt", case, cell), header=T, sep="\t")
    if (nrow(indel) > 0) {
	    indel$Case_ID <- case
	    indel$Cell_ID <- cell
 		  all_indel <- rbind(all_indel, indel)
    }
}

# reorder columns
all_snv <- as.data.frame(all_snv)
all_indel <- as.data.frame(all_indel)
all_snv <- all_snv[,c(12,13,1:11)]
all_indel <- all_indel[,c(12,13,1:11)]
# convert "TRUE" back to "T" (error in read.table)
all_snv$REF_Base[all_snv$REF_Base=="TRUE"] <- "T"
all_snv$ALT_Base[all_snv$ALT_Base=="TRUE"] <- "T"
all_indel$REF_Base[all_indel$REF_Base=="TRUE"] <- "T"
all_indel$ALT_Base[all_indel$ALT_Base=="TRUE"] <- "T"
# write combined lists
write.table(all_snv,"snv_list.combined.txt",sep="\t",quote=F,row.names=F)
write.table(all_indel,"indel_list.combined.txt",sep="\t",quote=F,row.names=F)


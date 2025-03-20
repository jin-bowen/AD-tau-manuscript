args=commandArgs(trailingOnly=TRUE)

library(scan2)

rdat_file = args[1]
id= args[2]
get.gbp.by.genome <- function(object) {
    if (object@genome.string == 'hs37d5') {
        # 93 contigs includes unplaced; see http://genomewiki.ucsc.edu/index.php/Hg19_Genome_size_statistics.
        total <- 3137161264
        chrx <- 155270560
        chry <- 59373566
        chrm <- 16571
        return((total - chrx - chry - chrm)*2 / 1e9) # = 5.845001134
    } else if (object@genome.string == 'hg38') {
        # 455 contigs; see http://genomewiki.ucsc.edu/index.php/Hg38_100-way_Genome_size_statistics
        total <- 3209286105
        chrx <- 156040895
        chry <- 57227415
        chrm <- 16569
        return((total - chrx - chry - chrm)*2 / 1e9) # = 5.992002452
    } else if (object@genome.string == 'mm10') {
        # 66 contigs; see http://genomewiki.ucsc.edu/index.php/Hg38_100-way_Genome_size_statistics
        total <- 2730871774
        chrx <- 171031299
        chry <- 91744698
        chrm <- 16299
        return((total - chrx - chry - chrm)*2 / 1e9) # = 4.936158956
    } else if (object@genome.string == 'CHM13v2.0') {
        total <- 3117292070
        chrx <- 154259566
        chry <- 62460029
        chrm <- 16569
        return((total - chrx - chry - chrm)*2 / 1e9)
    } else {
        warn(paste('gbp not yet implemented for genome', object@genome.string))
        warn("the mutation burden for this analysis is a placeholder!")
        warn("DO NOT USE!")
        # hopefully returning a negative number will alert people that something
        # has gone wrong so they don't ignore the warning messages above
        return(-1)
    }
}

load(rdat_file)


muttypes <- c('snv', 'indel')
mutburden <- setNames(lapply(muttypes, function(mt) {

  target.fdr = 0.01
  pre.geno.burden <- results@fdr.prior.data[[mt]]$burden[2]
  indel <- results@gatk[muttype == 'indel']
  indel[, condition.pass := static.filter == TRUE & (muttype == 'indel') & lysis.fdr <= target.fdr & mda.fdr <= target.fdr]
  
  if (mt == 'indel') {
  s <- indel[condition.pass == TRUE & muttype == mt]
  write.csv(s, paste0(id, '_',mt,'.tab'), row.names=FALSE)  
#  write.csv(indel, paste0(id, '_',mt,'_unfilter.tab'), row.names=FALSE)   
  } else { 
  s <- results@gatk[pass == TRUE & muttype == mt]
#  write.csv(s, paste0(id, '_',mt,'.tab'), row.names=FALSE)
  }
  
  
  g <- results@gatk[resampled.training.site == TRUE & muttype == mt]

  g[training.site == TRUE,training.pass := 
            (cigar.id.test & cigar.hs.test & dp.test & abc.test & min.sc.alt.test) &
            lysis.fdr <= target.fdr & mda.fdr <= target.fdr]


 
  dptab <- results@depth.profile$dptab
  dptab <- dptab[1:min(max(g$dp)+1, nrow(dptab)),]
  
  q=4
  qstouse <- c(1,2,4)
  qbreaks <- quantile(g$dp, prob=0:q/q)
  s$dpq <- cut(s$dp, qbreaks, include.lowest=T, labels=F)
  s$dpq[s$dpq==3] <- 2 # merge 25-75% into a single bin
  g$dpq <- cut(g$dp, qbreaks, include.lowest=T, labels=F)
  g$dpq[g$dpq==3] <- 2
  
  # select the subset of the depth profile passing the bulk depth requirement
  # cut down dptab to the max value in g$dp (+1 because 1 corresponds to dp=0)
  rowqs <- cut(0:(nrow(dptab)-1), qbreaks, include.lowest=T, labels=F)
  rowqs[rowqs==3] <- 2
  
  s <- s[dpq %in% qstouse]
  g <- g[dpq %in% qstouse]
  sfp <- results@static.filter.params[[mt]]
  
  
  gbp.per.genome <- get.gbp.by.genome(results) 

  if (mt == 'indel') {
  ncalls <- sapply(qstouse, function(q) sum(s[dpq == q]$condition.pass, na.rm=TRUE))
  } else { 
  ncalls <- sapply(qstouse, function(q) sum(s[dpq == q]$pass, na.rm=TRUE)) 
  }
  
  callable.bp <- sapply(split(dptab[,-(1:sfp$min.bulk.dp)], rowqs), sum)
  callable.sens <- sapply(qstouse, function(q) mean(g[bulk.dp >= sfp$min.bulk.dp & dpq == q & resampled.training.site == TRUE]$training.pass, na.rm=TRUE)) 
  ret <- data.frame(ncalls=ncalls, callable.sens=callable.sens, callable.bp=callable.bp)
  ret$callable.burden <- ret$ncalls / ret$callable.sens
  # dividing by 2 makes it haploid gb
  ret$rate.per.gb <- ret$callable.burden / ret$callable.bp * 1e9/2
  ret$burden <- ret$rate.per.gb * gbp.per.genome
  ret$somatic.sens <- ret$ncalls / ret$burden
  ret$pre.genotyping.burden <- pre.geno.burden
  
  ret$unsupported.filters <- sfp$max.bulk.alt > 0 
  mb <- ret[2,]

  fileConn <- file(paste0(id,'.output'))
   cat(sprintf("# %6s",id),  file=fileConn)
   cat(sprintf("# %6s: %6d somatic,   %0.1f%% sens,   %0.3f callable Gbp,   %0.1f muts/haploid Gbp\n",
                  mt, mb$ncalls, 100*mb$callable.sens, mb$callable.bp/1e9, mb$rate.per.gb), file=fileConn)
  close(fileConn)
}), muttypes)







args=commandArgs(trailingOnly=TRUE)

library(scan2)
tab_file = args[1]
input_yaml = args[2]
sample_name= args[3]
output_prefix= args[4]
mut=args[5]

mutation_tab <- read.csv(tab_file, header=T)
yaml <- yaml::read_yaml(input_yaml)
output_file <- paste0(output_prefix,'_',mut,'.vcf')
ref.genome <- yaml$ref
chrs <- yaml$chrs
fai <- read.table(paste0(ref.genome, '.fai'), sep='\t',stringsAsFactors=F)

f <- file(output_file, 'w')
vcf.header <- c(
        '##fileformat=VCFv4.0',
        '##source=SCAN2',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
vcf.header <- c(vcf.header,
            sprintf('##reference=%s', yaml$ref),
            sprintf('##contig=<ID=%s,length=%d>', fai[,1], fai[,2]))
vcf.header <- c(vcf.header,
            paste(c("#CHROM", "POS", 'ID', 'REF', 'ALT', 'QUAL', 'FILTER',
                'INFO', 'FORMAT', sample_name), collapse='\t'))

writeLines(vcf.header, con=f)
    s <- mutation_tab[!is.na(mutation_tab$pass) & !is.na(mutation_tab$chr) & mutation_tab$pass,]
    s <- do.call(rbind, lapply(chrs, function(chr) {
        ss <- s[s$chr==chr,]
        ss[order(ss$pos),]
    }))
    if (nrow(s) > 0) {
        writeLines(paste(s$chr, s$pos, s$dbsnp, s$refnt, s$altnt,
            '.', 'PASS', '.', 'GT', s[[sample_name]], sep='\t'),
            con=f)
    }
close(f)


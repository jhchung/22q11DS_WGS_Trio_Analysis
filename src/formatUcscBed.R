# Format UCSC bed file for use with vcftools

args <- commandArgs(trailingOnly = TRUE)

bedFile <- args[1]
outputFile <- args[2]

bedData <- read.table(bedFile, sep = "\t", header = TRUE)
bedData <- bedData[, c(1, 2, 3)]
bedData$chrom <- gsub("chr", "", bedData$chrom)

write.table(bedData, file = outputFile, col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
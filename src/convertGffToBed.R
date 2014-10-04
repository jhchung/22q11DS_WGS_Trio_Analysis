# R script to extract chromosome, start and end columns from a gff file and save as bed
options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)

gffFile <- args[1]
bedFile <- args[2]
gffData <- read.table(gffFile, sep = "\t", comment.char = "#", header = FALSE)
gffDataGene <- gffData[(gffData$V3 == "gene" | 
                     gffData$V3 == "ncRNA"), ]
gffDataExon <- gffData[(gffData$V3 == "exon"), ]


bedData <- gffDataGene[, c(1, 4, 5)]
bedData <- unique(bedData)
names(bedData) <- c("chrom", "chromStart", "chromEnd")
write.table(bedData, file = bedFile, row.names = FALSE, col.names = TRUE,
            sep = "\t", quote = FALSE)

bedExon <- gffDataExon[, c(1, 4, 5)]
bedExon <- unique(bedExon)
names(bedExon) <- c("chrom", "chromStart", "chromEnd")
write.table(bedExon, file = file.path(dirname(bedFile), paste("exon", basename(bedFile), sep = "")),
            row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
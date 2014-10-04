# Combine variant counts
library(plyr)
options(stringsAsFactors = FALSE)

# input_files <- c("C:/Users//Jonathan/Dropbox/WGS paper for Barrie/Sequencing Figures and tables for Results/Figure_5_TBX1_pathway_genes/Tbx1_pathway_genes_SNP/manual/BM1452.tbx1_genes.var.stat.txt",
#                  "C:/Users//Jonathan/Dropbox/WGS paper for Barrie/Sequencing Figures and tables for Results/Figure_5_TBX1_pathway_genes/Tbx1_pathway_genes_SNP/manual/BM1453.tbx1_genes.var.stat.txt")
args <- commandArgs(trailingOnly=TRUE)
output_file <- args[1]
input_files <- args[2:length(args)]

variant_counts <- llply(input_files, read.table, header = FALSE, sep = "\t")

merged_counts <- Reduce(function(...) merge(..., by = "V1", all=T), variant_counts)

write.table(merged_counts, file = output_file,
            quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE)
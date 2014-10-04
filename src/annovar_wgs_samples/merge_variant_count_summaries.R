library(plyr)
library(reshape)

parse_filename <- function(filename){
  if (!grepl("sample.BM", filename)){
    return("All")
  } else {
    filename <- basename(filename)
    filename_split <- strsplit(filename, "\\.")[[1]]
    bmid <- grep("BM", filename_split, value = TRUE)
    sample_extension <- grep("00", filename_split, value = TRUE)
    return(paste(bmid, sample_extension, sep = "."))
  }
}

add_sample_name_to_df <- function(data_frames, sample_names){
  for (i in 1:length(sample_names)){
    names(data_frames[[i]]) <- c("Variant type", sample_names[i])
  }
  return(data_frames)
}


args <- commandArgs(trailingOnly=TRUE)

snp_output <- args[1]
indel_output <- args[2]
input_files <- args[3:length(args)]

snp_files <- grep("wgs_phased_snp", input_files, value = TRUE)
indel_files <- grep("wgs_phased_indel", input_files, value = TRUE)

snp_dat <- lapply(snp_files, read.table, sep = "\t", header = FALSE)
indel_dat <- lapply(indel_files, read.table, sep = "\t", header = FALSE)

snp_names <- lapply(snp_files, parse_filename)
indel_names <- lapply(indel_files, parse_filename)

snp_dat <- add_sample_name_to_df(snp_dat, snp_names)
indel_dat <- add_sample_name_to_df(indel_dat, indel_names)

snp_merged <- merge_recurse(snp_dat)
indel_merged <- merge_recurse(indel_dat)

write.table(snp_merged, file = snp_output, sep = "\t",
            col.names = TRUE, row.names = FALSE, quote = FALSE)
write.table(indel_merged, file = indel_output, sep = "\t",
            col.names = TRUE, row.names = FALSE, quote = FALSE)
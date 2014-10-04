# Use this script to extract control sample names from FAM file.

args <- commandArgs(trailingOnly = TRUE)

famFile <- args[1]
outputDir <- args[2]
famData <- read.table(famFile, header = FALSE)
controlSamples <- famData[(famData$V6 == 1), ]
caseSamples <- famData[(famData$V6 == 2), ]

controlOutputFile <- file.path(outputDir , "wesControls.txt")
caseOutputFile <- file.path(outputDir, "wesTOF.txt")
write(controlSamples$V2, file = controlOutputFile, sep = "\n")
write(caseSamples$V2, file = caseOutputFile, sep = "\n")
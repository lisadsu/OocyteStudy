# Process output from HTSeq-count and turn it into a count table

source("http://bioconductor.org/biocLite.R")
library(Rsamtools)
library(DESeq)
library(edgeR)
library(sva)
library(ALL)
library(annotate)
library(GO.db)
library(Homo.sapiens)
library(goseq)
library(rtracklayer)
library(GenomicRanges)
library(ape)
library(rgl)
library(reshape2)
library(gplots)
library(ggplot2)
library(RColorBrewer)
library(pcaMethods)
library(plotrix)
library(scatterplot3d)
library(matrixStats)
library(clValid)

if (Sys.info()["sysname"] == "Linux") {
  source("/host/Users/Livia/Desktop/IVF/RnaSeqAnalysis/RawData/Scripts/R/UsefulFunctions.R")
  baseDataDirectory <- "/host/Users/Livia/Desktop/IVF/RnaSeqAnalysis/RawData/HumanHiSeq/HTSeq_Count_Out"
} else {
  source("C:/Users/Livia/Desktop/IVF/RnaSeqAnalysis/RawData/Scripts/R/UsefulFunctions.R")
  baseDataDirectory <- "C:/Users/Livia/Desktop/IVF/RnaSeqAnalysis/RawData/HumanHiSeq/HTSeq_Count_Out"
}


##################################################################################################
# Prepare count table (output from HTSeq-count)                                                  #
##################################################################################################

fileNameList <- c("S1", "S2", "S3", "S4", "S5", "S6", "S7", "S8", "S9", "S10", "S11", "S12", "S13", "S14", "S15", "S16", "S17", "S18", "S19", "S20", "S21", "S22")

# read in count tables for all samples
S1_count <- read.table(paste(baseDataDirectory, "/", fileNameList[1], "_paired_count.txt", sep = ""), header = FALSE, row.names=1, col.names=c("name", "S1"))
S2_count <- read.table(paste(baseDataDirectory, "/", fileNameList[2], "_paired_count.txt", sep = ""), header = FALSE, row.names=1, col.names=c("name", "S2"))
S3_count <- read.table(paste(baseDataDirectory, "/", fileNameList[3], "_paired_count.txt", sep = ""), header = FALSE, row.names=1, col.names=c("name", "S3"))
S4_count <- read.table(paste(baseDataDirectory, "/", fileNameList[4], "_paired_count.txt", sep = ""), header = FALSE, row.names=1, col.names=c("name", "S4"))
S5_count <- read.table(paste(baseDataDirectory, "/", fileNameList[5], "_paired_count.txt", sep = ""), header = FALSE, row.names=1, col.names=c("name", "S5"))
S6_count <- read.table(paste(baseDataDirectory, "/", fileNameList[6], "_paired_count.txt", sep = ""), header = FALSE, row.names=1, col.names=c("name", "S6"))
S7_count <- read.table(paste(baseDataDirectory, "/", fileNameList[7], "_paired_count.txt", sep = ""), header = FALSE, row.names=1, col.names=c("name", "S7"))
S8_count <- read.table(paste(baseDataDirectory, "/", fileNameList[8], "_paired_count.txt", sep = ""), header = FALSE, row.names=1, col.names=c("name", "S8"))
S9_count <- read.table(paste(baseDataDirectory, "/", fileNameList[9], "_paired_count.txt", sep = ""), header = FALSE, row.names=1, col.names=c("name", "S9"))
S10_count <- read.table(paste(baseDataDirectory, "/", fileNameList[10], "_paired_count.txt", sep = ""), header = FALSE, row.names=1, col.names=c("name", "S10"))
S11_count <- read.table(paste(baseDataDirectory, "/", fileNameList[11], "_paired_count.txt", sep = ""), header = FALSE, row.names=1, col.names=c("name", "S11"))
S12_count <- read.table(paste(baseDataDirectory, "/", fileNameList[12], "_paired_count.txt", sep = ""), header = FALSE, row.names=1, col.names=c("name", "S12"))
S13_count <- read.table(paste(baseDataDirectory, "/", fileNameList[13], "_paired_count.txt", sep = ""), header = FALSE, row.names=1, col.names=c("name", "S13"))
S14_count <- read.table(paste(baseDataDirectory, "/", fileNameList[14], "_paired_count.txt", sep = ""), header = FALSE, row.names=1, col.names=c("name", "S14"))
S15_count <- read.table(paste(baseDataDirectory, "/", fileNameList[15], "_paired_count.txt", sep = ""), header = FALSE, row.names=1, col.names=c("name", "S15"))
S16_count <- read.table(paste(baseDataDirectory, "/", fileNameList[16], "_paired_count.txt", sep = ""), header = FALSE, row.names=1, col.names=c("name", "S16"))
S17_count <- read.table(paste(baseDataDirectory, "/", fileNameList[17], "_paired_count.txt", sep = ""), header = FALSE, row.names=1, col.names=c("name", "S17"))
S18_count <- read.table(paste(baseDataDirectory, "/", fileNameList[18], "_paired_count.txt", sep = ""), header = FALSE, row.names=1, col.names=c("name", "S18"))
S19_count <- read.table(paste(baseDataDirectory, "/", fileNameList[19], "_paired_count.txt", sep = ""), header = FALSE, row.names=1, col.names=c("name", "S19"))
S20_count <- read.table(paste(baseDataDirectory, "/", fileNameList[20], "_paired_count.txt", sep = ""), header = FALSE, row.names=1, col.names=c("name", "S20"))
S21_count <- read.table(paste(baseDataDirectory, "/", fileNameList[21], "_paired_count.txt", sep = ""), header = FALSE, row.names=1, col.names=c("name", "S21"))
S22_count <- read.table(paste(baseDataDirectory, "/", fileNameList[22], "_paired_count.txt", sep = ""), header = FALSE, row.names=1, col.names=c("name", "S22"))

# combine tables into one
# first remove extra rows from end of HTSeq-count files
nRowToUse = (nrow(S2_count) - 5)
allCounts <- data.frame(matrix(ncol = (length(fileNameList)), nrow = nRowToUse))
names(allCounts) <- fileNameList
row.names(allCounts) <- row.names(S2_count)[1:nRowToUse]

allCounts[,"S1"] <- S1_count[1:nRowToUse,names(S1_count)[1]]
allCounts[,"S2"] <- S2_count[1:nRowToUse,names(S2_count)[1]]
allCounts[,"S3"] <- S3_count[1:nRowToUse,names(S3_count)[1]]
allCounts[,"S4"] <- S4_count[1:nRowToUse,names(S4_count)[1]]
allCounts[,"S5"] <- S5_count[1:nRowToUse,names(S5_count)[1]]
allCounts[,"S6"] <- S6_count[1:nRowToUse,names(S6_count)[1]]
allCounts[,"S7"] <- S7_count[1:nRowToUse,names(S7_count)[1]]
allCounts[,"S8"] <- S8_count[1:nRowToUse,names(S8_count)[1]]
allCounts[,"S9"] <- S9_count[1:nRowToUse,names(S9_count)[1]]
allCounts[,"S10"] <- S10_count[1:nRowToUse,names(S10_count)[1]]
allCounts[,"S11"] <- S11_count[1:nRowToUse,names(S11_count)[1]]
allCounts[,"S12"] <- S12_count[1:nRowToUse,names(S12_count)[1]]
allCounts[,"S13"] <- S13_count[1:nRowToUse,names(S13_count)[1]]
allCounts[,"S14"] <- S14_count[1:nRowToUse,names(S14_count)[1]]
allCounts[,"S15"] <- S15_count[1:nRowToUse,names(S15_count)[1]]
allCounts[,"S16"] <- S16_count[1:nRowToUse,names(S16_count)[1]]
allCounts[,"S17"] <- S17_count[1:nRowToUse,names(S17_count)[1]]
allCounts[,"S18"] <- S18_count[1:nRowToUse,names(S18_count)[1]]
allCounts[,"S19"] <- S19_count[1:nRowToUse,names(S19_count)[1]]
allCounts[,"S20"] <- S20_count[1:nRowToUse,names(S20_count)[1]]
allCounts[,"S21"] <- S21_count[1:nRowToUse,names(S21_count)[1]]
allCounts[,"S22"] <- S22_count[1:nRowToUse,names(S22_count)[1]]


exptDesign <- data.frame(row.names = colnames(allCounts), 
                         condition = c("bad1", "good1", "bad1", "bad1", "bad1", "good1", "bad1", "good1", "good1", "bad1", "bad2", "bad2", "bad2", "bad2", 
                                       "bad2", "bad2", "good2", "good2", "good2", "good2", "good2", "bad2"), 
                         libType = c(rep("pe", 22)))
condition = exptDesign$condition


# assemble count data set
conditionAll = c("bad", "good", "bad", "bad", "bad", "good", "good", "good", "good", "bad", "bad", "bad", "bad", "bad", 
                 "good", "bad", "good", "good", "good", "good", "good", "bad")

conditionDiff = c("bad1", "good1", "bad1", "bad1", "bad1", "good1", "good1", "good1", "good1", "bad1", "bad2", "bad2", "bad2", "bad2", 
                  "good2", "bad2", "good2", "good2", "good2", "good2", "good2", "bad2")
decDist = c(4.36, 0.7, 2.09, 3.19, 3.38, 1.08, 1.63, 0.50, 1.12, 3.90, 6.31, 4.04, 4.63, 2.87, 1.72, 2.55, 1.07, 0.94, 1.08, 1.46, 1.82, 2.57)
embryosSamePatient = c(3,4,6,7,8,10,11,12,13,14,15,16,17,18,20,21,22)
exptNum = c(1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2)
includePatient = c(FALSE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, TRUE)




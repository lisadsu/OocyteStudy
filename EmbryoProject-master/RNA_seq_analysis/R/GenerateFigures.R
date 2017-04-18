# load all the usual stuff

rm(list = ls())


if (Sys.info()["sysname"] == "Linux") {
  source("/host/Users/Livia/Desktop/IVF/RnaSeqAnalysis/RawData/Scripts/R/LoadHTSeqCountData_Human_HiSeq.R")
  # read in RPKM data
  conditionCombined = c(condition, 5,5,5)
  allRPKM = read.table("/host/Users/Livia/Desktop/IVF/RnaSeqAnalysis/RawData/HumanHiSeq/allHumanEmbryoRPKM.txt")
  #yanData = read.table("/host/Users/Livia/Desktop/IVF/RnaSeqAnalysis/RawData/HumanHiSeq/Yan2013.csv", header = TRUE, row.names = 1)
  combinedRPKMtable = read.table("/host/Users/Livia/Desktop/IVF/RnaSeqAnalysis/RawData/HumanHiSeq/combinedEmbryoRPKM.txt")
} else {
  source("C:/Users/Livia/Desktop/IVF/RnaSeqAnalysis/RawData/Scripts/R/LoadHTSeqCountData_Human_HiSeq.R")
  # read in RPKM data
  conditionCombined = c(condition, 5,5,5)
  allRPKM = read.table("C:/Users/Livia/Desktop/IVF/RnaSeqAnalysis/RawData/HumanHiSeq/allHumanEmbryoRPKM.txt")
  #yanData = read.table("C:/Users/Livia/Desktop/IVF/RnaSeqAnalysis/RawData/HumanHiSeq/Yan2013.csv", header = TRUE, row.names = 1)
  combinedRPKMtable = read.table("C:/Users/Livia/Desktop/IVF/RnaSeqAnalysis/RawData/HumanHiSeq/combinedEmbryoRPKM.txt")
}


conditionAll = c("bad", "good", "bad", "bad", "bad", "good", "good", "good", "good", "bad", "bad", "bad", "bad", "bad", 
                 "good", "bad", "good", "good", "good", "good", "good", "bad")
decDist = c(4.36, 0.7, 2.09, 3.19, 3.38, 1.08, 1.63, 0.50, 1.12, 3.90, 6.31, 4.04, 4.63, 2.87, 1.72, 2.55, 1.07, 0.94, 1.08, 1.46, 1.82, 2.57)
k1 = c(.2334, .3072, .2841, .2783, .2680, .3167, .2678, .3024, .3092, .2540, .2489, .2561, .2467, .3481, .2742, .2866, .3189, .2916, .3190, .3160, .3360, .3431)
n1 = c(.9846, .6002, .2841, .9604, .9848, .5634, .6229, .5988, .5586, .5216, 1.5256, 1.0629, 1.1369, .4737, .7248, .9019, .6379, .5296, .5620, .4926, .6242, .4885)
embryosSamePatient = c(3,4,6,7,8,10,11,12,13,14,15,16,17,18,20,21,22)


combat.edata = read.table(paste(baseDataDirectory, "/edgeR/ComBat_SVA/combat_edata.txt", sep = ""))
DEnames.qvals = read.table(paste(baseDataDirectory, "/edgeR/ComBat_SVA/DEnames_qvals.txt", sep = ""))
qValuesComBat = DEnames.qvals[,1]
names(qValuesComBat) = row.names(DEnames.qvals)
DEnames.qvals = qValuesComBat
logFC = rowMeans(combat.edata[,conditionAll[embryosSamePatient] == "good"]) - 
  rowMeans(combat.edata[,conditionAll[embryosSamePatient] == "bad"])

# re-plot phylo
hmcol = c("#009b00", "#009b00", "#00009b", "#00009b", "#00009b", "#00009b", "#00009b")
hc = hclust(dist(t(combat.edata)))
# hc = hclust(dist(t(log2(cpm(yDiff)+1))))
#plot(as.phylo(hc), tip.color=hmcol[lane-1])
plot(as.phylo(hc), tip.color=hmcol[ceiling(decDist[embryosSamePatient])])


# get GO data
if (Sys.info()["sysname"] == "Linux") {
  source("/host/Users/Livia/Desktop/IVF/RnaSeqAnalysis/RawData/Scripts/R/getGOdata.R")
} else {
  source("C:/Users/Livia/Desktop/IVF/RnaSeqAnalysis/RawData/Scripts/R/getGOdata.R")
}



# find genes most correlated to k1
k1.corr = apply(combat.edata, 1, cor, y = k1[embryosSamePatient])
k1.corr = sort(k1.corr, decreasing = TRUE)

n1.corr = apply(combat.edata, 1, cor, y = log2(n1[embryosSamePatient]))
n1.corr = sort(n1.corr, decreasing = TRUE)

decDist.corr = apply(combat.edata, 1, cor, y = decDist[embryosSamePatient])
decDist.corr = sort(decDist.corr, decreasing = TRUE)

embryoViability = conditionAll[embryosSamePatient]

dataToPlot = as.data.frame(cbind(combat.edata[names(k1.corr)[16],],k1[embryosSamePatient]))
names(dataToPlot) = c("expr", "k1")
ggplot(dataToPlot, aes(expr,k1, color=embryoViability)) + 
  geom_point(size=5) + geom_text(hjust=-.5, vjust=-.5,aes(label=as.character(1:22)[embryosSamePatient])) +
  theme(legend.text = element_text(size = 20), legend.title = element_text(size=20), 
        axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20))

# write.table(names(decDist.corr)[decDist.corr < -0.7], paste(baseDataDirectory, "/edgeR/mechanicalParamCorr/lowCorrDecDist.txt", sep = ""),
#             row.names = FALSE, quote = FALSE, col.names = FALSE)


###################################################################### 
# Plot log fold change box and scatter plots for selected genes      #
# Subplots for Figure 3D, and Figure 4A in paper #####################
###################################################################### 

# cell cycle (meiosis, mitosis, chromosome segregation)
genesToPlot = c("MAD2L1", "BUB1", "CDC25B", "SYCP3", "ANAPC4", "KIF2C", "CIT", "SKA1", "POGZ", "PHF13",  
                "SMC2", "SMC4", "SMARCAD1", "PPP2R1A", "TOP1", "TTN", "PTTG1", "CDK1", "CCNA2")

# genes involved in DNA repair, telomere maintenance
genesToPlot = c("TERF1", "XRCC6", "NBN",  "FANCB",  "KIN", "XAB2", 
                "ABL1", "HMGB2", "MORF4L1", "PARP2", "POLN", "PAPD7", "PCNA", "PRMT6") #"WRN", "WDR33","FANCI","FBXO18", 

# genes involved in regulation of transcription and chromatin modification
genesToPlot = c("CBX2", "ACTL6A", "USP16", "SMARCA2", "PHF20", "CENPI", "LRWD1", "PIM3", "SYCP3", 
                "TET3", "PRMT6", "DNMT3B", "HDAC1", "TERF1", "YY1")

# fertilization
genesToPlot = c("ZP3", "ZP4", "ITPR1", "STX2", "PLCB1", "PLCB3", "CD9", "PRKCG", "OOEP", "CD81", "PLCZ1")
# genesToPlot = c("ZP3", "ZP4", "ITPR1", "STX2", "PLCB1", "PLCD4", "CD9", "SRCIN1", "PRKCG", "OOEP", "PLCZ1")


dataToPlot = as.data.frame(matrix(ncol = 4, nrow = (length(genesToPlot)*17)))
names(dataToPlot) = c("sample", "condition", "gene", "value")

for (i in 1:length(genesToPlot)){
  
  dataToPlot[((i-1)*17 + 1):(i*17),] = cbind(names(as.data.frame(combat.edata)), conditionAll[embryosSamePatient], 
                                                   rep(genesToPlot[i],17), as.numeric(combat.edata[genesToPlot[i],]))
  
  
}

dataToPlot$value = as.numeric(dataToPlot$value)
dataToPlot$valueDiff = NA

# calculate logFC from average good embryo
for (i in 1:length(genesToPlot)){
  
  avgGood = median(dataToPlot[((i-1)*17 + 1):(i*17),"value"][conditionAll[embryosSamePatient] == "good"])
  dataToPlot[((i-1)*17 + 1):(i*17),"valueDiff"] = dataToPlot[((i-1)*17 + 1):(i*17),"value"] - avgGood
  
}

dataToPlot$condition[dataToPlot$condition == "good"] = "viable"
dataToPlot$condition[dataToPlot$condition == "bad"] = "nonviable"
embryo.colors = data.frame(good = "#009b00", bad = "#00009b")

ggplot(as.data.frame(dataToPlot)) + 
  geom_boxplot(aes(x=gene, y=valueDiff, col=condition, fill=condition), alpha = .3, position=position_dodge(width=.75), outlier.size = 0) +
  geom_point(aes(x=gene, y=valueDiff, col=condition, fill=condition), position = position_jitterdodge(dodge.width=.75)) +
  geom_abline(aes(intercept = 0, slope = 0), linetype = 2) +
  theme(legend.text = element_text(size = 16), legend.title = element_blank(), plot.title = element_text(size = 20, vjust = 2),
        axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20, vjust = 0.3),
        axis.text.x = element_text(angle=45, size = 16, vjust = 1, hjust = 1, colour = 'black'), axis.text.y = element_text(size = 20),
        panel.background = element_rect(fill = 'white', colour = 'black'))+#, legend.position = c(0,1), 
 #      legend.background = element_rect(colour = 'black'), legend.justification = c(0,1)) +
  xlab(NULL) +
  ylab("Log-fold-change") +
  ylim(-2, 2) +
  ggtitle("Genes Important for Fertilization") +
  scale_colour_manual(values = c("#00009b", "#009b00")) +
  scale_fill_manual(values = c("#00009b", "#009b00"))






######################################
# write data in format cytoscape will like
# 1. entrez gene IDs, 2. official gene symbol 3. expr_viable  4. expr_nonviable  5. logFC(good-bad)  6. qvalue

avgGood = rowMeans(combat.edata[namesSymbol,conditionAll[embryosSamePatient] == "good"])
avgBad = rowMeans(combat.edata[namesSymbol,conditionAll[embryosSamePatient] == "bad"])
logFC = avgGood - avgBad
qvals = DEnames.qvals[namesSymbol]

dataCytoscape = data.frame(namesEntrez, namesSymbol, avgGood, avgBad, (-1*logFC), qvals)
names(dataCytoscape) = c("entrezID", "geneSymbol", "goodExpr", "badExpr", "-logFC", "qval")
dataCytoscape = dataCytoscape[DEnames.qvals < .05,]
# dataCytoscape = data.frame(namesSymbol[DEnames.qvals < .05], logFC[DEnames.qvals < .05], DEnames.qvals[DEnames.qvals < .05])
write.table(dataCytoscape, file = paste(baseDataDirectory, "/edgeR/DEGenesTable_logFCinverted.txt", sep = ""), sep = "\t",
                         row.names = FALSE, quote = FALSE, col.names = TRUE)                           


write.table(namesEntrez, file = paste(baseDataDirectory, "/edgeR/Cytoscape/namesEntrezCytoscape.txt", sep = ""),
            row.names = FALSE, quote = FALSE, col.names = FALSE)




#########################################
# Figure 2D 
#########################################

embryoMech = as.data.frame(t(read.table("C:/Users/Livia/Dropbox/Embryo Mechanics outline shared/Data/embryoMechanics/humanEmbryoParams.txt")))
embryoDecDist = as.data.frame(t(read.table("C:/Users/Livia/Dropbox/Embryo Mechanics outline shared/Data/embryoMechanics/humanEmbryoDecDist.txt")))
names(embryoMech) = c("mN", "k1N", "n1N", "k0N")
names(embryoDecDist) = "decDist"
row.names(embryoMech) = 1:89
row.names(embryoDecDist) = 1:89
embryoMech = cbind(embryoMech, embryoDecDist)

embryoMech = embryoMech[,c("mN", "k1N", "n1N", "k0N")]#, "decDist")] # take out tN
embryoMech[,c("k1N", "n1N", "k0N")] = log(embryoMech[,c("k1N", "n1N", "k0N")]) # take log scale for n1 and k0
# 
colMeanVars = matrix(rep(apply(embryoMech[embryoMech[,"mN"] == 4,][,c("k1N", "n1N", "k0N")], 2, median), each = 89), nrow = 89)
embryoMech[,c("k1N", "n1N", "k0N")] = embryoMech[,c("k1N", "n1N", "k0N")] - colMeanVars # subtract avg good embryo

embryo.colors = data.frame(viable = "#009b00", nonviable = "#00009b")
embryoMech = melt(embryoMech, id = "mN")

embryoMech[,"mN"] = sub(1, "nonviable", embryoMech[,"mN"])
embryoMech[,"mN"] = sub(4, "viable", embryoMech[,"mN"])


ggplot(as.data.frame(embryoMech)) + 
  geom_boxplot(aes(x=variable, y=value, col=mN, fill=mN), alpha = .3, position=position_dodge(width=.75), outlier.colour=NULL) +
  geom_point(aes(x=variable, y=value, col=mN, fill=mN), position = position_jitterdodge(dodge.width=.75)) +
  geom_abline(aes(intercept = 0, slope = 0), linetype = 2) +
  theme(legend.text = element_text(size = 16), legend.title = element_blank(), plot.title = element_text(size = 20, vjust = 2),
        axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 16, vjust = 0.3),
        axis.text.x = element_text(angle=0, size = 16, vjust = 1, hjust = 1, colour = 'black'), axis.text.y = element_text(size = 20),
        panel.background = element_rect(fill = 'white', colour = 'black'), legend.position = c(0,1), 
        legend.background = element_rect(colour = 'black'), legend.justification = c(0,1)) + 
  xlab(NULL) +
  ylab("Change from Median Viable Embryo") +
#   ylim(-0.75, 4) +
#  scale_y_continuous(trans = "log") + 
  ggtitle("Mechanical Parameter Distributions") +
  scale_colour_manual(values = c("#00009b", "#009b00")) +
  scale_fill_manual(values = c("#00009b", "#009b00"))



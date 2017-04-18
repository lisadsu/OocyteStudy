# Run edgeR analysis of raw count data 

rm(list = ls())


if (Sys.info()["sysname"] == "Linux") {
  source("/host/Users/Livia/Desktop/IVF/RnaSeqAnalysis/RawData/Scripts/R/LoadHTSeqCountData_Human_HiSeq.R")
} else {
  source("C:/Users/Livia/Desktop/IVF/RnaSeqAnalysis/RawData/Scripts/R/LoadHTSeqCountData_Human_HiSeq.R")
}


# run edgeR
yAll = DGEList(counts = allCounts, group = conditionAll)
yAll = yAll[rowSums(allCounts > 20) > 4,]
yAll = calcNormFactors(yAll)
yAll = estimateDisp(y = yAll, design = model.matrix(~0 + factor(conditionAll)))


# create model matrix and null matrix
patient = c(1,2,3,3,4,5,5,6,7,6,8,9,8,8,9,8,10,10,11,12,10,12) - 1
batch = c(1,1,2,2,3,3,4,5,4,4,5,4,6,6,7,6,7)

exptInfo = data.frame(conditionAll[embryosSamePatient], exptNum[embryosSamePatient])
mod = model.matrix(~as.factor(conditionAll)[embryosSamePatient], data = exptInfo)
mod0 = model.matrix(~1, data = exptInfo)
mod1 = model.matrix(~as.factor(patient[embryosSamePatient]))


# use ComBat to adjust for batch effects
combat.edata = as.data.frame(ComBat(dat = log2(cpm(yAll)[,embryosSamePatient]+1), batch = batch, mod = mod, 
                      par.prior = TRUE, prior.plots = TRUE))

# filter out genes with < 1 cpm average expression
combat.edata = combat.edata[rowMeans(combat.edata) > 0,]

# calculate pvalues and qvalues
pValuesComBat = f.pvalue(as.matrix(combat.edata), mod, mod0)
qValuesComBat = p.adjust(pValuesComBat, method = "BH")

# write out adjusted qvalue data
write.table(as.character(names(qValuesComBat[qValuesComBat < .01])), paste(baseDataDirectory, "/edgeR/ComBat_SVA/qVals_lt_pt01_ComBat.txt", sep = ""), row.names=FALSE, quote=FALSE)
# write.table(as.character(names(qValuesSv[qValuesSv < .1])), paste(baseDataDirectory, "/edgeR/ComBat_SVA/qVals_lt_pt05_SVA.txt", sep = ""), row.names=FALSE, quote=FALSE)
# write.table(combat.edata, paste(baseDataDirectory, "/edgeR/ComBat_SVA/combat_edata.txt", sep = ""))
# write.table(DEnames.qvals, paste(baseDataDirectory, "/edgeR/ComBat_SVA/DEnames_qvals.txt", sep = ""))


# make a few different PCA plots
# Figure 3B
x = combat.edata
x = x[!is.na(rowSums(x)),]
pca = prcomp(t(x))
fac = factor(conditionDiff[embryosSamePatient])
colours = brewer.pal(nlevels(fac), "Paired") 
xyplot(PC2 ~ PC1, groups = fac, data = as.data.frame(pca$x), pch = 16, cex = 2, aspect = "fill", 
       col = colours, main = draw.key(key = list(rect = list(col = colours), text = list(levels(fac)), rep = FALSE)))

colours = c("#00009b", "#009b00")
plot3d(pca$x[,"PC1"], pca$x[,"PC2"], pca$x[,"PC3"], col = colours[as.numeric(as.factor(conditionAll[embryosSamePatient]))], 
       size = 20, box = FALSE, xlab = "", ylab = "", zlab = "")
decorate3d(xlim = c(-30, 30), xlab = "", ylab = "", zlab = "", box = FALSE)

ggplot(as.data.frame(pca$x), aes(PC1,PC2, color=conditionAll[embryosSamePatient])) + geom_point(size=5) + 
      geom_text(hjust=-.5, vjust=-.5, aes(label=as.character(1:22)[embryosSamePatient]))



# Figure 3A
# phylo plot
hmcol = c("#009b00", "#009b00", "#00009b", "#00009b", "#00009b", "#00009b", "#00009b")
hc = hclust(dist(t(combat.edata)))
plot(as.phylo(hc), tip.color=hmcol[ceiling(decDist[embryosSamePatient])])


# calculate average logFC
avgGood = rowMeans(combat.edata[,conditionAll[embryosSamePatient] == "good"])
avgBad = rowMeans(combat.edata[,conditionAll[embryosSamePatient] == "bad"]) 
logFC = rowMeans(combat.edata[,conditionAll[embryosSamePatient] == "good"]) - 
                        rowMeans(combat.edata[,conditionAll[embryosSamePatient] == "bad"])
names(logFC) = row.names(combat.edata)


# only keep genes with known identifiers in Entrez catalog
# For the reverse map:
x <- org.Hs.egSYMBOL2EG
# Get the entrez gene identifiers that are mapped to a gene symbol
mapped_genes <- mappedkeys(x)
# Convert to a list
xx <- as.list(x[mapped_genes])

# get entrezID for each gene symbol name
# remove gene symbol names for which no entrez ID is found, and genes with na values
namesSymbol = names(qValuesCombat)
namesEntrez = unlist(sapply(namesSymbol, getEntrezID, mapping = xx))


# list of DE genes 
qValuesComBat = qValuesComBat[names(namesEntrez)]
qValuesComBat = qValuesComBat[!is.na(qValuesComBat)]
exprValues = combat.edata[names(qValuesComBat),]
logFC = logFC[names(qValuesComBat)]
DEnames.qvals = qValuesComBat[qValuesComBat < .01]

namesEntrezFilt = namesEntrez[names(qValuesComBat)]
qValuesComBatEntrez = qValuesComBat
names(qValuesComBatEntrez) = namesEntrezFilt

# reverse map from nums to names
namesEntrezFiltReverse = names(namesEntrezFilt)
names(namesEntrezFiltReverse) = namesEntrezFilt



# get lists of upregulated and downregulated genes
upGenesDE = names(logFC[(logFC > 0) & (qValuesComBat < .01)])
downGenesDE = names(logFC[(logFC < 0) & (qValuesComBat < .01)])


# Figure 3C
# plot log-fold-change for all genes in catalog
# DE genes are highlighted in blue
exprsToPlot = cbind(rowMeans(exprValues), logFC)
colorMat = rep(1, length(logFC))
colorMat[qValuesComBat < .01] = colorMat[qValuesComBat < .01] + 1
colorsForPlot = c("#636363", "#3182bd")
plot(exprsToPlot, main="Smear Plot", xlab = "mean expr", ylab = "log fold change (viable vs nonviable) ", col = colorsForPlot[colorMat])
abline(h = 0, lwd = 2, col = "black", lty = "dashed")
abline(h = 1, lwd = 1, col = "black")
abline(h = -1, lwd = 1, col = "black")



# package goseq
# get gene ontology information for DE genes
pwf = nullp(qValuesComBatEntrez < .01, "hg18", "knownGene")
GO.wall = goseq(pwf, "hg18", "knownGene")
go = getgo(names(qValuesComBatEntrez[qValuesComBatEntrez < .01]), "hg18", "knownGene")


# list of cell cycle GO categories
c6 = c("GO:0030163", "GO:0042176", "GO:0070271",  "GO:0030163", "GO:0008104")

# finds indices in list of ALL genes which contain desired GO terms
# so we take the intersection of genes containing desired GO terms, and genes that are differentially expressed at q < 0.01
genes.ind = sapply(go, findGenesWithGoTerms, desiredGoList = c6)
genes = intersect(namesEntrezFiltReverse[genes.ind], names(DEnames.qvals))




# go through all gene ontology terms (1 through 99,999) and find genes associated with that term
# see if the DE genes among those genes are mostly up or down regulated 

percent.misregulated = rep(0,75500)
numGenesInCategory = rep(0,75500)

for (i in 1:32000){
  
  iterChar = as.character(i)
  addString = paste(as.character(rep(0,5-nchar(iterChar))), sep = "", collapse = "")
  currGOstring = paste("GO:00", addString, iterChar, sep="")
  
  genes.ind = sapply(go, findGenesWithGoTerms, desiredGoList = currGOstring)
  genes = namesSymbol[genes.ind]
  numGenesInCategory[i] = length(genes)
  names(percent.misregulated)[i] = currGOstring
  percent.misregulated[i] = length(DEnames.qvals[genes][DEnames.qvals[genes] < .01])/length(DEnames.qvals[genes])
  
#   # now take average and standard deviation logFC for all those genes
#   logFC = rowMeans(combat.edata[genes,conditionAll[embryosSamePatient] == "good"]) - 
#     rowMeans(combat.edata[genes,conditionAll[embryosSamePatient] == "bad"])

  if (i %% 1000 == 0){
    print(i)
    print(currGOstring)
    print(percent.misregulated[i])
  }

  
}

# sort percent.misregulated
names(numGenesInCategory) = names(percent.misregulated)
percent.mis.sort = sort(percent.misregulated[numGenesInCategory > 50], decreasing = TRUE)
numGenes.sort = numGenesInCategory[names(percent.mis.sort)]




# scatterplot with 2 genes only
allGeneNames = row.names(combat.edata)
embryoViability = conditionAll[embryosSamePatient]
namesSymbol[DEnames.qvals < .05 & namesSymbol < "MB" & namesSymbol > "MA"]
allGeneNames[allGeneNames < "J" & allGeneNames > "I"]
ggplot(as.data.frame(t(combat.edata)), aes(AURKA,CCNA2, color=embryoViability)) + 
  geom_point(size=5) + geom_text(hjust=-.5, vjust=-.5,aes(label=as.character(1:22)[embryosSamePatient])) +
  theme(legend.text = element_text(size = 20), legend.title = element_text(size=20), 
        axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20))



# bar chart with names of all genes in a particular list / GO category
genesToPlot = c("IGF1R", "MAPK1", "MAP2K1", "PRKACA", "YWHAG", "CDK1","MAD2L1", "SKP2", "CCNA2", "CCNE1",
                "PPP2R5A", "ESPL1", "RAD21",  "EME1", "PTTG1", "BTRC", "PPP1CC", "SLK", "PIM1",
                "CDC25C", "FBXO5", "ANAPC4", "SMC2", "SMC3",   "BUB1")
genesToPlot = c("B4GALT1", "CALML4", "CCT2", "CCT5", "CCT6A", "CD9", "CDK1", "CRCP", "ITPR1",
                 "NOX5", "OOEP", "PARK7", "PLCB1", "PRKCA", "PRKCG", "PRKCI", "PLCZ1", "STX2", "WEE2", "ZP3", "ZP4")


goodMeans = rowMeans(exprValues[genesToPlot,embryoViability == "good"])
badMeans = rowMeans(exprValues[genesToPlot,embryoViability == "bad"])
goodStds = sqrt(rowVars(as.matrix(exprValues[genesToPlot, embryoViability == "good"])))
badStds = sqrt(rowVars(as.matrix(exprValues[genesToPlot, embryoViability == "bad"])))

edata.means = as.data.frame(cbind(goodMeans, badMeans))
edata.stds = as.data.frame(cbind(goodStds, badStds))

mp = barplot2(as.matrix(t(edata.means)), beside=TRUE, las = 2, col = c("#009b00", "#00009b"), ylim = c(0,15),
              plot.ci = TRUE, ci.l =t(edata.means) - t(edata.stds), ci.u = t(edata.means) + t(edata.stds))

legend(1,14, c("viable", "nonviable"), fill = c("#009b00", "#00009b"))
box()


########
# Filter genes by highest within-group mean/variance
conditionNum = c(0,2,0,0,0,2,2,2,2,0,1,1,1,1,1,1,3,3,3,3,3,1)
BCV_b1 = as.data.frame(sqrt(rowVars(cpm(yDiff)[,conditionNum == 0])) / rowMeans(cpm(yDiff)[,conditionNum == 0]))
BCV_g1 = as.data.frame(sqrt(rowVars(cpm(yDiff)[,conditionNum == 2])) / rowMeans(cpm(yDiff)[,conditionNum == 2]))
BCV_b2 = as.data.frame(sqrt(rowVars(cpm(yDiff)[,conditionNum == 1])) / rowMeans(cpm(yDiff)[,conditionNum == 1]))
BCV_g2 = as.data.frame(sqrt(rowVars(cpm(yDiff)[,conditionNum == 3])) / rowMeans(cpm(yDiff)[,conditionNum == 3]))

# average difference between good and bad divided by standard deviation of that batch
d1 = as.data.frame(abs(rowMeans(cpm(yDiff)[,conditionNum == 0]) - rowMeans(cpm(yDiff)[,conditionNum == 2])) / 
                     (rowMeans(cpm(yDiff)[,conditionNum == 0 | conditionNum == 2])))
d2 = as.data.frame(abs(rowMeans(cpm(yDiff)[,conditionNum == 1]) - rowMeans(cpm(yDiff)[,conditionNum == 3])) / 
                     (rowMeans(cpm(yDiff)[,conditionNum == 1 | conditionNum == 3])))

meanBCV = as.data.frame(rowMeans(cbind(BCV_b1, BCV_g1, BCV_b2, BCV_g2)))
x = as.data.frame(log2(cpm(yDiff)[order(meanBCV, decreasing = FALSE),]+1))
x = x[1:1000,]
pca = prcomp(t(x))
fac = factor(conditionDiff)
colours = brewer.pal(nlevels(fac), "Paired") 
ggplot(as.data.frame(pca$x), aes(PC1,PC2, color=conditionDiff)) + geom_point(size=5) + geom_text(hjust=-.5, vjust=-.5,aes(label=as.character(1:22)))




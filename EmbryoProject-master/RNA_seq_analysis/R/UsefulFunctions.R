# Useful R functions written by me


# function to get max fold-changes in all genes with max expression of at least 10 fpkm
# input fpkm(sigGenes), or fpkm(genes(cuff)), 10
getFoldChange <- function(fpkmTable, numSamples, foldChange = 2) {

	numGenes <- nrow(fpkmTable)/numSamples

	# initialize vectors
	# maxChangeList <- NULL
	# if we are comparing 2 samples, then save values for group 1 and group 2
	if (numSamples == 2) {
		diffGenes <- data.frame(matrix(ncol=3,nrow=numGenes))
		names(diffGenes) <- c("gene_id","fold_change","samp_with_high_val")		
	} else {
		diffGenes <- data.frame(matrix(ncol=2,nrow=numGenes))
		names(diffGenes) <- c("gene_id","fold_change")
	}

	
	
	for (i in 1:numGenes) {

		# pick out just cols with sample name and fpkm val
		currGeneTable <- fpkmTable[(numSamples*(i-1)+1):(numSamples*i),c("gene_id","sample_name","fpkm")]

		# only keep samples with fpkm > 0
		nonzeroValsTable <- currGeneTable[currGeneTable[,"fpkm"] > 0,]
		
		#print(i)
		#print(currGeneTable)
		#print(nonzeroValsTable)
		#readline("Press <Enter> to continue")

		# as long as there is at least one nonzero value which is greater than 10, calculate max fold change
		if (nrow(nonzeroValsTable) > 0) {

			# if min is < 1, just replace it with 1 to avoid crazy high fold change vals
			if (max(nonzeroValsTable[,"fpkm"]) > 10) {
				maxChange <- max(nonzeroValsTable[,"fpkm"])/max(1,min(nonzeroValsTable[,"fpkm"]))
			} else {
				maxChange <- 0
			}

			# if num samples = 2, save which sample name has the higher val
			if (numSamples == 2) {

				if (nrow(nonzeroValsTable) == 1) {
					 highSampleName <- nonzeroValsTable[1,"sample_name"]
				} else if (nrow(nonzeroValsTable) > 1) { 
					if (nonzeroValsTable[1,"fpkm"] > nonzeroValsTable[2,"fpkm"]) {
						highSampleName <- nonzeroValsTable[1,"sample_name"]
					} else {
						highSampleName <- nonzeroValsTable[2,"sample_name"]
					}
				} 
			} 

		} else {

			maxChange <- 0
			highSampleName <- "none"

		}

		# record max fold-change and gene name
		# maxChangeList[i] <- maxChange
		diffGenes[i,"gene_id"] <- currGeneTable[1,"gene_id"]
		diffGenes[i,"fold_change"] <- maxChange
		diffGenes[i,"samp_with_high_val"] <- as.character(highSampleName)

		#print(highSampleName)
		#print(diffGenes[i,])
		#readline("Press <Enter> to continue")

	}

	diffGenes <- diffGenes[(diffGenes[,"fold_change"] > foldChange),]

	return(diffGenes)
 }

# function to get correlation coefficients between fpkm and decDist for each gene
# input: data.frame sorted first by gene_id, then by sample_name, and numSamples
calcGeneExpCor <- function(dataTable, numSamples) {

	# calculate number of genes
	numGenes <- nrow(dataTable)/numSamples

	#initialize output structure
	corList <- data.frame(matrix(ncol=2,nrow=numGenes))
	names(corList) <- c("gene_id", "cor_coef")

	for (i in 1:numGenes) {

		corList[i,"gene_id"] <- dataTable[(numSamples*(i-1)+1),"gene_id"]
		corList[i,"cor_coef"] <- cor(dataTable[(numSamples*(i-1)+1):(numSamples*i),"fpkm"],dataTable[(numSamples*(i-1)+1):(numSamples*i),"normDist"])

	}

	return(corList)

}

# take 2 cols from data.frame, and make them row and col names of a new data.frame, with vals from a 3rd col inside the frame
convertDataFrameCol <- function(dataFrameIn, rowName, colName, dataName) {

	# make output data frame
	numCol <- length(unique(dataFrameIn[, c(colName)]))
	numRow <- nrow(dataFrameIn)/numCol
	dataFrameOut <- data.frame(matrix(ncol=numCol, nrow=numRow))
	names(dataFrameOut) <- dataFrameIn[1:numCol,c(colName)]

	rowNameList <- NULL

	for (i in 1:numRow) {

		dataFrameOut[i,] <- as.numeric(dataFrameIn[(numCol*(i-1)+1):(numCol*i), c(dataName)])
		rowNameList[i] <- dataFrameIn[numCol*i, c(rowName)]

	}

	row.names(dataFrameOut) <- rowNameList
	return(dataFrameOut)

}

# find which strings two arrays of strings have in common (to compare sigGenes between DESeq and edgeR output)
findCommonGenes <- function(listA, listB) {
  
  commonGenes = NULL
  
  # iterate over all rows in listA
  for (i in 1:(nrow(listA))) {
        
    currGeneName = listA[i,1]
    
    if (currGeneName %in% listB[,1]){
      commonGenes = append(commonGenes, as.character(currGeneName))
    }
    
  }
  
  return(commonGenes)
  
}

# calculate length of each gene in a GTF file for RPKM calculation
calc_length <- function(x) {
  sum(elementMetadata(x)$widths)
}


# combine two gene count tables, keeping only genes in common
combineGeneExprData <- function(countTable1, countTable2){
  
  numGenesTable1 = nrow(countTable1)
  combinedCountTable = data.frame(matrix(ncol = (ncol(countTable1) + ncol(countTable2)), nrow = 1))
  firstRowDone = 0
  names(combinedCountTable) = c(names(countTable1), names(countTable2))
  geneNamesList = NULL
  
  namesTable1 = row.names(countTable1)
  namesTable2 = row.names(countTable2)
    
  for (i in 1:numGenesTable1) {
    
    # see if current gene name in first count table is found within row names of second count table
    currGeneName = namesTable1[i]
    matchIndices = grep(currGeneName, namesTable2)
    
    # if it is found, take that row from second count table and append it to current row of 
    if (length(matchIndices) > 0){
      
      bestIndex = matchIndices[nchar(namesTable2)[matchIndices] == nchar(currGeneName)]
      
      if (length(bestIndex) > 0){
        if (firstRowDone == 0){
          combinedCountTable[1,] = c(countTable1[i,],countTable2[bestIndex,])
          firstRowDone = 1
        } else {
          combinedCountTable = rbind(combinedCountTable, c(countTable1[i,],countTable2[bestIndex,]))
        }
        
        geneNamesList = c(geneNamesList, currGeneName)
      }
      
    }

  }
  
  row.names(combinedCountTable) = geneNamesList
  return(combinedCountTable)
  
}


# PCA plot for things other than DESeq objects
makePCAplot <- function (x, condition, ntop = 500) 
{
  rv = rowVars(x)
  select = order(rv, decreasing = TRUE)[seq_len(ntop)]
  pca = prcomp(t(x[select, ]))
  fac = factor(condition)
  if (length(fac) >= 3) 
    colours = brewer.pal(nlevels(fac), "Paired")
  else colours = c("#1F78B4", "#33A02C")
  xyplot(PC2 ~ PC1, groups = fac, data = as.data.frame(pca$x), 
         pch = 16, cex = 2, aspect = "iso", col = colours, main = draw.key(key = list(rect = list(col = colours), 
                                                                                      text = list(levels(fac)), rep = FALSE)))
}



addChar <- function(vectorElement, char) paste(vectorElement, char, sep="")

# find genes with at least one GO term in desiredGoList
# outputs TRUE if gene's GO terms overlap with at least one in desiredGoList
findGenesWithGoTerms = function(currGeneGo, desiredGoList){
  if (isEmpty(intersect(currGeneGo, desiredGoList))) {
    return(FALSE)
  } else {
    return(TRUE)
  }
}


# for (go in GO.wall$category[1:20]){
#   print(GOTERM[[go]])
#   cat(" ------------------------------------------------------ \n")
# }

getEntrezID = function(geneSymbol, mapping) {
  return(mapping[[geneSymbol]])
}







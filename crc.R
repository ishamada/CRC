# Read in the raw read counts
rawCounts <- read.delim("/home/islam/Downloads/E-GEOD-50760-raw-counts (1).tsv")
head(rawCounts)

# Read in the sample mappings
sampleData <- read.delim("/home/islam/Downloads/E-GEOD-50760-experiment-design (1).tsv")
head(sampleData)

# Also save a copy for later
sampleData_v2 <- sampleData

# Convert count data to a matrix of appropriate form that DEseq2 can read
geneID <- rawCounts$Gene.ID
sampleIndex <- grepl("SRR\\d+", colnames(rawCounts))
rawCounts <- as.matrix(rawCounts[,sampleIndex])
rownames(rawCounts) <- geneID
head(rawCounts)

# Convert sample variable mappings to an appropriate form that DESeq2 can read
head(sampleData)
rownames(sampleData) <- sampleData$Run
keep <- c("Sample.Characteristic.biopsy.site.", "Sample.Characteristic.individual.")
sampleData <- sampleData[,keep]
colnames(sampleData) <- c("tissueType", "individualID")
sampleData$individualID <- factor(sampleData$individualID) # note
head(sampleData)

# Put the columns of the count data in the same order as rows names of the sample mapping, then make sure it worked
rawCounts <- rawCounts[,unique(rownames(sampleData))] # note
all(colnames(rawCounts) == rownames(sampleData))

# rename the tissue types
rename_tissues <- function(x){
  x <- switch(as.character(x), "normal"="normal-looking surrounding colonic epithelium", "primary tumor"="primary colorectal cancer",  "colorectal cancer metastatic in the liver"="metastatic colorectal cancer to the liver")
  return(x)
}
sampleData$tissueType <- unlist(lapply(sampleData$tissueType, rename_tissues))

# Order the tissue types so that it is sensible and make sure the control sample is first: normal sample -> primary tumor -> metastatic tumor
sampleData$tissueType <- factor(sampleData$tissueType, levels=c("normal-looking surrounding colonic epithelium", "primary colorectal cancer", "metastatic colorectal cancer to the liver"))

library(DESeq2)

# Create the DEseq2DataSet object
deseq2Data <- DESeqDataSetFromMatrix(countData=rawCounts, colData=sampleData, design= ~ individualID + tissueType)

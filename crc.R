# Read in the raw read counts
rawCountsOrg <- read.delim("/home/islam/Downloads/E-GEOD-50760-raw-counts (1).tsv")
head(rawCounts)

# Read in the sample mappings
sampleData <- read.delim("/home/islam/Downloads/E-GEOD-50760-experiment-design (1).tsv")
head(sampleData)


# tmp <- (which(rawCountstemp$Gene.Name != "" ))

# write.table(tmp, file='test.tsv')

# length(unique(rawCountstemp$Gene.Name))


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


dim(deseq2Data)
dim(deseq2Data[rowSums(counts(deseq2Data)) > 5, ])


# Perform pre-filtering of the data
deseq2Data <- deseq2Data[rowSums(counts(deseq2Data)) > 5, ]

BiocManager::install("BiocParallel")

library(BiocParallel)

register(MulticoreParam(4))

deseq2Data <- DESeq(deseq2Data) # , parallel = TRUE)


# Extract differential expression results
# For "tissueType" perform primary vs normal comparison
deseq2Results <- results(deseq2Data, contrast=c("tissueType", "primary colorectal cancer", "normal-looking surrounding colonic epithelium"))

# View summary of results
summary(deseq2Results)


plotMA(deseq2Results)


# Load libraries
# install.packages(c("ggplot2", "scales", "viridis"))
library(ggplot2)
library(scales) # needed for oob parameter
library(viridis)

# Coerce to a data frame
deseq2ResDF <- as.data.frame(deseq2Results)

# Examine this data frame
head(deseq2ResDF)

# Set a boolean column for significance
deseq2ResDF$significant <- ifelse(deseq2ResDF$padj < .1, "Significant", NA)

# Plot the results similar to DEseq2
ggplot(deseq2ResDF, aes(baseMean, log2FoldChange, colour=significant)) + geom_point(size=1) + scale_y_continuous(limits=c(-3, 3), oob=squish) + scale_x_log10() + geom_hline(yintercept = 0, colour="tomato1", size=2) + labs(x="mean of normalized counts", y="log fold change") + scale_colour_manual(name="q-value", values=("Significant"="red"), na.value="grey50") + theme_bw()

# Let's add some more detail
ggplot(deseq2ResDF, aes(baseMean, log2FoldChange, colour=padj)) + geom_point(size=1) + scale_y_continuous(limits=c(-3, 3), oob=squish) + scale_x_log10() + geom_hline(yintercept = 0, colour="darkorchid4", size=1, linetype="longdash") + labs(x="mean of normalized counts", y="log fold change") + scale_colour_viridis(direction=-1, trans='sqrt') + theme_bw() + geom_density_2d(colour="black", size=2)


length(which(deseq2ResDF$significant == "Significant"))


# Extract counts for the gene otop2
otop2Counts <- plotCounts(deseq2Data, gene="ENSG00000183034", intgroup=c("tissueType", "individualID"), returnData=TRUE)

# Plot the data using ggplot2
colourPallette <- c("#7145cd","#bbcfc4","#90de4a","#cd46c1","#77dd8e","#592b79","#d7c847","#6378c9","#619a3c","#d44473","#63cfb6","#dd5d36","#5db2ce","#8d3b28","#b1a4cb","#af8439","#c679c0","#4e703f","#753148","#cac88e","#352b48","#cd8d88","#463d25","#556f73")
ggplot(otop2Counts, aes(x=tissueType, y=count, colour=individualID, group=individualID)) + geom_point() + geom_line() + theme_bw() + theme(axis.text.x=element_text(angle=15, hjust=1)) + scale_colour_manual(values=colourPallette) + guides(colour=guide_legend(ncol=3)) + ggtitle("OTOP2")


deseq2ResDF["ENSG00000183034",]
rawCounts["ENSG00000183034",]
normals=row.names(sampleData[sampleData[,"tissueType"]=="normal-looking surrounding colonic epithelium",])
primaries=row.names(sampleData[sampleData[,"tissueType"]=="primary colorectal cancer",])
rawCounts["ENSG00000183034",normals]
rawCounts["ENSG00000183034",primaries]


# Transform count data using the variance stablilizing transform
deseq2VST <- vst(deseq2Data)

# Convert the DESeq transformed object to a data frame
deseq2VST <- assay(deseq2VST)
deseq2VST <- as.data.frame(deseq2VST)
deseq2VST$Gene <- rownames(deseq2VST)
head(deseq2VST)

# Keep only the significantly differentiated genes where the fold-change was at least 3
sigGenes <- rownames(deseq2ResDF[deseq2ResDF$padj <= .05 & abs(deseq2ResDF$log2FoldChange) > 3,])
deseq2VST <- deseq2VST[deseq2VST$Gene %in% sigGenes,]

which (as.matrix(rawCountsOrg$Gene.ID) == sigGenes)



rawCountsOrg["ENSG00000065675"]

rc <- as.matrix(rawCountsOrg)
rc[sigGenes,]

rownames(rc) <- (rawCountsOrg$Gene.ID)

top222 <- rc[sigGenes,2]

table(top222)

# Convert the VST counts to long format for ggplot2
library(reshape2)

# First compare wide vs long version
deseq2VST_wide <- deseq2VST
deseq2VST_long <- melt(deseq2VST, id.vars=c("Gene"))

head(deseq2VST_wide)
head(deseq2VST_long)

# Now overwrite our original data frame with the long format
deseq2VST <- melt(deseq2VST, id.vars=c("Gene"))

# Make a heatmap
heatmap <- ggplot(deseq2VST, aes(x=variable, y=Gene, fill=value)) + geom_raster() + scale_fill_viridis(trans="sqrt") + theme(axis.text.x=element_text(angle=65, hjust=1), axis.text.y=element_blank(), axis.ticks.y=element_blank())
heatmap


#============================

# Convert the significant genes back to a matrix for clustering
deseq2VSTMatrix <- dcast(deseq2VST, Gene ~ variable)
rownames(deseq2VSTMatrix) <- deseq2VSTMatrix$Gene
deseq2VSTMatrix$Gene <- NULL

# Compute a distance calculation on both dimensions of the matrix
distanceGene <- dist(deseq2VSTMatrix)
distanceSample <- dist(t(deseq2VSTMatrix))

# Cluster based on the distance calculations
clusterGene <- hclust(distanceGene, method="average")
clusterSample <- hclust(distanceSample, method="average")

# Construct a dendogram for samples
install.packages("ggdendro")
library(ggdendro)
sampleModel <- as.dendrogram(clusterSample)
sampleDendrogramData <- segment(dendro_data(sampleModel, type = "rectangle"))
sampleDendrogram <- ggplot(sampleDendrogramData) + geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + theme_dendro()

# Re-factor samples for ggplot2
deseq2VST$variable <- factor(deseq2VST$variable, levels=clusterSample$labels[clusterSample$order])

# Construct the heatmap. note that at this point we have only clustered the samples NOT the genes
heatmap <- ggplot(deseq2VST, aes(x=variable, y=Gene, fill=value)) + geom_raster() + scale_fill_viridis(trans="sqrt") + theme(axis.text.x=element_text(angle=65, hjust=1), axis.text.y=element_blank(), axis.ticks.y=element_blank())
heatmap

# Combine the dendrogram and the heatmap
install.packages("gridExtra")
library(gridExtra)
grid.arrange(sampleDendrogram, heatmap, ncol=1, heights=c(1,5))


# ====== fixing ========

# Load in libraries necessary for modifying plots
#install.packages("gtable")
library(gtable)
library(grid)

# Modify the ggplot objects
sampleDendrogram_1 <- sampleDendrogram + scale_x_continuous(expand=c(.0085, .0085)) + scale_y_continuous(expand=c(0, -10))
heatmap_1 <- heatmap + scale_x_discrete(expand=c(0, 0)) + scale_y_discrete(expand=c(-1, 1))

# Convert both grid based objects to grobs
sampleDendrogramGrob <- ggplotGrob(sampleDendrogram_1)
heatmapGrob <- ggplotGrob(heatmap_1)

# Check the widths of each grob
sampleDendrogramGrob$widths
heatmapGrob$widths

# Add in the missing columns
sampleDendrogramGrob <- gtable_add_cols(sampleDendrogramGrob, heatmapGrob$widths[7], 6)
sampleDendrogramGrob <- gtable_add_cols(sampleDendrogramGrob, heatmapGrob$widths[8], 7)

# Make sure every width between the two grobs is the same
maxWidth <- unit.pmax(sampleDendrogramGrob$widths, heatmapGrob$widths)
sampleDendrogramGrob$widths <- as.list(maxWidth)
heatmapGrob$widths <- as.list(maxWidth)

# Arrange the grobs into a plot
finalGrob <- arrangeGrob(sampleDendrogramGrob, heatmapGrob, ncol=1, heights=c(2,5))

# Draw the plot
grid.draw(finalGrob)

# ========================================

# Re-order the sample data to match the clustering we did
sampleData_v2$Run <- factor(sampleData_v2$Run, levels=clusterSample$labels[clusterSample$order])

# Construct a plot to show the clinical data
colours <- c("#743B8B", "#8B743B", "#8B3B52")
sampleClinical <- ggplot(sampleData_v2, aes(x=Run, y=1, fill=Sample.Characteristic.biopsy.site.)) + geom_tile() + scale_x_discrete(expand=c(0, 0)) + scale_y_discrete(expand=c(0, 0)) + scale_fill_manual(name="Tissue", values=colours) + theme_void()

# Convert the clinical plot to a grob
sampleClinicalGrob <- ggplotGrob(sampleClinical)

# Make sure every width between all grobs is the same
maxWidth <- unit.pmax(sampleDendrogramGrob$widths, heatmapGrob$widths, sampleClinicalGrob$widths)
sampleDendrogramGrob$widths <- as.list(maxWidth)
heatmapGrob$widths <- as.list(maxWidth)
sampleClinicalGrob$widths <- as.list(maxWidth)

# Arrange and output the final plot
finalGrob <- arrangeGrob(sampleDendrogramGrob, sampleClinicalGrob, heatmapGrob, ncol=1, heights=c(2,1,5))
grid.draw(finalGrob)

#====================================

library(gage)

library(DESeq2)

tumor_v_normal_DE <- results(deseq2Data, contrast=c("tissueType", "primary colorectal cancer", "normal-looking surrounding colonic epithelium"))



# set up KEGG database
kg.hsa <- kegg.gsets(species="hsa")


kegg.sigmet.gs <- kg.hsa$kg.sets[kg.hsa$sigmet.idx]
kegg.dise.gs <- kg.hsa$kg.sets[kg.hsa$dise.idx]

# set up GO database
go.hs <- go.gsets(species="human")
go.bp.gs <- go.hs$go.sets[go.hs$go.subs$BP]
go.mf.gs <- go.hs$go.sets[go.hs$go.subs$MF]
go.cc.gs <- go.hs$go.sets[go.hs$go.subs$CC]

go.bp.gs$`GO:0000097 sulfur amino acid biosynthetic process`[1]


# load in libraries to annotate data
library(AnnotationDbi)
library(org.Hs.eg.db)

# annotate the deseq2 results with additional gene identifiers
tumor_v_normal_DE$symbol <- mapIds(org.Hs.eg.db, keys=row.names(tumor_v_normal_DE), column="SYMBOL", keytype="ENSEMBL", multiVals="first")
tumor_v_normal_DE$entrez <- mapIds(org.Hs.eg.db, keys=row.names(tumor_v_normal_DE), column="ENTREZID", keytype="ENSEMBL", multiVals="first")
tumor_v_normal_DE$name <- mapIds(org.Hs.eg.db, keys=row.names(tumor_v_normal_DE), column="GENENAME", keytype="ENSEMBL", multiVals="first")


# grab the log fold changes for everything
tumor_v_normal_DE.fc <- tumor_v_normal_DE$log2FoldChange
names(tumor_v_normal_DE.fc) <- tumor_v_normal_DE$entrez


# Run enrichment analysis on all log fc
fc.kegg.sigmet.p <- gage(tumor_v_normal_DE.fc, gsets = kegg.sigmet.gs)
fc.kegg.dise.p <- gage(tumor_v_normal_DE.fc, gsets = kegg.dise.gs)
fc.go.bp.p <- gage(tumor_v_normal_DE.fc, gsets = go.bp.gs)
fc.go.mf.p <- gage(tumor_v_normal_DE.fc, gsets = go.mf.gs)
fc.go.cc.p <- gage(tumor_v_normal_DE.fc, gsets = go.cc.gs)


# covert the kegg results to data frames
fc.kegg.sigmet.p.up <- as.data.frame(fc.kegg.sigmet.p$greater)
fc.kegg.dise.p.up <- as.data.frame(fc.kegg.dise.p$greater)

fc.kegg.sigmet.p.down <- as.data.frame(fc.kegg.sigmet.p$less)
fc.kegg.dise.p.down <- as.data.frame(fc.kegg.dise.p$less)

# convert the go results to data frames
fc.go.bp.p.up <- as.data.frame(fc.go.bp.p$greater)
fc.go.mf.p.up <- as.data.frame(fc.go.mf.p$greater)
fc.go.cc.p.up <- as.data.frame(fc.go.cc.p$greater)

fc.go.bp.p.down <- as.data.frame(fc.go.bp.p$less)
fc.go.mf.p.down <- as.data.frame(fc.go.mf.p$less)
fc.go.cc.p.down <- as.data.frame(fc.go.cc.p$less)



dim(fc.go.bp.p.up)



# Load the AnnotationDbi package and the GO.db package
library(AnnotationDbi)
library(GO.db)

# Choose a specific GO term you are interested in, for example "GO:0042254" (cellularization)
go_id <- "GO:0042254"

# Get the gene IDs associated with the GO term using the GO.db package
genes <- select(org.Hs.eg.db, keys =  "GO:0006268" , columns = c("ENTREZID","ENSEMBL","SYMBOL"), keytype = "GO")

columns(org.Hs.eg.db)
columns(GO.db)

genes <- select(org.Hs.eg.db, keys =  "PLA2G6" , columns = c("ENTREZID","ENSEMBL","SYMBOL"), keytype = "SYMBOL")



# View the list of gene IDs
genes$ENTREZID



packageVersion("org.Hs.eg.db")
?org.Hs.eg.db


################################################################################
################ Exercise 1 ####################################################

#!!!!!!!!!! "fc.go.bp.p.up" and "go.bp.gs" should already be in your R environment, see tutorial above !!!!!!!!!#

# get genes associated with pathways
a <- function(sigPathways, geneSetDB){
  
  # grab the names of all significant pathways
  sigPathwaysID <- rownames(sigPathways)
  
  # subset the geneSet list to only those in the significant pathways
  sigPathwaysGenes <- geneSetDB[which(names(geneSetDB) %in% sigPathwaysID)]
  numberSigPathways <- length(sigPathwaysGenes)
  sigPathwaysGenes <- unlist(sigPathwaysGenes)
  
  # count the number of times a gene occurs in these significant pathways
  sigPathwaysTable <- plyr::count(sigPathwaysGenes)
  
  # annotate these final results with the gene symbol and some extra information
  sigPathwaysTable$symbol <- mapIds(org.Hs.eg.db, keys=as.character(sigPathwaysTable$x), column="SYMBOL", keytype="ENTREZID", multiVals="first")
  sigPathwaysTable$proportion <- sigPathwaysTable$freq/numberSigPathways
  
  # sort and return the results
  sigPathwaysTable <- sigPathwaysTable[order(-sigPathwaysTable$proportion),]
  
  return(sigPathwaysTable)
}

# run the function defined above
geneTable1 <- a(na.omit(fc.go.bp.p.up[fc.go.bp.p.up$q.val <= .05,]), go.bp.gs)
geneTable2 <- a(na.omit(fc.kegg.sigmet.p.up[fc.kegg.sigmet.p.up$q.val <= .05,]), kegg.sigmet.gs)
geneTable3 <- a((fc.kegg.dise.p.up[fc.kegg.dise.p.up,]), kegg.dise.gs)
geneTable4 <- a(na.omit(fc.go.mf.p.up[fc.go.mf.p.up$q.val <= .05,]), go.mf.gs)
geneTable5 <- a(na.omit(fc.go.cc.p.up[fc.go.cc.p.up$q.val <= .05,]), go.cc.gs)

#========== pathway viz ====

library(pathview)


# View the hsa03430 pathway from the pathway analysis
fc.kegg.sigmet.p.up[grepl("hsa03430", rownames(fc.kegg.sigmet.p.up), fixed=TRUE),]

# Overlay the expression data onto this pathway
pathview(gene.data=tumor_v_normal_DE.fc, species="hsa", pathway.id="hsa03430")


# View the hsa03430 pathway from the pathway analysis
fc.kegg.sigmet.p.up[grepl("hsa03460", rownames(fc.kegg.sigmet.p.up), fixed=TRUE),]

# Overlay the expression data onto this pathway
pathview(gene.data=tumor_v_normal_DE.fc, species="hsa", pathway.id="hsa03460", kegg.native=FALSE)







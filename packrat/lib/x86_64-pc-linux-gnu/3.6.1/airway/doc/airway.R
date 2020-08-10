## ------------------------------------------------------------------------
suppressPackageStartupMessages( library( "GEOquery" ) )
suppressPackageStartupMessages( library( "airway" ) )
dir <- system.file("extdata",package="airway")
geofile <- file.path(dir, "GSE52778_series_matrix.txt")
gse <- getGEO(filename=geofile)
pdata <- pData(gse)[,grepl("ch1",names(pData(gse)))]
names(pdata) <- c("treatment","tissue","ercc_mix","cell","celltype")
pdataclean <- data.frame(treatment=sub("treatment: (.*)","\\1",pdata$treatment),
                         cell=sub("cell line: (.*)","\\1",pdata$cell),
                         row.names=rownames(pdata))
pdataclean$dex <- ifelse(grepl("Dex",pdataclean$treatment),"trt","untrt")
pdataclean$albut <- ifelse(grepl("Albut",pdataclean$treatment),"trt","untrt")
pdataclean$SampleName <- rownames(pdataclean)
pdataclean$treatment <- NULL

## ------------------------------------------------------------------------
srafile <- file.path(dir, "SraRunInfo_SRP033351.csv")
srp <- read.csv(srafile)
srpsmall <- srp[,c("Run","avgLength","Experiment","Sample","BioSample","SampleName")]

## ------------------------------------------------------------------------
coldata <- merge(pdataclean, srpsmall, by="SampleName")
rownames(coldata) <- coldata$Run
coldata <- coldata[coldata$albut == "untrt",]
coldata$albut <- NULL
coldata

## ----eval=FALSE----------------------------------------------------------
#  write.csv(coldata, file="sample_table.csv")

## ----eval=FALSE----------------------------------------------------------
#  library( "GenomicFeatures" )
#  txdb <- makeTranscriptDbFromBiomart( biomart="ensembl", dataset="hsapiens_gene_ensembl")
#  exonsByGene <- exonsBy( txdb, by="gene" )

## ----eval=FALSE----------------------------------------------------------
#  sampleTable <- read.csv( "sample_table.csv", row.names=1 )
#  fls <- file.path("aligned",rownames(sampleTable), ".bam")
#  library( "Rsamtools" )
#  bamLst <- BamFileList( fls, yieldSize=2000000 )

## ----eval=FALSE----------------------------------------------------------
#  library( "BiocParallel" )
#  register( MulticoreParam( workers=8 ) )
#  library( "GenomicAlignments" )
#  airway <- summarizeOverlaps( features=exonsByGene, reads=bamLst,
#                              mode="Union", singleEnd=FALSE,
#                              ignore.strand=TRUE, fragments=TRUE )

## ----eval=FALSE----------------------------------------------------------
#  colData(airway) <- DataFrame( sampleTable )

## ----eval=FALSE----------------------------------------------------------
#  library( "annotate" )
#  miame <- list(pmid2MIAME("24926665"))
#  miame[[1]]@url <- "http://www.ncbi.nlm.nih.gov/pubmed/24926665"
#  # because R's CHECK doesn't like non-ASCII characters in data objects
#  # or in vignettes. the actual char was used in the first argument
#  miame[[1]]@abstract <- gsub("micro","micro",abstract(miame[[1]]))
#  miame[[1]]@abstract <- gsub("beta","beta",abstract(miame[[1]]))
#  metadata(airway) <- miame
#  save(airway, file="airway.RData")

## ------------------------------------------------------------------------
library("airway")
data(airway)
airway
as.data.frame(colData(airway))
summary(colSums(assay(airway))/1e6)
metadata(rowRanges(airway))

## ------------------------------------------------------------------------
sessionInfo()


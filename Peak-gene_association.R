#############################################################################
####################                                #########################
####################    PEAK LOCATION ANALYSIS      #########################
####################                                #########################
#############################################################################

##-----------------------------------------------------------------------------
## TRANSCRIPT DATABASE CROMOHALOBACTER SALEXIGENS 
##-----------------------------------------------------------------------------

# BiocManager::install("GenomicFeatures")

library(GenomicFeatures)

# Create a TxDB object using the GTF file with genomic annotations
# obtained from Ensembl Bacteria

txDB_Chb_Sxg <- makeTxDbFromGFF(file = "../Chromohalobacter_salexigens_dsm_3043_gca_000055785.ASM5578v1.56.gtf",
                                format = "gtf")
saveDb(txDB_Chb_Sxg, file = "txDB_Chb_Sxg.sqlite")

# Save the names for all the genes in the genome and write a table with them

GENE_NAMES_LIST <- genes(txDB_Chb_Sxg)
GENE_NAMES_LIST <- GENE_NAMES_LIST@ranges@NAMES

write.table(GENE_NAMES_LIST, "LISTAS DE GENES/CHROMOHALOBACTER_ALL_GENES.txt",
            quote = FALSE, row.names = FALSE)


##-----------------------------------------------------------------------------
## PROMOTER ANNOTATION
##-----------------------------------------------------------------------------

library(ChIPseeker)
library(GenomicRanges)

# Loading the database just created

DB_Chromohalobacter <- loadDb("txDB_Chb_Sxg.sqlite")

# Get promoters from the database

promoters <- getPromoters(TxDb = DB_Chromohalobacter, upstream = 1000, downstream = 1000)

# Read file with Chromohalobacter Salexigens genome annotations

GOCHSAL<- read.csv("../GENOMA CHROMOHALOBACTER 2019 -newKO.csv", sep = ";")

##------------------------------------------------------------------------------
## PEAK ANALYSIS FUR1 0.6
##------------------------------------------------------------------------------

# Read peak files

FUR1_AP_0.6 <- readPeakFile("../IGV/FUR1/0.6/FUR1_Peaks_0.6_peaks.narrowPeak")

# Annotate peaks

FUR1_AP_0.6_PeakAnno <- annotatePeak(peak = FUR1_AP_0.6, 
                                     tssRegion=c(-1000, 1000),
                                     TxDb=DB_Chromohalobacter)

# Visualización con distintas gráficas

par(mfrow = c(1,5))
covplot(FUR1_AP_0.6, weightCol="V5")
plotDistToTSS(FUR1_AP_0.6_PeakAnno)
plotAnnoPie(FUR1_AP_0.6_PeakAnno, ti)
plotAnnoBar(FUR1_AP_0.6_PeakAnno)
upsetplot(FUR1_AP_0.6_PeakAnno)
plotPeakProf2(peak = FUR1_AP_0.6, upstream = rel(0.2), downstream = rel(0.2),
              conf = 0.95, by = "gene", type = "body", nbin = 800,
              TxDb = DB_Chromohalobacter, weightCol = "V5",ignore_strand = F)

FUR1_AP_0.6_PeakAnno <- as.data.frame(FUR1_AP_0.6_PeakAnno)

# Filtering, retaining upstream gene regulators

filteredGenes_FUR1_AP_0.6 <- subset(FUR1_AP_0.6_PeakAnno,
                                    upstream = TRUE, downstream = FALSE, regulon = TRUE)
filteredGenes_FUR1_AP_0.6 <- as.data.frame(filteredGenes_FUR1_AP_0.6)

# Filter peaks assigned to non-promoter regions and write a table with 
# associated genes for functional enrichment

target_genes_FUR1_AP_0.6 <- as.data.frame(filteredGenes_FUR1_AP_0.6$geneId[filteredGenes_FUR1_AP_0.6$annotation == "Promoter"])
colnames(target_genes_FUR1_AP_0.6) <- "target_genes_FUR1_AP_0.6"
write.table(target_genes_FUR1_AP_0.6, file = "LISTAS DE GENES/genes_FUR1_0.6.txt",
            quote = FALSE, row.names = FALSE)

##------------------------------------------------------------------------------
## PEAK ANALYSIS FUR 1 2.5
##------------------------------------------------------------------------------

# Read peak files

FUR1_AP_2.5 <- readPeakFile("../IGV/FUR1/2.5/FUR1_Peaks_2.5_peaks.narrowPeak")

# Annotate peaks

FUR1_AP_2.5_PeakAnno <- annotatePeak(peak = FUR1_AP_2.5, 
                                     tssRegion=c(-1000, 1000),
                                     TxDb=DB_Chromohalobacter)

# Data visualization

covplot(FUR1_AP_2.5, weightCol="V5")
par(mfrow = c(1,2))
plotDistToTSS(FUR1_AP_2.5_PeakAnno)
plotAnnoPie(FUR1_AP_2.5_PeakAnno)
plotAnnoBar(FUR1_AP_2.5_PeakAnno)
upsetplot(FUR1_AP_2.5_PeakAnno)
plotPeakProf2(peak = FUR1_AP_2.5, upstream = rel(0.2), downstream = rel(0.2),
              conf = 0.95, by = "gene", type = "body", nbin = 800,
              TxDb = DB_Chromohalobacter, weightCol = "V5",ignore_strand = F)

FUR1_AP_2.5_PeakAnno <- as.data.frame(FUR1_AP_2.5_PeakAnno)

# Filtering, retaining upstream gene regulators

filteredGenes_FUR1_AP_2.5 <- subset(FUR1_AP_2.5_PeakAnno, upstream = TRUE,
                                    downstream = FALSE, regulon = TRUE)
filteredGenes_FUR1_AP_2.5 <- as.data.frame(filteredGenes_FUR1_AP_2.5)

# Filter peaks assigned to non-promoter regions and write a table with 
# associated genes for gene-operon association 

target_genes_FUR1_AP_2.5 <- as.data.frame(filteredGenes_FUR1_AP_2.5$geneId[filteredGenes_FUR1_AP_2.5$annotation == "Promoter"])
colnames(target_genes_FUR1_AP_2.5) <- "target_genes_FUR1_AP_2.5"
write.table(target_genes_FUR1_AP_2.5, file = "LISTAS DE GENES/genes_FUR1_2.5.txt",
            quote = FALSE, row.names = FALSE, sep = "\t")

##------------------------------------------------------------------------------
## PEAK ANALYSIS FUR 2 AP 0.6
##------------------------------------------------------------------------------

# Read the peak file obtained from the macs2 analysis

AP_0.6 <- readPeakFile("../IGV/FUR2/FUR2 0.6/peaks_FUR2_0.6_peaks.narrowPeak",
                       header = FALSE)

# Get annotations from the peaks from the transcript database

AP_0.6_PeakAnno <- annotatePeak(peak = AP_0.6, 
                                tssRegion=c(-1000, 1000),
                                TxDb=DB_Chromohalobacter)

# Data visualization

par(mfrow = c(1,2))
covplot(AP_0.6, weightCol="V5")
plotDistToTSS(AP_0.6_PeakAnno)
plotAnnoPie(AP_0.6_PeakAnno)
plotAnnoBar(AP_0.6_PeakAnno)
upsetplot(AP_0.6_PeakAnno)
plotPeakProf2(peak = AP_0.6, upstream = rel(0.2), downstream = rel(0.2),
              conf = 0.95, by = "gene", type = "body", nbin = 800,
              TxDb = DB_Chromohalobacter, weightCol = "V5",ignore_strand = F)

AP_0.6_PeakAnnotation1 <- as.data.frame(AP_0.6_PeakAnno)

# Filtering, retaining upstream genes and the gene regulators

filteredGenes <- subset(AP_0.6_PeakAnno, upstream = TRUE, downstream = FALSE, regulon = TRUE)
filteredGenes0.6_df <- as.data.frame(filteredGenes)

# Filter peaks assigned to non-promoter regions and write a table with 
# associated genes for gene-operon association

targetGenes0.6 <- as.data.frame(filteredGenes0.6_df$geneId[filteredGenes0.6_df$annotation == "Promoter"])
colnames(targetGenes0.6) <- "target_genes_FUR2_AP_0.6"
write.table(targetGenes0.6, file = "LISTAS DE GENES/genes_FUR2_0.6.txt",
            quote = FALSE, row.names = FALSE, sep = "\t")

##------------------------------------------------------------------------------
## PEAK ANALYSIS FUR 2 AP 2.5 
##------------------------------------------------------------------------------

# Read peak files

AP_2.5 <- readPeakFile("../IGV/FUR2/FUR2 2.5/AP_peaks_2.5_peaks.narrowPeak",
                           header = FALSE)

# Peak annotation 

AP_2.5_PeakAnno <- annotatePeak(peak = AP_2.5, 
                                    tssRegion=c(-1000, 1000),
                                    TxDb=DB_Chromohalobacter)

# Visualización con distintas gráficas

par(mfrow = c(1,2))
covplot(AP_2.5, weightCol="V5")
plotDistToTSS(AP_2.5_PeakAnno)
plotAnnoPie(AP_2.5_PeakAnno)
plotAnnoBar(AP_2.5_PeakAnno)
upsetplot(AP_2.5_PeakAnno)
plotPeakProf2(peak = AP_2.5, upstream = rel(0.2), downstream = rel(0.2),
              conf = 0.95, by = "gene", type = "body", nbin = 800,
              TxDb = DB_Chromohalobacter, weightCol = "V5",ignore_strand = F)

AP_2.5_PeakAnnotation1 <- as.data.frame(AP_2.5_PeakAnno)

# Filtering, retaining upstream genes and the gene regulators

filteredGenes_2.5 <- subset(AP_2.5_PeakAnno, upstream = TRUE, downstream = FALSE, regulon = TRUE)
filteredGenes2.5_df <- as.data.frame(filteredGenes_2.5)

# Filter peaks assigned to non-promoter regions and write a table with 
# associated genes for gene-operon association

target_genes_2.5 <- as.data.frame(filteredGenes2.5_df$geneId[filteredGenes2.5_df$annotation == "Promoter"])
colnames(target_genes_2.5) <- "target_genes_FUR2_AP_2.5"
write.table(target_genes_2.5, file = "LISTAS DE GENES/genes_FUR2_2.5.txt",
            quote = FALSE, row.names = FALSE, sep = "\t")

##------------------------------------------------------------------------------
## Generate images
##------------------------------------------------------------------------------

library(VennDiagram)

# Venn diagram Fur 1 0.6 vs 2.5

venn.data.Fur1 <- list("Fur1 0.6M NaCl" = filteredGenes_FUR1_AP_0.6$geneId, 
                   "Fur1 2.5M NaCl" = filteredGenes_FUR1_AP_2.5$geneId)

venn.diagram(
  venn.data.Fur1,
  filename = "Fur1.png",
  output = TRUE,
  imagetype = "png",
  col = "transparent",
  fill = c("lightgreen", "red"),  
  alpha = 0.5, 
  label.col = "black", 
  cex = 1,
  cat.pos = 0,
  fontfamily = "sans",
  fontface = "bold",
  height = 1500,
  width = 3000
)

# Venn diagram Fur 2 0.6M NaCl vs 2.5M NaCl

venn.data.Fur2 <- list("Fur2 0.6M NaCl" = filteredGenes0.6_df$geneId, 
                       "Fur2 2.5M NaCl" = filteredGenes2.5_df$geneId)

venn.diagram(
  venn.data.Fur2,
  filename = "Fur2.png",
  output = TRUE,
  imagetype = "png",
  col = "transparent",
  fill = c("lightgreen", "red"),  
  alpha = 0.5, 
  label.col = "black", 
  cat.cex = 1,
  cat.pos = 0,
  fontfamily = "sans",
  fontface = "bold",
  height = 1500,
  width = 3000,
)

## Peak location relative to TSS for all conditions

peakList <- list("Fur1 0.6 NaCl" = FUR1_AP_0.6_PeakAnno,
                 "Fur1 2.5 NaCl" = FUR1_AP_2.5_PeakAnno,
                 "Fur2 0.6 NaCl" = AP_0.6_PeakAnno,
                 "Fur2 2.5 NaCl" = AP_2.5_PeakAnno)

plotDistToTSS(peakList, title = "Peak location relative to TSS")




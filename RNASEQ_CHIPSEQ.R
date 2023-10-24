############################################################################
####################                             ###########################
####################       Fur1 and Fur2         ###########################
####################   gene regulation profile   ###########################
####################                             ###########################
############################################################################


## Read files with RNA_Seq results

# Fur1 0.6M NaCl

RNASEQ_FUR1_0.6 <- read.csv(file = "Fur1_0_6_vs_Wt_06_results.csv")

# Fur1 2.5M NaCl

RNASEQ_FUR1_2.5 <- read.csv(file = "Fur1_2_5_vs_Wt_25_results.csv")

# Fur2 0.6M NaCl

RNASEQ_FUR2_0.6 <- read.csv(file = "Fur2_0_6_vs_Wt_06_results.csv")

# Fur2 2.5M NaCl

RNASEQ_FUR2_2.5 <- read.csv(file = "Fur2_2_5_vs_Wt_25_results.csv")

## Read files with ChIP-Seq results 
 
# Fur1 0.6M NaCl

genes_CHIP_FUR1_0.6 <- read.csv(file = "../OPERONES/Genes Operones/genes_operone_FUR1_0.6.txt",
                                header = TRUE)

# Fur1 2.5M NaCl

genes_CHIP_FUR1_2.5 <- read.csv(file = "../OPERONES/Genes Operones/genes_operone_FUR1_2.5.txt",
                                header = TRUE)

# Fur2 0.6M NaCl

genes_CHIP_FUR2_0.6 <- read.csv(file = "../OPERONES/Genes Operones/genes_operone_FUR2_0.6.txt",
                                header = TRUE)

# Fur2 2.5M NaCl

genes_CHIP_FUR2_2.5 <- read.csv(file = "../OPERONES/Genes Operones/genes_operone_FUR2_2.5.txt",
                                header = TRUE)

## Find coincidences between datasets and write tables with the results of the 
## intersected datasets

# Fur1 0.6M NaCl

coincidences_FUR1_0.6 <- subset(RNASEQ_FUR1_0.6, `name` %in% genes_CHIP_FUR1_0.6$x)
write.table(coincidences_FUR1_0.6, file = "Coincidencias datos tratados/coincidences_FUR1_0.6.csv", row.names = FALSE, sep = "\t")

# Fur1 2.5M NaCl

coincidences_FUR1_2.5 <- subset(RNASEQ_FUR1_2.5, `name` %in% genes_CHIP_FUR1_2.5$x)
write.table(coincidences_FUR1_2.5, file = "Coincidencias datos tratados/coincidences_FUR1_2.5.csv", row.names = FALSE, sep = "\t")

# Fur2 0.6M NaCl

coincidences_FUR2_0.6 <- subset(RNASEQ_FUR2_0.6, `name` %in% genes_CHIP_FUR2_0.6$x)
write.table(coincidences_FUR2_0.6, file = "Coincidencias datos tratados/coincidences_FUR2_0.6.csv", row.names = FALSE, sep = "\t")

# Fur2 2.5M NaCl

coincidences_FUR2_2.5 <- subset(RNASEQ_FUR2_2.5, `name` %in% genes_CHIP_FUR2_2.5$x)
write.table(coincidences_FUR2_2.5, file = "Coincidencias datos tratados/coincidences_FUR2_2.5.csv", row.names = FALSE, sep = "\t")

## Create barplots comparing upregulated and downregulated genes for each set 

library(ggplot2)

## Fur1 0.6M NaCl

# Count the number of upregulated and downregulated genes

upregulated_Fur1_0.6 <- sum(coincidences_FUR1_0.6$FoldChange > 0)
downregulated_Fur1_0.6 <- sum(coincidences_FUR1_0.6$FoldChange <0)

# Create a vector with both values

diff_expr_Fur1_0.6 <- c(upregulated_Fur1_0.6, downregulated_Fur1_0.6)

## Fur1 2.5M NaCl

upregulated_Fur1_2.5 <- sum(coincidences_FUR1_2.5$FoldChange > 0)
downregulated_Fur1_2.5 <- sum(coincidences_FUR1_2.5$FoldChange <0)

diff_expr_Fur1_2.5 <- c(upregulated_Fur1_2.5, downregulated_Fur1_2.5)

## Fur2 0.6M NaCl

upregulated_Fur2_0.6 <- sum(coincidences_FUR2_0.6$FoldChange > 0)
downregulated_Fur2_0.6 <- sum(coincidences_FUR2_0.6$FoldChange <0)

diff_expr_Fur2_0.6 <- c(upregulated_Fur2_0.6, downregulated_Fur2_0.6)

## Fur2 2.5M NaCl

upregulated_Fur2_2.5 <- sum(coincidences_FUR2_2.5$FoldChange > 0)
downregulated_Fur2_2.5 <- sum(coincidences_FUR2_2.5$FoldChange <0)

diff_expr_Fur2_2.5 <- c(upregulated_Fur2_2.5, downregulated_Fur2_2.5)

diff_expr_Fur1and2_total <- c(diff_expr_Fur1_0.6,
                              diff_expr_Fur2_0.6,
                              diff_expr_Fur1_2.5,
                              diff_expr_Fur2_2.5)
## Create barplots ##

par(mfrow = c(1,4))

## Fur1 0.65M NaCl

barplot(diff_expr_Fur1_0.6, names.arg = c("Upregulated", "Downregulated"),
        col = c("purple","lightgreen" ),
        main = "Fur1 0.6M NaCl regulated genes",
        cex.names = 1.2)

## Fur2 0.6M NaCl

barplot(diff_expr_Fur2_0.6, names.arg = c("Upregulated", "Downregulated"),
        col = c("purple","lightgreen" ),
        main = "Fur2 0.6M NaCl regulated genes",
        cex.names = 1.2)

## Fur1 2.5M NaCl

barplot(diff_expr_Fur1_2.5, names.arg = c("Upregulated", "Downregulated"),
        col = c("purple","lightgreen" ),
        main = "Fur1 2.5M NaCl regulated genes",
        cex.names = 1.2)

## Fur2 2.5M NaCl

barplot(diff_expr_Fur2_2.5, names.arg = c("Upregulated", "Downregulated"),
        col = c("purple","lightgreen" ),
        main = "Fur2 2.5M NaCl regulated genes",
        cex.names = 1.2)

dev.off()

# Complete view of the differential expression for Fur1 and Fur2 0.6 and 2.5 M NaCl

barplot(diff_expr_Fur1and2_total, names.arg = c("Upregulated",
                                                "Downregulated",
                                                "Upregulated",
                                                "Downregulated",
                                                "Upregulated",
                                                "Downregulated",
                                                "Upregulated",
                                                "Downregulated"),
        col = c("purple","lightgreen" ),
        ylim = c(0,250),
        cex.names = 1)



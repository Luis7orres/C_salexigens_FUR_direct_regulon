#############################################################################
####################                                #########################
####################    GENE-OPERON ASSOCIATION     #########################
####################                                #########################
#############################################################################

# Read file with predicted operons

chs_predicted_operons <- read.csv(
  "LIST-OF-OPERONS-OF-Chromohalobacter-salexigens-DSM-3043-ProOpDB.csv")

colnames(chs_predicted_operons) <- chs_predicted_operons[1,]
chs_predicted_operons <- chs_predicted_operons[-1,]

### FUR 1

## 0.6

# Read the files containing the genes associated with peaks 

FUR1_0.6_genes <- read.table("../TxDB/LISTAS DE GENES/genes_FUR1_0.6.txt",
                             header = TRUE)

# Create a subset using the genes from this file to find if they are present
# in an operon 

operon_genes_FUR1_0.6 <- subset(x = chs_predicted_operons,
                                `Locus Name` %in% FUR1_0.6_genes$target_genes_FUR_AP_0.6)
operon_genes_FUR1_0.6 <- operon_genes_FUR1_0.6$Operon
operon_genes_FUR1_0.6 <- subset(chs_predicted_operons, Operon %in% operon_genes_FUR1_0.6)

# Save the gene IDs for functional enrichment

genes_operone_FUR1_0.6<- operon_genes_FUR1_0.6$`Locus Name`

# Write table with results

write.table(operon_genes_FUR1_0.6, "Operon_genes_FUR1_0.6.csv", sep = ";", quote = TRUE, row.names = FALSE)
write.table(operon_genes_FUR1_0.6$`Locus Name`, "Genes Operones/genes_operone_FUR1_0.6.txt", 
            quote = FALSE, row.names = FALSE)

# 2.5

# Read the files containing the genes associated with peaks 

FUR1_2.5_genes <- read.table("../TxDB/LISTAS DE GENES/genes_FUR1_2.5.txt",
                             header = TRUE)

# Create a subset using the genes from this file to find if they are present
# in an operon 

operon_genes_FUR1_2.5 <- subset(x = chs_predicted_operons,
                                `Locus Name` %in% FUR1_2.5_genes$target_genes_FUR_AP_2.5)
operon_genes_FUR1_2.5 <- operon_genes_FUR1_2.5$Operon
operon_genes_FUR1_2.5 <- subset(chs_predicted_operons, Operon %in% operon_genes_FUR1_2.5)

# Save the gene IDs for functional enrichment

genes_operone_FUR1_2.5<- operon_genes_FUR1_2.5$`Locus Name`

# Write table with results

write.table(operon_genes_FUR1_2.5, "Operon_genes_FUR1_2.5.csv", sep = ";", quote = TRUE, row.names = FALSE)
write.table(operon_genes_FUR1_2.5$`Locus Name`, "Genes Operones/genes_operone_FUR1_2.5.txt", 
            quote = FALSE, row.names = FALSE)

### FUR2

## 0.6 

# Read the files containing the genes associated with peaks 

FUR2_0.6_genes <- read.table("../TxDB/LISTAS DE GENES/genes_FUR2_0.6.txt",
                             header = FALSE)

# Create a subset using the genes from this file to find if they are present
# in an operon 

operon_genes_FUR2_0.6 <- subset(x = chs_predicted_operons,
                                `Locus Name` %in% FUR2_0.6_genes$V1)
operon_genes_FUR2_0.6 <- operon_genes_FUR2_0.6$Operon
operon_genes_FUR2_0.6 <- subset(chs_predicted_operons, Operon %in% operon_genes_FUR2_0.6)

# Save the gene IDs for functional enrichment

genes_operone_FUR2_0.6 <- operon_genes_FUR2_0.6$`Locus Name`

# Write tables with results

write.table(operon_genes_FUR2_0.6, "Operon_genes_FUR2_0.6.csv", sep = ";", quote = TRUE, row.names = FALSE)
write.table(operon_genes_FUR2_0.6$`Locus Name`, "Genes Operones/genes_operone_FUR2_0.6.txt", 
            quote = FALSE, row.names = FALSE)

## 2.5

# Read the files containing the genes associated with peaks 

FUR2_2.5_genes <- read.table("../TxDB/LISTAS DE GENES/genes_FUR2_2.5.txt",
                             header = FALSE)
# Create a subset using the genes from this file to find if they are present
# in an operon 

operon_genes_FUR2_2.5 <- subset(x = chs_predicted_operons,
                                `Locus Name` %in% FUR2_2.5_genes$V1)
operon_genes_FUR2_2.5 <- operon_genes_FUR2_2.5$Operon
operon_genes_FUR2_2.5 <- subset(chs_predicted_operons, Operon %in% operon_genes_FUR2_2.5)

# Save the gene IDs for functional enrichment

genes_operone_FUR2_2.5 <- operon_genes_FUR2_2.5$`Locus Name`

# Write tables with results

write.table(operon_genes_FUR2_2.5, "Operon_genes_FUR2_2.5.csv", sep = ";", quote = TRUE, row.names = FALSE)
write.table(operon_genes_FUR2_2.5$`Locus Name`, "Genes Operones/genes_operone_FUR2_2.5.txt", 
            quote = FALSE, row.names = FALSE)

##-----------------------------------------------------------------------------
## Functional enrichment
##-----------------------------------------------------------------------------

# Read file containing gene anotations

genome_chbsxg <- read.csv("../GENOMA CHROMOHALOBACTER 2019 -newKO.csv", sep = ";")

## FUR 1

# 0.6

# Find gene anotations for the operon genes

operon_enrichment_Fur1_0.6 <- subset(x= genome_chbsxg, `Old.locus_tag` %in% genes_operone_FUR1_0.6)

# Write a table with the results of the intersection

write.table(operon_enrichment_Fur1_0.6, "Anotaciones operones/operon_enrichment_Fur1_0.6.csv",
            sep = "\t", quote = TRUE, row.names = FALSE)

# 2.5

# Find gene anotations for the operon genes

operon_enrichment_Fur1_2.5 <- subset(x= genome_chbsxg, `Old.locus_tag` %in% genes_operone_FUR1_2.5)

# Write a table with the results of the intersection

write.table(operon_enrichment_Fur1_2.5, "Anotaciones operones/operon_enrichment_Fur1_2.5.csv",
            sep = "\t", quote = TRUE, row.names = FALSE)

## FUR2

# 0.6

# Find gene anotations for the operon genes

operon_enrichment_Fur2_0.6 <- subset(x= genome_chbsxg, `Old.locus_tag` %in% genes_operone_FUR2_0.6)

# Write a table with the results of the intersection

write.table(operon_enrichment_Fur2_0.6, "Anotaciones operones/operon_enrichment_Fur2_0.6.csv",
            sep = "\t", quote = TRUE, row.names = FALSE)

# 2.5

# Find gene anotations for the operon genes

operon_enrichment_Fur2_2.5 <- subset(x= genome_chbsxg, `Old.locus_tag` %in% genes_operone_FUR2_2.5)

# Write a table with the results of the intersection

write.table(operon_enrichment_Fur2_2.5, "Anotaciones operones/operon_enrichment_Fur2_2.5.csv",
            sep = "\t", quote = TRUE, row.names = FALSE)


##-----------------------------------------------------------------------------
## OVERLAPS IN REGULATED GENES
##-----------------------------------------------------------------------------

library(ggVennDiagram)
library(ggplot2)
library(gridExtra)

## Create lists with the sets of genes to be compared

# All datasets

operon_genes_list <- list("Fur1 0.6M NaCl" = genes_operone_FUR1_0.6,
                          "Fur1 2.5M NaCl" = genes_operone_FUR1_2.5,
                          "Fur2 0.6M NaCl" = genes_operone_FUR2_0.6,
                          "Fur2 2.5M NaCl" = genes_operone_FUR2_2.5)

# Fur1 vs Fur2 0.6M NaCl

operon_genes_list2 <-list("Fur1 0.6M NaCl" = genes_operone_FUR1_0.6,
                          "Fur2 0.6M NaCl" = genes_operone_FUR2_0.6)

# Fur1 vs Fur2 2.5M NaCl

operon_genes_list3 <- list("Fur1 2.5M NaCl" = genes_operone_FUR1_2.5,
                           "Fur2 2.5M NaCl" = genes_operone_FUR2_2.5)

# Fur2 0.6M NaCl vs 2.5M NaCl 

operon_genes_list4 <- list("Fur2 0.6M NaCl" = genes_operone_FUR2_0.6,
                           "Fur2 2.5M NaCl" = genes_operone_FUR2_2.5)

# Fur1 0.6M NaCl vs 2.5M NaCl

operon_genes_list5 <- list("Fur1 0.6M NaCl" = genes_operone_FUR1_0.6,
                           "Fur1 2.5M NaCl" = genes_operone_FUR1_2.5)

## Annotate the genes common in the four datasets and write txt files with results

# All datasets

common_genes1 <- as.data.frame(Reduce(intersect, operon_genes_list))
colnames(common_genes1) <- "common genes Fur1 and Fur2 all salinity conditions"
write.table(common_genes1, file = "Fur1and2coincidencesallconditions.txt",
            row.names = FALSE, quote = FALSE, sep = "\t")

# Fur1 vs Fur2 0.6M NaCl

common_genes2 <- as.data.frame(Reduce(intersect, operon_genes_list2))
colnames(common_genes2) <- "common genes Fur1 and Fur2 0.6M NaCl"
write.table(common_genes2, file = "Fur1and2coincidences0.6NaCl.txt",
            row.names = FALSE, quote = FALSE, sep = "\t")

# Fur1 vs Fur2 2.5M NaCl

common_genes3 <- as.data.frame(Reduce(intersect, operon_genes_list3))
colnames(common_genes3) <- "common genes Fur1 and Fur2 2.5M NaCl"
write.table(common_genes3, file = "Fur1and2coincidences2.5NaCl.txt",
            row.names = FALSE, quote = FALSE, sep = "\t")

# Fur2 0.6M NaCl vs 2.5M NaCl 

common_genes4 <- as.data.frame(Reduce(intersect, operon_genes_list4))
colnames(common_genes4) <- "common genes Fur2 0.6M and 2.5M NaCl"
write.table(common_genes4, file = "Fur2coincidences.txt",
            row.names = FALSE, quote = FALSE, sep = "\t")

# Fur1 0.6M NaCl vs 2.5M NaCl

common_genes5 <- as.data.frame(Reduce(intersect, operon_genes_list5))
colnames(common_genes5) <- "common genes Fur1 0.6M and 2.5M NaCl"
write.table(common_genes5, file = "Fur1coincidences.txt",
            row.names = FALSE, quote = FALSE, sep = "\t")

## Plot Venn diagram for visualizing coincidences

# All datasets

ggVennDiagram(operon_genes_list, label = "count", set_size = 3) + 
  scale_fill_gradient(low = "white", high = "green")

# Fur1 vs Fur2 0.6M NaCl

ggVennDiagram(operon_genes_list2, label = "count") +
  scale_fill_gradient(low = "white", high = "green")

# Fur1 vs Fur2 2.5M NaCl

ggVennDiagram(operon_genes_list3, label = "count") +
  scale_fill_gradient(low = "white", high = "green")

# Fur2 0.6M NaCl vs 2.5M NaCl 

ggVennDiagram(operon_genes_list4, label = "count") +
  scale_fill_gradient(low = "white", high = "green")

# Fur1 0.6M NaCl vs 2.5M NaCl

ggVennDiagram(operon_genes_list5, label = "count") +
  scale_fill_gradient(low = "white", high = "green")


## Assign functional annotations for genes that coincide for each TF and condition

# All datasets

common_genes1_annotation <- subset(x= genome_chbsxg, `Old.locus_tag` %in%
                                     common_genes1$`common genes`)

# Fur1 vs Fur2 0.6M NaCl

common_genes2_annotation <- subset(x= genome_chbsxg, `Old.locus_tag` %in%
                                     common_genes2$`common genes`)

# Fur1 vs Fur2 2.5M NaCl

common_genes3_annotation <- subset(x= genome_chbsxg, `Old.locus_tag` %in%
                                     common_genes3$`common genes`)

# Fur2 0.6M NaCl vs 2.5M NaCl 

common_genes4_annotation <- subset(x= genome_chbsxg, `Old.locus_tag` %in%
                                     common_genes4$`common genes`)

# Fur1 0.6M NaCl vs 2.5M NaCl

common_genes5_annotation <- subset(x= genome_chbsxg, `Old.locus_tag` %in%
                                     common_genes5$`common genes`)

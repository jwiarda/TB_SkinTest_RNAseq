#################################################
##                                             ##
## 82-Cattle-16 Project RNA-seq DGE analysis   ##
## PPD skin testing                            ##
## Unstimulated whole blood direct collections ##
##                                             ##
#################################################

#############################
# Set up required packages: #
#############################

### Install required packages:
#source("https://bioconductor.org/biocLite.R") # will install most updated compatible versions
#biocLite("Rsubread")
#biocLite("edgeR")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("edgeR")

#install.packages("dplyr")
#install.packages("tidyverse")
#install.packages("VennDiagram")
#install.packages("gplots")
#install.packages("pca3d")
#install.packages("reshape2")

### Load required packages: 
library(edgeR) # using version 1.30.0
library(Rsubread) # using version 3.22.0
library(dplyr) # using version 0.8.0.1
library(tidyverse) # using version 1.2.1
library(VennDiagram) # using version 1.6.20
library(gplots) # using version 3.0.1.1
library(pca3d) # using version 0.10
library(reshape2) # using version 1.4.3

########################################
# Load and organize count data into R: #
########################################

### Upload count data and target file:

# Load counts data back into R:
counts <- read.delim("/Users/jayne.wiarda/Desktop/TB/Thesis_Project/82_Cattle_RNAseq/Skin_test_RNAseq/82Cattle_gene_counts_skin_test.txt")
colnames(counts) 

########################
# EdgeR DGE Analysis:: #
########################

### Samples are from pre-skin testing (d0), 1 day post-skin testing (d1), and 3 days post-skin testing (d3)


# Define treatment factor:
group <- factor(rep(c("d0", "d1", "d3"), times = 5))
group

# Define the experimental design and generalized linear model fit:
design <- model.matrix(~ 0 + group) # defining using group-means parameterization
design

# Create the DGE object:
dim(counts)
DGE <- DGEList(counts = counts[,2:16], genes = counts[,1], group = group)
DGEold <- DGE # use this later

### Calculate the normalization factors (TMM method):

DGE <- calcNormFactors(DGE)
DGE$samples # shows adjusted normalization factors

### Obtain raw gene count densities for quality check:
rawDGE <- log10(DGE$counts[,1:ncol(DGE$counts)] + 1) # log transform the raw data

### Filter out lowly expressed genes:

# Filter out genes with less than 0.5 counts per million (cpm) in at least 5 samples:
keep <- rowSums(cpm(DGE) > 1) >= 5

table(keep) # shows how many genes are kept vs. filtered out

DGE <- DGE[keep, , keep.lib.sizes = FALSE] # removes lowly expressed genes from DGE object

write.table(DGE$genes, "/Users/jayne.wiarda/Desktop/TB/Thesis_Project/82_Cattle_RNAseq/Skin_test_RNAseq/skin_test_filtered_genes.txt", sep = "\t")

head(DGE$counts, n = 5) # check that lowly expressed genes are filtered out

DGE$samples # shows new library sizes

# Plot the old vs new library sizes after filtering out low genes:
mean(DGEold$samples$lib.size)
mean(DGE$samples$lib.size) # average library size remains similar

### Plot raw and post-filtering gene count densities for quality check:

filteredDGE <- log10(DGE$counts[,1:ncol(DGE$counts)] + 1) # log transform the filtered data

plot(density(rawDGE[,1]),
     main = "Gene Count Densities Unfiltered vs Filtered Reads",
     lty = 1,
     lwd = 2,
     xlab = "log10 gene counts",
     ylab = "Density",
     col = "blue4",
     ylim = c(0.0,1.0))
lines(density(filteredDGE[,1]), 
      col = "lightskyblue",
      lty = 1,
      lwd = 2)
legend("topright", 
       legend = c("Unfiltered Reads", "Filtered Reads"), 
       col = c("blue4", "lightskyblue"), 
       lty = 1, 
       lwd =2)

### Create MDS plot of all samples:
col <- rep(c("red", "blue", "green"), times = 5)
pch <- c(21,21,21,22,22,22,23,23,23,24,24,24,25,25,25)

plotMDS(DGE, 
        top = 500,
        col = col [as.factor(group)],
        pch = 16,
        cex = 2)
legend("bottomright", 
       legend=c("d0", "d1", "d3"), 
       col = col, 
       pch = c(16,16,16),
       cex = 2,
       pt.cex = 2)

plotMDS(DGE, #plot samples from each animal as a different symbol
        top = 500,
        col = col [as.factor(group)],
        pch = pch,
        cex = 2)

# Estimate common dispersion (required to estimate tagwise dispersion):
DGE <- estimateGLMCommonDisp(DGE, design, verbose = TRUE) # output will provide common dispersion & biological coefficient of variation (BCV)

# Estimate trended dispersion:
DGE <- estimateGLMTrendedDisp(DGE, design)

# Estimate tagwise dispersion (will be used for DE analysis):
DGE <- estimateGLMTagwiseDisp(DGE, design)

### Plot the estimated dispersions:
plotBCV(DGE)

### Test for DGE using GLM:

# Define GLM:
fit <- glmFit(DGE, design)

### Define comparison contrasts & obtain DEGs:
colnames(design) # locate which columns the groups of interest for comparison are in
# for reference:
# genes with positive logFC will be upregulated in groups with contrast = 1 compared to group with -1 contrast
# genes with negative logFC will be downregulated in groups with contrast = 1 compared to group with -1 contrast

### d0 vs d1:
# Perform likelihood ratio test & identify DEGs:
d0v1lrt <- glmLRT(fit, contrast = c(-1, 1, 0)) # perform likelihood ratio test 
summary(d0v1DGE <- decideTestsDGE(d0v1lrt)) # reports upregulated, downregulated, and not DE genes
d0v1topTags <- topTags(d0v1lrt, n = Inf) # table of genes sorted by FDR value (FDR < .05 means DE)
d0v1DE <- subset(d0v1topTags$table, FDR < .05) # create a new table of only DE (FDR < .05) genes
d0v1Up <- subset(d0v1DE, logFC > 0)
d0v1Down <- subset(d0v1DE, logFC < 0)
head(d0v1DE)
# Create plot smear:
d0v1DEnames <- rownames(d0v1DE)
plotSmear(d0v1lrt, de.tags = d0v1DEnames, cex = .75, ylim = c(-6,6))
abline(h = c(-1, 1), col = "blue") 
write.table(d0v1DE, "/Users/jayne.wiarda/Desktop/TB/Thesis_Project/82_Cattle_RNAseq/Skin_test_RNAseq/d0v1_allDE.txt", sep = "\t")
write.table(d0v1Up, "/Users/jayne.wiarda/Desktop/TB/Thesis_Project/82_Cattle_RNAseq/Skin_test_RNAseq/d0v1_Up.txt", sep = "\t")
write.table(d0v1Down, "/Users/jayne.wiarda/Desktop/TB/Thesis_Project/82_Cattle_RNAseq/Skin_test_RNAseq/d0v1_Down.txt", sep = "\t")

### d0 vs d3:
# Perform likelihood ratio test & identify DEGs:
d0v3lrt <- glmLRT(fit, contrast = c(-1, 0, 1)) # perform likelihood ratio test 
summary(d0v3DGE <- decideTestsDGE(d0v3lrt)) # reports upregulated, downregulated, and not DE genes
d0v3topTags <- topTags(d0v3lrt, n = Inf) # table of genes sorted by FDR value (FDR < .05 means DE)
d0v3DE <- subset(d0v3topTags$table, FDR < .05) # create a new table of only DE (FDR < .05) genes
d0v3Up <- subset(d0v3DE, logFC > 0)
d0v3Down <- subset(d0v3DE, logFC < 0)
head(d0v3DE)
# Create plot smear:
d0v3DEnames <- rownames(d0v3DE)
plotSmear(d0v3lrt, de.tags = d0v3DEnames, cex = .75, ylim = c(-6,6))
abline(h = c(-1, 1), col = "blue") 
write.table(d0v3DE, "/Users/jayne.wiarda/Desktop/TB/Thesis_Project/82_Cattle_RNAseq/Skin_test_RNAseq/d0v3_allDE.txt", sep = "\t")
write.table(d0v3Up, "/Users/jayne.wiarda/Desktop/TB/Thesis_Project/82_Cattle_RNAseq/Skin_test_RNAseq/d0v3_Up.txt", sep = "\t")
write.table(d0v3Down, "/Users/jayne.wiarda/Desktop/TB/Thesis_Project/82_Cattle_RNAseq/Skin_test_RNAseq/d0v3_Down.txt", sep = "\t")

### d1 vs d3:
# Perform likelihood ratio test & identify DEGs:
d1v3lrt <- glmLRT(fit, contrast = c(0, -1, 1)) # perform likelihood ratio test 
summary(d1v3DGE <- decideTestsDGE(d1v3lrt)) # reports upregulated, downregulated, and not DE genes
d1v3topTags <- topTags(d1v3lrt, n = Inf) # table of genes sorted by FDR value (FDR < .05 means DE)
d1v3DE <- subset(d1v3topTags$table, FDR < .05) # create a new table of only DE (FDR < .05) genes
d1v3Up <- subset(d1v3DE, logFC > 0)
d1v3Down <- subset(d1v3DE, logFC < 0)
head(d1v3DE)
# Create plot smear:
d1v3DEnames <- rownames(d1v3DE)
plotSmear(d1v3lrt, de.tags = d1v3DEnames, cex = .75, ylim = c(-6,6))
abline(h = c(-1, 1), col = "blue") 
write.table(d1v3DE, "/Users/jayne.wiarda/Desktop/TB/Thesis_Project/82_Cattle_RNAseq/Skin_test_RNAseq/d1v3_allDE.txt", sep = "\t")
write.table(d1v3Up, "/Users/jayne.wiarda/Desktop/TB/Thesis_Project/82_Cattle_RNAseq/Skin_test_RNAseq/d1v3_Up.txt", sep = "\t")
write.table(d1v3Down, "/Users/jayne.wiarda/Desktop/TB/Thesis_Project/82_Cattle_RNAseq/Skin_test_RNAseq/d1v3_Down.txt", sep = "\t")

### Venn Diagram comparison of DE genes:
vd <- venn.diagram(x = list("d0v1" = d0v1Up$genes, "d0v3" = d0v3Up$genes), filename = NULL, fill = c("blue", "green"))
grid.draw(vd, recording = TRUE)
vd <- venn.diagram(x = list("d0v1" = d0v1Down$genes, "d0v3" = d0v3Down$genes), filename = NULL, fill = c("blue", "green"))
grid.draw(vd, recording = TRUE)

######################################################################################################################

### Now run same analysis but omit animal #846:

### Upload count data and target file:

# Load counts data back into R:
#counts <- read.delim("/Users/jayne.wiarda/Desktop/TB/Thesis_Project/82_Cattle_RNAseq/Skin_test_RNAseq/82Cattle_gene_counts_skin_test.txt")
colnames(counts) 
counts <- counts %>% select(-(2:4))
colnames(counts) 

########################
# EdgeR DGE Analysis:: #
########################

### Samples are from pre-skin testing (d0), 1 day post-skin testing (d1), and 3 days post-skin testing (d3)


# Define treatment factor:
group <- factor(rep(c("d0", "d1", "d3"), times = 4))
group

# Define the experimental design and generalized linear model fit:
design <- model.matrix(~ 0 + group) # defining using group-means parameterization
design

# Create the DGE object:
dim(counts)
DGE <- DGEList(counts = counts[,2:13], genes = counts[,1], group = group)
DGEold <- DGE # use this later

### Calculate the normalization factors (TMM method):

DGE <- calcNormFactors(DGE)
DGE$samples # shows adjusted normalization factors

### Obtain raw gene count densities for quality check:
rawDGE <- log10(DGE$counts[,1:ncol(DGE$counts)] + 1) # log transform the raw data

### Filter out lowly expressed genes:

# Filter out genes with less than 0.5 counts per million (cpm) in at least 5 samples:
keep <- rowSums(cpm(DGE) > 1) >= 4

table(keep) # shows how many genes are kept vs. filtered out

DGE <- DGE[keep, , keep.lib.sizes = FALSE] # removes lowly expressed genes from DGE object

write.table(DGE$genes, "/Users/jayne.wiarda/Desktop/TB/Thesis_Project/82_Cattle_RNAseq/Skin_test_RNAseq/no846_skin_test_filtered_genes.txt", sep = "\t")

head(DGE$counts, n = 5) # check that lowly expressed genes are filtered out

DGE$samples # shows new library sizes

# Plot the old vs new library sizes after filtering out low genes:
mean(DGEold$samples$lib.size)
mean(DGE$samples$lib.size) # average library size remains similar

### Plot raw and post-filtering gene count densities for quality check:

filteredDGE <- log10(DGE$counts[,1:ncol(DGE$counts)] + 1) # log transform the filtered data

plot(density(rawDGE[,1]),
     main = "Gene Count Densities Unfiltered vs Filtered Reads",
     lty = 1,
     lwd = 2,
     xlab = "log10 gene counts",
     ylab = "Density",
     col = "blue4",
     ylim = c(0.0,1.0))
lines(density(filteredDGE[,1]), 
      col = "lightskyblue",
      lty = 1,
      lwd = 2)
legend("topright", 
       legend = c("Unfiltered Reads", "Filtered Reads"), 
       col = c("blue4", "lightskyblue"), 
       lty = 1, 
       lwd =2)

### Create MDS plot of all samples:
col <- rep(c("red", "blue", "green"), times = 4)
pch <- c(21,21,21,22,22,22,23,23,23,24,24,24)

plotMDS(DGE, 
        top = 500,
        col = col [as.factor(group)],
        pch = 16,
        cex = 2)
legend("bottomright", 
       legend=c("d0", "d1", "d3"), 
       col = col, 
       pch = c(16,16,16),
       cex = 2,
       pt.cex = 2)

plotMDS(DGE, #plot samples from each animal as a different symbol
        top = 500,
        col = col [as.factor(group)],
        pch = pch,
        cex = 2)

# Estimate common dispersion (required to estimate tagwise dispersion):
DGE <- estimateGLMCommonDisp(DGE, design, verbose = TRUE) # output will provide common dispersion & biological coefficient of variation (BCV)

# Estimate trended dispersion:
DGE <- estimateGLMTrendedDisp(DGE, design)

# Estimate tagwise dispersion (will be used for DE analysis):
DGE <- estimateGLMTagwiseDisp(DGE, design)

### Plot the estimated dispersions:
plotBCV(DGE)

### Test for DGE using GLM:

# Define GLM:
fit <- glmFit(DGE, design)

### Define comparison contrasts & obtain DEGs:
colnames(design) # locate which columns the groups of interest for comparison are in
# for reference:
# genes with positive logFC will be upregulated in groups with contrast = 1 compared to group with -1 contrast
# genes with negative logFC will be downregulated in groups with contrast = 1 compared to group with -1 contrast

### d0 vs d1:
# Perform likelihood ratio test & identify DEGs:
d0v1lrt <- glmLRT(fit, contrast = c(-1, 1, 0)) # perform likelihood ratio test 
summary(d0v1DGE <- decideTestsDGE(d0v1lrt)) # reports upregulated, downregulated, and not DE genes
d0v1topTags <- topTags(d0v1lrt, n = Inf) # table of genes sorted by FDR value (FDR < .05 means DE)
d0v1DE <- subset(d0v1topTags$table, FDR < .05) # create a new table of only DE (FDR < .05) genes
d0v1Up <- subset(d0v1DE, logFC > 0)
d0v1Down <- subset(d0v1DE, logFC < 0)
head(d0v1DE)
# Create plot smear:
d0v1DEnames <- rownames(d0v1DE)
plotSmear(d0v1lrt, de.tags = d0v1DEnames, cex = .75, ylim = c(-6.5,6.5))
abline(h = c(-1, 1), col = "blue") 
write.table(d0v1DE, "/Users/jayne.wiarda/Desktop/TB/Thesis_Project/82_Cattle_RNAseq/Skin_test_RNAseq/no846_d0v1_allDE.txt", sep = "\t")
write.table(d0v1Up, "/Users/jayne.wiarda/Desktop/TB/Thesis_Project/82_Cattle_RNAseq/Skin_test_RNAseq/no846_d0v1_Up.txt", sep = "\t")
write.table(d0v1Down, "/Users/jayne.wiarda/Desktop/TB/Thesis_Project/82_Cattle_RNAseq/Skin_test_RNAseq/no846_d0v1_Down.txt", sep = "\t")

### d0 vs d3:
# Perform likelihood ratio test & identify DEGs:
d0v3lrt <- glmLRT(fit, contrast = c(-1, 0, 1)) # perform likelihood ratio test 
summary(d0v3DGE <- decideTestsDGE(d0v3lrt)) # reports upregulated, downregulated, and not DE genes
d0v3topTags <- topTags(d0v3lrt, n = Inf) # table of genes sorted by FDR value (FDR < .05 means DE)
d0v3DE <- subset(d0v3topTags$table, FDR < .05) # create a new table of only DE (FDR < .05) genes
d0v3Up <- subset(d0v3DE, logFC > 0)
d0v3Down <- subset(d0v3DE, logFC < 0)
head(d0v3DE)
# Create plot smear:
d0v3DEnames <- rownames(d0v3DE)
plotSmear(d0v3lrt, de.tags = d0v3DEnames, cex = .75, ylim = c(-6.5,6.5))
abline(h = c(-1, 1), col = "blue") 
write.table(d0v3DE, "/Users/jayne.wiarda/Desktop/TB/Thesis_Project/82_Cattle_RNAseq/Skin_test_RNAseq/no846_d0v3_allDE.txt", sep = "\t")
write.table(d0v3Up, "/Users/jayne.wiarda/Desktop/TB/Thesis_Project/82_Cattle_RNAseq/Skin_test_RNAseq/no846_d0v3_Up.txt", sep = "\t")
write.table(d0v3Down, "/Users/jayne.wiarda/Desktop/TB/Thesis_Project/82_Cattle_RNAseq/Skin_test_RNAseq/no846_d0v3_Down.txt", sep = "\t")

### d1 vs d3:
# Perform likelihood ratio test & identify DEGs:
d1v3lrt <- glmLRT(fit, contrast = c(0, -1, 1)) # perform likelihood ratio test 
summary(d1v3DGE <- decideTestsDGE(d1v3lrt)) # reports upregulated, downregulated, and not DE genes
d1v3topTags <- topTags(d1v3lrt, n = Inf) # table of genes sorted by FDR value (FDR < .05 means DE)
d1v3DE <- subset(d1v3topTags$table, FDR < .05) # create a new table of only DE (FDR < .05) genes
d1v3Up <- subset(d1v3DE, logFC > 0)
d1v3Down <- subset(d1v3DE, logFC < 0)
head(d1v3DE)
# Create plot smear:
d1v3DEnames <- rownames(d1v3DE)
plotSmear(d1v3lrt, de.tags = d1v3DEnames, cex = .75, ylim = c(-6.5,6.5))
abline(h = c(-1, 1), col = "blue") 
write.table(d1v3DE, "/Users/jayne.wiarda/Desktop/TB/Thesis_Project/82_Cattle_RNAseq/Skin_test_RNAseq/no846_d1v3_allDE.txt", sep = "\t")
write.table(d1v3Up, "/Users/jayne.wiarda/Desktop/TB/Thesis_Project/82_Cattle_RNAseq/Skin_test_RNAseq/no846_d1v3_Up.txt", sep = "\t")
write.table(d1v3Down, "/Users/jayne.wiarda/Desktop/TB/Thesis_Project/82_Cattle_RNAseq/Skin_test_RNAseq/no846_d1v3_Down.txt", sep = "\t")

### Venn Diagram comparison of DE genes:
vd <- venn.diagram(x = list("d0v1" = d0v1Up$genes, "d0v3" = d0v3Up$genes), filename = NULL, fill = c("blue", "green"))
grid.draw(vd, recording = TRUE)
vd <- venn.diagram(x = list("d0v1" = d0v1Down$genes, "d0v3" = d0v3Down$genes), filename = NULL, fill = c("blue", "green"))
grid.draw(vd, recording = TRUE)

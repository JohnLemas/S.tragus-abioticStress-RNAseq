#######################################

#__Date:__ January 21, 2024
#__Author:__ Erin Osborne Nishimura, John M. Lemas
#__Script:__ SalTr_DESeq2.R
#__Project:__ To assess expression differences between salt and cold treated Salsola tragus plants. 
#__Requires:__ 
# 
# + R (4.2.2)
# + DESeq2 (1.38.1)   
# + corrplot (0.92)
# + RColorBrewer (1.1-3)
# + pheatmap (1.0.12)
# + apeglm (1.20.0)

######################################

###########  LOAD PACKAGE  #####################

# Do this every time.
# Load packages:
library(DESeq2)
library(corrplot)
library(RColorBrewer)
library(pheatmap)
library(apeglm)
library(RColorBrewer)
library(ggplot2)
library(kableExtra)
library(emmeans)
#Newbies from the vignette:
library(pcaExplorer)
library(vsn)
library(RNAseqQC)
library(ensembldb)
library(dplyr)
library(purrr)
library(tidyr)
library(tibble)
library(magrittr)
library(EnhancedVolcano)
library(extrafont)

font_import()
y
loadfonts()
#################################################

setwd("~/SalTrGenome/RNA_Seq/RT_StressTrials/RNA-seq_data/")

# Read in the counts data:
countsData <- read.table(file = "03_output/03_feature/STAR_counts_v6.txt", header = FALSE, row.names = 1, skip = 2) # 

# Read in the metadata:
metadata_full <- read.table(file = "01_input/S.tragus_metadata_v2.txt", header = FALSE)

# Insert header for metadata:
colnames(metadata_full) <- c("fasta1", "fasta2", "names1", "names2", "type", "trt", "time", "rep")

# insert header for counts data:
colnames(countsData) <- c("chr", "start", "stop", "strand", "length", as.vector(metadata_full$names1))

# Subset counts data for analysis:

countColNum <- dim(countsData)[2]

head(countsData[,6:countColNum])
dim(countsData[,6:countColNum])

# Ensure that this will capture all of the data!
# We only want the sample names and the counts for each gene. Save just the subset as an object called cts:
cts <- as.matrix(countsData[,6:countColNum])

# Below are just the shoot samples with only one set of shoot control samples at 24 HAT

cts <- cts[,c(4:6,28:30,34:36)]

# Subset the metadata for the analysis:
metadata1 <- metadata_full[c(4:6,28:30,34:36),]

# Create row headers with selected labels:
rownames(metadata1)<- metadata1$names1

# Subset selected factors into a new object:
coldata <- metadata1[,c("trt", "time", "rep")]
coldata$trt <- as.factor(coldata$trt)
coldata$time <- as.factor(coldata$time)
coldata$rep <- as.factor(coldata$rep)

#Check that the names match  --> Should be TRUE
all(rownames(coldata) == colnames(cts))

# Create dds object from matrix:
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ time)
# Filter for low counts:
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# Designate factor levels:

dds$time <- factor(dds$time, levels = c("0", "24"))
dds$trt <- factor(dds$trt, levels = c("Contl", "NaCl", "Cold"))
dds$rep <- factor(dds$rep, levels = c("1", "2", "3"))

# PERFORM DESEQ2 analysis:

dds <- DESeq(dds)

plotDispEsts(dds)

# Save for reference:
pdf("03_output/05_DESeq2/Figure_files/Root/dispEst.pdf", height = 6, width = 8)
par(mfrow=c(1,1))

plotDispEsts(dds)

dev.off()
dev.off()

# Both of these datasets are good supplemental tables for papers

head(counts(dds, normalized = TRUE))
head(counts(dds, normalized = FALSE))

counts_norm <- counts(dds, normalized = TRUE)
counts_nonNorm <- counts(dds, normalized = FALSE)

write(counts_norm, file = "03_output/05_DESeq2/03_cts_files/Root_24_cts_norm.txt", sep = "\t")
write(counts_nonNorm, file = "03_output/05_DESeq2/03_cts_files/Root_24_cts_nonNorm.txt", sep = "\t")

#' upload both of these tables for supplimental information in the manuscript


# Plot some QC figures:
plot_total_counts(dds)
plot_library_complexity(dds)

pdf("03_output/05_DESeq2/Figure_files/Root/total_counts.pdf", height = 6, width = 8)
par(mfrow=c(1,1))

plot_total_counts(dds)

dev.off()
dev.off()

pdf("03_output/05_DESeq2/Figure_files/Root/library_complx.pdf", height = 6, width = 8)
par(mfrow=c(1,1))

plot_library_complexity(dds)

dev.off()
dev.off()

# Need to calibrate for log fold shrinkage:

# I followed this work through: https://www.biostars.org/p/448959/#484944
# originally I did not have output from resultNames(res) so I used this work around
# here is the edited scripting to get lfcShrink to work:

dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)
dds <- nbinomWaldTest(dds)

relevel(dds$time, ref = "0")
relevel(dds$trt, ref = "Contl")
# Now when I run resultNames(dds) the reference level is defined as above

resultsNames(dds)

# Now we need to re-reun the Wald test:

dds <- nbinomWaldTest(dds)
resultsNames(dds)

resLFC_0_24 <- lfcShrink(dds, coef = 2, type = "apeglm")

#vignette output from variance stabalization:

ntd <- normTransform(dds)

meanSdPlot(assay(ntd), ylab = "Standard Deviation", xlab = "Rank (mean)")

# Save for reference:
pdf("03_output/05_DESeq2/Figure_files/Root/sdVmean_normTrans.pdf", height = 6, width = 8)
par(mfrow=c(1,1))

meanSdPlot(assay(ntd), ylab = "Standard Deviation", xlab = "Rank (mean)")

dev.off()
dev.off()

vsd <- vst(dds)

meanSdPlot(assay(vsd))

# Save for reference:
pdf("03_output/05_DESeq2/Figure_files/Root/sdVmean.pdf", height = 6, width = 8)
par(mfrow=c(1,1))

meanSdPlot(assay(vsd))

dev.off()
dev.off()


#' variance stabalization transformation (vst; as above with meanSDPlot).
#' the assumption is that greater means equal greater variance, but for 
#' count data we want to had similar variance for all counts regardless of mean
#' so we tranform it so that increasing mean rank doesn't increase the variance.
#' So why the two options? The two previous plots look quite different. The
#' first normalizes before transforming? Why does it have a bump at the start?

# The vignette also has some pca stuff in it! Just use the vsd object:

plot_pca(vsd, PC_x = 1, PC_y = 2, color_by = "rep", shape_by = "trt")

pca <- plot_pca(vsd, PC_x = 1, PC_y = 2, color_by = "rep", shape_by = "trt")

pdf("03_output/05_DESeq2/Figure_files/Root/PCA.pdf", height = 6, width = 8)
par(mfrow=c(1,1))

plot_pca(vsd, PC_x = 1, PC_y = 2, color_by = "rep", shape_by = "trt")

dev.off()
dev.off()

# Notes from Jake:
# create a table with the loadings for each PC
loadings <- as.data.frame(pca$loadings)

#export loadings into a file
write.csv(x = loadings, file = "03_output/05_DESeq2/Figure_files/Root/PCA_loadings.csv")

# I wonder if this would look better with all treatments and timepoints.
# Maybe I'll generate some QC plots for the whole exp 

# I need to split this up into shoot and root samples because the expression
# between tissue types is throwing everything off. 

# This looks much better after splitting shoot and root samples. 

##### CORRELATION MATRICES:

#Take r-stabilized log transformations of all the normalized count data.
rld <- rlog(dds, blind=FALSE)

# Calculate the distances between each sample
sampleDists <- dist(t(assay(rld)))


sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- colnames(rld) # Add some labels
#colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255) # Go back and alter colors later

# Draw the heatmap
par(mfrow=c(1,1))
p <- pheatmap(sampleDistMatrix,
              clustering_distance_rows=sampleDists,
              clustering_distance_cols=sampleDists,
              col=colors) # Plot the heatmap

# Save the CORRELATION MATRIX as a .pdf file


pdf("03_output/05_DESeq2/Figure_files/Root/heatmap.pdf", height = 6, width = 7)
par(mfrow=c(1,1))
p

dev.off()
dev.off() 

############## DIFFERENTIAL EXPRESSION ANALYSIS #####################

# calculate the statistically significant differences
# Set your contrast expression: Contrast - a character vector with exactly three elements: the name of a factor in the design formula, the name of the numerator level for the fold change, and the name of the denominator level for the fold change (simplest case)
# Set your lfcThreshold to a value of your choosing. This will set whether you want to exclude any genes as significant because their fold change is too low. Default is 0. That means, take all statistically significant genes. In this dataset, I selected 0.5 because there were A LOT of significant genes and I'd like to reduce the number.
# Set your alpha value. This is your p-value cut-off for significance. The default is 0.1. 

res_0_24 <- results(dds,
                    lfcThreshold = 0.5,
                    alpha = 0.05,
                    contrast=c("time", "24", "0"))

summary(res_0_24)


# Lets check out the anitporter that I was looking at before:

plotCounts(dds, gene=which(rownames(res_0_24) == "SalTrChr04Ag088810"), intgroup="time")

plot_gene("SalTrChr04Ag088810", dds, x_var = "time", color_by = "trt")

pdf("03_output/05_DESeq2/Figure_files/Root/Antiporter.pdf", height = 6, width = 8)
par(mfrow=c(1,1))

plot_gene("SalTrChr04Ag088810", dds, x_var = "trt")

dev.off()
dev.off()

########################
# NHX1 vacuolar NA+/H+ antiporter:

plotCounts(dds, gene=which(rownames(res_0_24) == "SalTrChr04Ag088810"), intgroup=c("time", "trt"))

antiporter_countsplot <- plotCounts(dds, gene=which(rownames(res_0_24) == "SalTrChr04Ag088810"), intgroup=c("time", "trt"))

antiporter_counts_table <- plotCounts(dds, gene=which(rownames(res_0_24) == "SalTrChr04Ag088810"), intgroup=c("time", "trt"), returnData = TRUE)

antiporter_table <- data.frame(sample = metadata1$names1,
                               treatment = c("0 HAT", "0 HAT", "0 HAT",
                                             "24 HAT", "24 HAT", "24 HAT",
                                             "24 HAT", "24 HAT", "24 HAT"),
                               count = antiporter_counts_table$count,
                               time = antiporter_counts_table$time,
                               trt = antiporter_counts_table$trt)

ggplot(data = antiporter_table, mapping = aes(
  x = trt, y = count)) + 
  geom_point() + theme_light() + 
  xlab("Hours after treatment") + ylab("Normalized gene count") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        plot.caption = element_text(hjust = 0.01),
        text = element_text(family = "Times New Roman", size = 12)) +
  labs(title = "Salt treated root tissue", subtitle = "SalTrChr04Ag088810\nNHX1: Vacuolar Na+/H+ Exchanger",
       caption = "Tukey adjusted pairwise comparison: Contl - NaCl: p = 0.1; Contl - Cold: p = 0.001;\nNaCl - Cold: p = 0.012")

# save this gene count plot

antiporter_boxplot <- ggplot(data = antiporter_table, mapping = aes(
  x = trt, y = count)) + 
  geom_point() + theme_light() + 
  xlab("Hours after treatment") + ylab("Normalized gene count") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        plot.caption = element_text(hjust = 0.01),
        text = element_text(family = "Times New Roman", size = 12)) +
  labs(title = "Salt treated root tissue", subtitle = "SalTrChr04Ag088810\nNHX1: Vacuolar Na+/H+ Exchanger",
       caption = "Tukey adjusted pairwise comparison: Contl - NaCl: p = 0.1; Contl - Cold: p = 0.001;\nNaCl - Cold: p = 0.012")


pdf(file = "03_output/05_DESeq2/Figure_files/Root/NHX1_boxplot.pdf", width = 6, height = 6)

par(mfrow=c(1,1))

antiporter_boxplot

dev.off()
dev.off()

# Pairwise comparison for each timepoint

antiporter_anova <- anova(lm(count ~ trt, data = antiporter_table))

antiporter_emmeans <- emmeans(lm(count ~ trt, data = antiporter_table),
                              ~ trt)

antiporter_pairs <- pairs(
  antiporter_emmeans, adjust = "Tukey"
)

kbl(antiporter_pairs, 
    caption = "NA+/H+ Echanger ANOVA", 
    booktabs = T,
    digits = 3,
) %>%
  kable_styling(latex_options = c("striped", "hold_position"))


################ Gene of Interest ID from lowest adj. p-value:

# Lets look for the gene with the lowest p value and do some research:

plotCounts(dds, gene=which.min(res_0_24$padj), intgroup="time")
# Gene: SalTrChr03Bg286720

pdf("03_output/05_DESeq2/Figure_files/Root/padj_min_dotplot.pdf", height = 6, width = 8)
par(mfrow=c(1,1))

plot_gene("SalTrChr03Bg286720", dds, x_var = "time", color_by = "trt")

dev.off()
dev.off()

##########################
# This one shows upregulation for both treatments!

plot_gene("SalTrChr03Ag060690", dds, x_var = "time", color_by = "trt")

pdf("03_output/05_DESeq2/Figure_files/Root/SalTrChr03Ag060690_dotplot.pdf", height = 6, width = 8)
par(mfrow=c(1,1))

plot_gene("SalTrChr03Ag060690", dds, x_var = "time", color_by = "trt")

dev.off()
dev.off()

# Looks very interesting! I'm blasting the protein sequence of this gene into 
# Interpro to see what I can dig up.

##########################
# Gene AmCMO from halophyte paper

plot_gene("SalTrChr02Ag029170", dds, x_var = "time", color_by = "trt")

# Whoa cool. Looks more important for cold stress response

plot_gene("SalTrChr03Ag053910", dds, x_var = "time", color_by = "trt")

# That one also appears more important for cold treated tissue

plot_gene("SalTrChr03Ag065710", dds, x_var = "time", color_by = "trt")

# Thats a good one too. Like the second one but with a larger y axis

plot_gene("SalTrChr03Ag065720", dds, x_var = "time", color_by = "trt")

# Even higher on that one. Nice

plot_gene("SalTrChr05Ag102850", dds, x_var = "time", color_by = "trt")

# Only salt again

plot_gene("SalTrChr08Ag192970", dds, x_var = "time", color_by = "trt")

# Super tight super highly expressed and only cold. Cool.

################################################

# Select the significant subset of genes that are up-regulated
Up_0_24 <- subset(res_0_24, padj < 0.1 & log2FoldChange > 0.5)
Up_0_24 <- Up_0_24[order(Up_0_24$padj),] #order them

head(Up_0_24) # Check them
dim(Up_0_24)

# Select the significant subset of genes that are down-regulated
Down_0_24 <- subset(res_0_24, padj < 0.1 & log2FoldChange < -0.5)
Down_0_24 <- Down_0_24[order(Down_0_24$padj),]

head(Down_0_24)
dim(Down_0_24)

# Save these lists to output files:
write(rownames(Up_0_24), file = "03_output/05_DESeq2/Up_root_0_24.txt", sep = "\n")
write(rownames(Down_0_24), file = "03_output/05_DESeq2/Down_root_0_24.txt", sep = "\n")
write(rownames(res_0_24), file = "03_output/05_DESeq2/All_root_0_24.txt", sep = "\n")

###########################################

pdf("03_output/05_DESeq2/Figure_files/Root/Evolcano_root_0_24.pdf",
    width = 8, height = 8)

par(mfrow = c(1,1))

EnhancedVolcano(toptable = res_0_24, lab = rownames(res_0_24), 
                x = 'log2FoldChange',
                y = "pvalue",
                title = 'Differentially Expressed Genes 24 HAT',
                subtitle = 'C: Cold and salt treated root tissue samples',
                caption = bquote(~Log[2]~ "fold change cutoff = ±2; p-value cutoff: p < 0.05"),
                xlab = bquote(~Log[2]~ 'fold change'),
                selectLab = c(
                  "SalTrChr03Bg286720",
                  "SalTrChr03Bg288430",
                  "SalTrChr02Ag043960",
                  "SalTrChr06Bg341440",
                  "SalTrChr02Ag034730",
                  "SalTrChr07Bg361280"
                ),
                pCutoff = 5e-2,
                FCcutoff = 2,
                cutoffLineType = 'twodash',
                cutoffLineWidth = 0.8,
                cutoffLineCol = 'grey',
                hline = 10e-10,
                hlineType = 'twodash',
                hlineCol = 'black',
                hlineWidth = 0.8,
                pointSize = 4.0,
                labSize = 4.5,
                labFace = 'bold',
                boxedLabels = TRUE,
                colAlpha = 0.5,
                legendPosition = 'none',
                drawConnectors = TRUE,
                widthConnectors = 0.5,
                colConnectors = 'black')

dev.off()
dev.off()


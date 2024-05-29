#######################################

#__Date:29 January 2024
#__Author:__ John Lemas
#__Script:__ 
#__Project:__
#__Developed with the following versions:__ 
# 
# + R (4.3.1)
# + DESeq2 (1.42.1)   
# + corrplot (0.92)
# + RColorBrewer (1.1-3)
# + pheatmap (1.0.12)
# + apeglm (1.24.0)
# + tidyerse (2.0.0)

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
#library(tidyverse) # causes some minor conflicts, but you can install this later for creating plots

font_import()
y
loadfonts()

#################################################

# Set the working directory:

setwd("~/SalTrGenome/RNA_Seq/RT_StressTrials/RNA-seq_data/")

# Import the counts data file

countsData <- read.table(file = "03_output/03_feature/STAR_counts_v6.txt", header = TRUE, row.names = 1) # 

# Read in the metadata
metadata_full <- read.table(file = "01_input/S.tragus_metadata_v2.txt", header = FALSE)



# Add headers

colnames(metadata_full) <- c("fasta1", "fasta2", "names1", "names2", "type", "trt", "time", "rep")

# Now add headers for the counts file. the final columns should be the samples in order from the metadata

as.vector(metadata_full$names1)

# Name countsData columns headers:

colnames(countsData) <- c("chr", "start", "stop", "strand", "length", as.vector(metadata_full$names1))

################### COUNT MATRIX INPUT ###################

# In this section we will prepare our input data for analysis.
# In the instruction for DESeq2, it states: "We read in a count matrix, which we will name cts, and the sample information table, which we will name coldata."

# Subset the countsData 
# Here we can use dim() to grab the number of columns from the counts file. 
# dim() spits out rows first [1] and columns second [2]

countColNum <- dim(countsData)[2]

head(countsData[,6:countColNum])
dim(countsData[,6:countColNum])

# Ensure that this will capture all of the data!
# We only want the sample names and the counts for each gene. Save just the subset as an object called cts:
cts <- as.matrix(countsData[,6:countColNum])

# I need to filter the cts object for just the contrast that I am running in this script. featureCounts output
# includes all 36 samples, so filtering the metadata files wont cut it. I'll have to do
# it manually in R. Enter a vector of indexes to capture the specific contrast- ie c(1:12, 19:24, 31:36)
# Edit the following based on the contrast for this run:

cts_contrast <- cts[,c(1:6,25:36)]

# Now instead of indexing the metadata accordingly I can import the metadata file I need for this run that was generated in Unix:

metadata_contrast <- metadata_full[c(1:6,25:36),]

# Next we need to make an information called coltable. We can make this out of the metadata table.
# Reorganize the metadata table so the names2 column are now row headers

rownames(metadata_contrast) <- metadata_contrast$names1

# Adjust the following based on the contrasts included if needed:

coldata <- metadata_contrast[,c("type", "trt", "time", "rep")]
coldata$type <- as.factor(coldata$type)
coldata$trt <- as.factor(coldata$trt)
coldata$time <- as.factor(coldata$time)
coldata$rep <- as.factor(coldata$rep)
str(coldata)

# Double check that row names and column names are identical:

rownames(coldata)
colnames(cts)

all(rownames(coldata) == colnames(cts_contrast))


# Next we will create an ddsHTSeq object out of cts and coldata:
# This will set a base design for your experiment. Please insert the design:

dds <- DESeqDataSetFromMatrix(countData = cts_contrast,
                              colData = coldata,
                              design = ~ time + type)


################# PRE-FILTERING -- FILTER FOR PRESENT GENES: ################# 
# Not necessary, but helps keep things fast. 
# Exclude all samples that have less than 10 reads:
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

#' taken from the DESeq2 Vignette. you can adjust this value for your study

# Exercise: How many did we exclude?
dim(dds)

str(dds)

################### NOTE ON FACTOR LEVELS ###################
# Alther the below factor levels based on what the design for this run is:

dds$type <- factor(dds$type, levels = c("Shoot","Root"))
dds$time <- factor(dds$time, levels = c("0", "24"))
dds$trt <- factor(dds$trt, levels = c("Contl", "NaCl", "Cold"))

# PERFORM DESEQ2 analysis:

dds <- DESeq(dds)

plotDispEsts(dds)

pdf("03_output/05_DESeq2/Figure_files/All_24/dispEst.pdf", height = 6, width = 8)
par(mfrow=c(1,1))

plotDispEsts(dds)

dev.off()
dev.off()

dds$sizeFactor

# Library Complexity and counts
plot_total_counts(dds)
plot_library_complexity(dds)

pdf("03_output/05_DESeq2/Figure_files/All_24/total_counts.pdf", height = 7, width = 7)
par(mfrow=c(1,1))

plot_total_counts(dds)

dev.off()
dev.off()

pdf("03_output/05_DESeq2/Figure_files/All_24/library_complx.pdf", height = 7, width = 7)
par(mfrow=c(1,1))

plot_library_complexity(dds)

dev.off()
dev.off()

# Both of these datasets are good supplemental tables for papers

head(counts(dds, normalized = TRUE))
head(counts(dds, normalized = FALSE))

counts_norm <- counts(dds, normalized = TRUE)
counts_nonNorm <- counts(dds, normalized = FALSE)

write.table(counts_norm, file = "03_output/05_DESeq2/03_cts_files/All_24_cts_norm.txt", sep = "\t",
            row.names = TRUE, col.names = TRUE)
write.table(counts_nonNorm, file = "03_output/05_DESeq2/03_cts_files/All_24_cts_nonNorm.txt", sep = "\t",
            row.names = TRUE, col.names = TRUE)

# Lets also include the summary of the normalized counts so that we can see where the expression
# levels are for each sample.

counts_sum <- summary(counts(dds, normalized = TRUE))
write.table(counts_sum, file = "03_output/05_DESeq2/03_cts_files/All_24_cts_norm_summary.txt", sep = "\t",
            row.names = TRUE, col.names = TRUE)

#' upload both of these tables for supplimental information in the manuscript
############## PERFORM LOG FOLD CHANGE SHRINKAGE FOR VISUALIZATION #####################

# An input requirement of the lfcShrink function is a coef term. This is pulled from the resultsNames of dds:
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

# ok now we can run lfcShrink(dds = dds, coef = [2, 3, or 4], type = "ageglm")
# coef = 2 means 0 vs 3; coef = 3 means 0 vs 8; coef = 4 means 0 vs 24

resLFC_0_24 <- lfcShrink(dds, coef = 2, type = "apeglm")

# Set your alpha value. This is your p-value cut-off for significance. The default is 0.1. 

res_0_24 <- results(dds,
                   lfcThreshold = 0.5,
                   alpha = 0.05,
                   contrast=c("time", "24", "0"))

#vignette notes:
# This is pretty and I like it but what does it mean?

ntd <- normTransform(dds)

meanSdPlot(assay(ntd), ylab = "Standard Deviation", xlab = "Rank (mean)")

# Lets save it just for funsies:

sdVmean <- meanSdPlot(assay(ntd), ylab = "Standard Deviation", xlab = "Rank (mean)")

pdf("03_output/05_DESeq2/Figure_files/Cold_shoot/sdVmean_normTrans.pdf", height = 8, width = 7)
par(mfrow=c(1,1))

sdVmean

dev.off()
dev.off() 

# I wish I knew how to change the colors!

vsd <- vst(dds)

meanSdPlot(assay(vsd))

pdf("03_output/05_DESeq2/Figure_files/Cold_shoot/sdVmean.pdf", height = 8, width = 7)
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

plot_pca(vsd, PC_x = 1, PC_y = 2, color_by = "trt", shape_by = "time")

pdf("03_output/05_DESeq2/Figure_files/All_24/PCA.pdf", height = 8, width = 7)
par(mfrow=c(1,1))

plot_pca(vsd, PC_x = 1, PC_y = 2, color_by = "trt", shape_by = "time")

dev.off()
dev.off() 

# create a table with the loadings for each PC

pca <- plot_pca(vsd, PC_x = 1, PC_y = 2, color_by = "trt", shape_by = "time")

loadings <- as.data.frame(pca$loadings)

#export loadings into a file
write.csv(x = loadings, file = "03_output/05_DESeq2/Figure_files/All_24/PCA_loadings.csv")

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


pdf("03_output/05_DESeq2/Figure_files/All_24/corr_matrix_contrast.pdf", height = 6, width = 7)
par(mfrow=c(1,1))
p

dev.off()
dev.off() 

# MA plots

par(mfrow=c(1,2))

plotMA(res_0_24, main="0 versus 24 HAT\nunshrunken", ylim = c(-7,7), 
       ylab = "LFC (24 HAT / 0 HAT)",
       xlab = "means of normalized counts")

plotMA(resLFC_0_24, main="0 versus 24 HAT\nshrunken", ylim = c(-7,7), 
       ylab = "LFC (24 HAT / 0 HAT)",
       xlab = "means of normalized counts")

# Save these plots in one panel:

pdf("03_output/05_DESeq2/Figure_files/All_24/MAplots.pdf", height = 6, width = 8)

par(mfrow=c(1,2))

plotMA(res_0_24, main="0 versus 24 HAT\nunshrunken", ylim = c(-7,7), 
       ylab = "LFC (24 HAT / 0 HAT)",
       xlab = "means of normalized counts")

plotMA(resLFC_0_24, main="0 versus 24 HAT\nshrunken", ylim = c(-7,7), 
       ylab = "LFC (24 HAT / 0 HAT)",
       xlab = "means of normalized counts")

dev.off()
dev.off() 


############## DIFFERENTIAL EXPRESSION ANALYSIS #####################

# calculate the statistically significant differences
# Set your contrast expression: Contrast - a character vector with exactly three elements: the name of a factor in the design formula, the name of the numerator level for the fold change, and the name of the denominator level for the fold change (simplest case)
# Set your lfcThreshold to a value of your choosing. This will set whether you want to exclude any genes as significant because their fold change is too low. Default is 0. That means, take all statistically significant genes. In this dataset, I selected 0.5 because there were A LOT of significant genes and I'd like to reduce the number.
# Set your alpha value. This is your p-value cut-off for significance. The default is 0.1. 


##################  Exploring and exporting results ################## 

# NHX1 Na+/H+ antiporter

plotCounts(dds, gene=which(rownames(res_0_24) == "SalTrChr04Ag088810"), intgroup=c("type", "time", "trt"))

plot_gene("SalTrChr04Ag088810", dds = dds, x_var = "trt", color_by = "type")

antiporter_countsplot <- plotCounts(dds, gene=which(rownames(res_0_24) == "SalTrChr04Ag088810"), intgroup=c("type", "time", "trt"))

antiporter_counts_table <- plotCounts(dds, gene=which(rownames(res_0_24) == "SalTrChr04Ag088810"), intgroup=c("type", "time", "trt"), returnData = TRUE)

antiporter_table <- data.frame(sample = metadata_contrast$names1,
                               treatment = c("0 HAT", "0 HAT", "0 HAT",
                                             "24 HAT", "24 HAT", "24 HAT"),
                               count = antiporter_counts_table$count,
                               time = antiporter_counts_table$time,
                               type = antiporter_counts_table$type,
                               trt = antiporter_counts_table$trt)

ggplot(data = antiporter_table) + 
  geom_point(mapping = aes(
    x = trt, y = count, color = type)) +
  theme_light() + ggtitle("NA+/H+ Echanger") +
  xlab("Hours after treatment") + ylab("Normalized gene count") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        plot.caption = element_text(hjust = 0.01),
        text = element_text(family = "Times New Roman", size = 12)) +
  labs(title = "All tissues 24 HAT", subtitle = "SalTrChr04Ag088810\nNHX1: Vacuolar Na+/H+ Exchanger",
       color = "Tissue type",
       shape = "Treatment") +
  scale_color_brewer(palette = "Dark2") + 
  guides(size = "none", alpha = "none")

# save this gene count plot

antiporter_boxplot <- ggplot(data = antiporter_table) + 
  geom_point(mapping = aes(
    x = trt, y = count, color = type)) +
  theme_light() + ggtitle("NA+/H+ Echanger") +
  xlab("Hours after treatment") + ylab("Normalized gene count") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        plot.caption = element_text(hjust = 0.01),
        text = element_text(family = "Times New Roman", size = 12)) +
  labs(title = "All tissues 24 HAT", subtitle = "SalTrChr04Ag088810\nNHX1: Vacuolar Na+/H+ Exchanger",
       color = "Tissue type",
       shape = "Treatment") +
  scale_color_brewer(palette = "Dark2") + 
  guides(size = "none", alpha = "none")

pdf(file = "03_output/05_DESeq2/Figure_files/All_24/NHX1_boxplot.pdf", width = 6, height = 6)

par(mfrow=c(1,1))

antiporter_boxplot

dev.off()
dev.off()

# Pairwise comparison for each timepoint

antiporter_anova <- anova(lm(count ~ type + trt, data = antiporter_table))

antiporter_emmeans <- emmeans(lm(count ~ type + trt, data = antiporter_table),
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

####### lowest adjusted p-vaue:

plotCounts(dds, gene=which.min(res_0_24$padj), intgroup=c("type", "time", "trt"))
# SalTrChr03Bg288430 heat shock factor protein HSF30-like isoform X1 IPR027725

plotCounts(dds, gene=which(rownames(res_0_24) == "SalTrChr03Bg288430"), intgroup=c("type", "time", "trt"))

HSF_countsplot <- plotCounts(dds, gene=which(rownames(res_0_24) == "SalTrChr03Bg288430"), intgroup=c("type", "time", "trt"))

HSF_counts_table <- plotCounts(dds, gene=which(rownames(res_0_24) == "SalTrChr03Bg288430"), intgroup=c("type", "time", "trt"), returnData = TRUE)

HSF_table <- data.frame(sample = metadata_contrast$names1,
                               treatment = c("0 HAT", "0 HAT", "0 HAT",
                                             "24 HAT", "24 HAT", "24 HAT"),
                               count = HSF_counts_table$count,
                               time = HSF_counts_table$time,
                               type = HSF_counts_table$type,
                               trt = HSF_counts_table$trt)

ggplot(data = HSF_table) + 
  geom_point(mapping = aes(
    x = trt, y = count, color = type)) +
  theme_light() + ggtitle("Heat Shock TF HSF-30 Isoform") +
  xlab("Hours after treatment") + ylab("Normalized gene count") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        plot.caption = element_text(hjust = 0.01),
        text = element_text(family = "Times New Roman", size = 12)) +
  labs(title = "All tissues 24 HAT", subtitle = "SalTrChr03Bg288430\nHeat Shock TF: HSF-30 Isoform",
       color = "Tissue type") +
  scale_color_brewer(palette = "Dark2")

# save this gene count plot

HSF_boxplot <- ggplot(data = HSF_table) + 
  geom_point(mapping = aes(
    x = trt, y = count, color = type)) +
  theme_light() + ggtitle("Heat Shock TF HSF-30 Isoform") +
  xlab("Hours after treatment") + ylab("Normalized gene count") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        plot.caption = element_text(hjust = 0.01),
        text = element_text(family = "Times New Roman", size = 12)) +
  labs(title = "All tissues 24 HAT", subtitle = "SalTrChr03Bg288430\nHeat Shock TF: HSF-30 Isoform",
       color = "Tissue type") +
  scale_color_brewer(palette = "Dark2")

pdf(file = "03_output/05_DESeq2/Figure_files/All_24/HSF_boxplot.pdf", width = 6, height = 6)

par(mfrow=c(1,1))

HSF_boxplot

dev.off()
dev.off()

##### KNOWN GENES:

# BADH no true differential expression. Ex:
plot_gene("SalTrChr09Bg424260", dds, x_var = "time", color_by = "type")

plotCounts(dds, gene=which(rownames(res_0_24) == "SalTrChr09Bg424260"), intgroup=c("type", "time", "trt"))

#Pro1:
plot_gene("SalTrChr03Bg285640", dds, x_var = "trt", color_by = "type")

plotCounts(dds, gene=which(rownames(res_0_24) == "SalTrChr03Bg285640"), intgroup=c("type", "time", "trt"))

# MDHR:
plot_gene("SalTrChr01Ag014560", dds, x_var = "trt", color_by = "type")

plot_gene("SalTrChr01Ag027080", dds, x_var = "trt", color_by = "type")

plotCounts(dds, gene=which(rownames(res_0_24) == "SalTrChr01Ag027080"), intgroup=c("type", "time", "trt"))

####################
plotCounts(dds, gene=which(rownames(res_0_24) == "SalTrChr01Ag027080"), intgroup=c("type", "time", "trt"))

antiporter_countsplot <- plotCounts(dds, gene=which(rownames(res_0_24) == "SalTrChr01Ag027080"), intgroup=c("type", "time", "trt"))

antiporter_counts_table <- plotCounts(dds, gene=which(rownames(res_0_24) == "SalTrChr01Ag027080"), intgroup=c("type", "time", "trt"), returnData = TRUE)

antiporter_table <- data.frame(sample = metadata_contrast$names1,
                               treatment = c("0 HAT", "0 HAT", "0 HAT",
                                             "24 HAT", "24 HAT", "24 HAT"),
                               count = antiporter_counts_table$count,
                               time = antiporter_counts_table$time,
                               type = antiporter_counts_table$type,
                               trt = antiporter_counts_table$trt)

ggplot(data = antiporter_table) + 
  geom_point(mapping = aes(
    x = trt, y = count, color = type)) +
  theme_light() +
  xlab("Hours after treatment") + ylab("Normalized gene count") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        plot.caption = element_text(hjust = 0.01),
        text = element_text(family = "Times New Roman", size = 12)) +
  labs(title = "All tissues 24 HAT", subtitle = "SalTrChr01Ag027080\nAscorbate regeneration and ROS scavaging",
       color = "Tissue type") +
  scale_color_brewer(palette = "Dark2")

# save this gene count plot

antiporter_boxplot <- ggplot(data = antiporter_table) + 
  geom_point(mapping = aes(
    x = trt, y = count, color = type)) +
  theme_light() +
  xlab("Hours after treatment") + ylab("Normalized gene count") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        plot.caption = element_text(hjust = 0.01),
        text = element_text(family = "Times New Roman", size = 12)) +
  labs(title = "All tissues 24 HAT", subtitle = "SalTrChr01Ag027080\nAscorbate regeneration and ROS scavaging",
       color = "Tissue type") +
  scale_color_brewer(palette = "Dark2")

pdf(file = "03_output/05_DESeq2/Figure_files/All_24/MDHR_boxplot.pdf", width = 6, height = 6)

par(mfrow=c(1,1))

antiporter_boxplot

dev.off()
dev.off()

#########################

plot_gene("SalTrChr04Ag086160", dds, x_var = "trt", color_by = "type")

plotCounts(dds, gene=which(rownames(res_0_24) == "SalTrChr04Ag086160"), intgroup=c("type", "time", "trt"))

plot_gene("SalTrChr08Ag189450", dds, x_var = "trt", color_by = "type")

plotCounts(dds, gene=which(rownames(res_0_24) == "SalTrChr08Ag189450"), intgroup=c("type", "time", "trt"))

plot_gene("SalTrChr01Bg243720", dds, x_var = "trt", color_by = "type")

plotCounts(dds, gene=which(rownames(res_0_24) == "SalTrChr01Bg243720"), intgroup=c("type", "time", "trt"))

# V-ATPase:
plot_gene("SalTrChr09Ag208250", dds, x_var = "trt", color_by = "type")

plotCounts(dds, gene=which(rownames(res_0_24) == "SalTrChr09Ag208250"), intgroup=c("type", "time", "trt"))

plot_gene("SalTrChr09Bg423040", dds, x_var = "trt", color_by = "type")

plotCounts(dds, gene=which(rownames(res_0_24) == "SalTrChr09Bg423040"), intgroup=c("type", "time", "trt"))

# ASR1:

# GSTU:
plot_gene("SalTrChr07Ag154120", dds, x_var = "trt", color_by = "type")

plotCounts(dds, gene=which(rownames(res_0_24) == "SalTrChr07Ag154120"), intgroup=c("type", "time", "trt"))

plot_gene("SalTrChr09Ag212040", dds, x_var = "trt", color_by = "type")

plotCounts(dds, gene=which(rownames(res_0_24) == "SalTrChr09Ag212040"), intgroup=c("type", "time", "trt"))

# pAPX:
plot_gene("SalTrChr01Ag020470", dds, x_var = "trt", color_by = "type")

plot_gene("SalTrChr08Ag177890", dds, x_var = "trt", color_by = "type")

plotCounts(dds, gene=which(rownames(res_0_24) == "SalTrChr08Ag177890"), intgroup=c("type", "time", "trt"))

plot_gene("SalTrChr02Bg246130", dds, x_var = "trt", color_by = "type")

plot_gene("SalTrChr08Bg390460", dds, x_var = "trt", color_by = "type")

# SDR1:
plot_gene("SalTrChr02Ag049280", dds, x_var = "trt", color_by = "type")

plot_gene("SalTrChr03Ag067930", dds, x_var = "trt", color_by = "type")

plot_gene("SalTrChr04Ag079160", dds, x_var = "trt", color_by = "type")

plotCounts(dds, gene=which(rownames(res_0_24) == "SalTrChr04Ag079160"), intgroup=c("type", "time", "trt"))

plot_gene("SalTrChr05Ag118290", dds, x_var = "trt", color_by = "type")

plot_gene("SalTrChr03Bg286200", dds, x_var = "trt", color_by = "type")

plot_gene("SalTrChr09Bg415010", dds, x_var = "trt", color_by = "type")

plotCounts(dds, gene=which(rownames(res_0_24) == "SalTrChr09Bg415010"), intgroup=c("type", "time", "trt"))


# USP:

# NAC:
plot_gene("SalTrChr07Ag148570", dds, x_var = "trt", color_by = "type")

plot_gene("SalTrChr01Bg241430", dds, x_var = "trt", color_by = "type")

# CAX1:
plot_gene("SalTrChr07Ag163170", dds, x_var = "time", color_by = "type")

plotCounts(dds, gene=which(rownames(res_0_24) == "SalTrChr07Ag163170"), intgroup=c("type", "time", "trt"))

############ make a boxplot and save it!

antiporter_countsplot <- plotCounts(dds, gene=which(rownames(res_0_24) == "SalTrChr07Ag163170"), intgroup=c("type", "time", "trt"))

antiporter_counts_table <- plotCounts(dds, gene=which(rownames(res_0_24) == "SalTrChr07Ag163170"), intgroup=c("type", "time", "trt"), returnData = TRUE)

antiporter_table <- data.frame(sample = metadata_contrast$names1,
                               treatment = c("0 HAT", "0 HAT", "0 HAT",
                                             "24 HAT", "24 HAT", "24 HAT"),
                               count = antiporter_counts_table$count,
                               time = antiporter_counts_table$time,
                               type = antiporter_counts_table$type,
                               trt = antiporter_counts_table$trt)

ggplot(data = antiporter_table) + 
  geom_point(mapping = aes(
    x = trt, y = count, color = type)) +
  theme_light() + ggtitle("Ca2+/H+ Echanger") +
  xlab("Hours after treatment") + ylab("Normalized gene count") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        plot.caption = element_text(hjust = 0.01),
        text = element_text(family = "Times New Roman", size = 12)) +
  labs(title = "All tissues 24 HAT", subtitle = "SalTrChr07Ag163170\nCAX1: Vacuolar Ca2+/H+ Exchanger",
       color = "Tissue type") +
  scale_color_brewer(palette = "Dark2")

# save this gene count plot

antiporter_boxplot <- ggplot(data = antiporter_table) + 
  geom_point(mapping = aes(
    x = trt, y = count, color = type)) +
  theme_light() + ggtitle("Ca2+/H+ Echanger") +
  xlab("Hours after treatment") + ylab("Normalized gene count") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        plot.caption = element_text(hjust = 0.01),
        text = element_text(family = "Times New Roman", size = 12)) +
  labs(title = "All tissues 24 HAT", subtitle = "SalTrChr07Ag163170\nCAX1: Vacuolar Ca2+/H+ Exchanger",
       color = "Tissue type") +
  scale_color_brewer(palette = "Dark2")

pdf(file = "03_output/05_DESeq2/Figure_files/All_24/CAX1_boxplot.pdf", width = 6, height = 6)

par(mfrow=c(1,1))

antiporter_boxplot

dev.off()
dev.off()

############################

plot_gene("SalTrChr06Ag141700", dds, x_var = "time", color_by = "type")

plotCounts(dds, gene=which(rownames(res_0_24) == "SalTrChr06Ag141700"), intgroup=c("type", "time", "trt"))

plot_gene("SalTrChr06Bg348910", dds, x_var = "time", color_by = "type")

plotCounts(dds, gene=which(rownames(res_0_24) == "SalTrChr06Bg348910"), intgroup=c("type", "time", "trt"))


# sAPX:
plot_gene("SalTrChr01Ag020470", dds, x_var = "trt", color_by = "type")

plotCounts(dds, gene=which(rownames(res_0_24) == "SalTrChr01Ag020470"), intgroup=c("type", "time", "trt"))

plot_gene("SalTrChr01Ag025000", dds, x_var = "trt", color_by = "type")

plotCounts(dds, gene=which(rownames(res_0_24) == "SalTrChr01Ag025000"), intgroup=c("type", "time", "trt"))

plot_gene("SalTrChr08Ag177890", dds, x_var = "trt", color_by = "type")

plot_gene("SalTrChr02Bg246130", dds, x_var = "trt", color_by = "type")

plot_gene("SalTrChr08Bg390460", dds, x_var = "trt", color_by = "type")

# PRxQ:
plot_gene("SalTrChr03Ag051460", dds, x_var = "trt", color_by = "type")

plot_gene("SalTrChr07Ag159830", dds, x_var = "trt", color_by = "type")

plot_gene("SalTrChr09Ag197210", dds, x_var = "trt", color_by = "type")

plot_gene("SalTrChr09Ag207800", dds, x_var = "trt", color_by = "type")

plot_gene("SalTrChr07Bg379770", dds, x_var = "trt", color_by = "type")


# MnSOD:
plot_gene("SalTrChr07Bg376420", dds, x_var = "trt", color_by = "type")


# SOS1:
plot_gene("SalTrChr04Ag088810", dds, x_var = "trt", color_by = "type")

plotCounts(dds, gene=which(rownames(res_0_24) == "SalTrChr04Ag088810"), intgroup=c("type", "time", "trt"))

plot_gene("SalTrChr04Ag089370", dds, x_var = "trt", color_by = "type")

plotCounts(dds, gene=which(rownames(res_0_24) == "SalTrChr04Ag089370"), intgroup=c("type", "time", "trt"))


# LEA1:
plot_gene("SalTrChr03Ag056540", dds, x_var = "trt", color_by = "type")

plot_gene("SalTrChr08Ag187250", dds, x_var = "trt", color_by = "type")

plot_gene("SalTrChr02Bg253260", dds, x_var = "trt", color_by = "type")

plot_gene("SalTrChr08Bg401630", dds, x_var = "trt", color_by = "type")

plotCounts(dds, gene=which(rownames(res_0_24) == "SalTrChr08Bg401630"), intgroup=c("type", "time", "trt"))

plot_gene("SalTrChr08Bg401670", dds, x_var = "trt", color_by = "type")

# TIP1:
plot_gene("SalTrChr02Ag027890", dds, x_var = "trt", color_by = "type")

plotCounts(dds, gene=which(rownames(res_0_24) == "SalTrChr02Ag027890"), intgroup=c("type", "time", "trt"))

plot_gene("SalTrChr02Bg244930", dds, x_var = "trt", color_by = "type")

plot_gene("SalTrChr03Ag067010", dds, x_var = "trt", color_by = "type")

plot_gene("SalTrChr08Ag176520", dds, x_var = "trt", color_by = "type")

plot_gene("SalTrChr01Bg238960", dds, x_var = "trt", color_by = "type")

plot_gene("SalTrChr02Bg252720", dds, x_var = "trt", color_by = "type")

plot_gene("SalTrChr03Bg287180", dds, x_var = "trt", color_by = "type")

plot_gene("SalTrChr09Bg429080", dds, x_var = "trt", color_by = "type")

plotCounts(dds, gene=which(rownames(res_0_24) == "SalTrChr09Bg429080"), intgroup=c("type", "time", "trt"))

# BZR1:
plot_gene("SalTrChr07Bg375530", dds, x_var = "time", color_by = "type")

plotCounts(dds, gene=which(rownames(res_0_24) == "SalTrChr07Bg375530"), intgroup=c("type", "time", "trt"))

######### Best genes of interest from above to plot in one pane:

plot_gene("SalTrChr04Ag088810", dds = dds, x_var = "trt", color_by = "type")# NHX1

plot_gene("SalTrChr01Ag027080", dds, x_var = "trt", color_by = "type")# MDHR

plot_gene("SalTrChr01Bg243720", dds, x_var = "trt", color_by = "type")# MDHR

plot_gene("SalTrChr09Bg423040", dds, x_var = "trt", color_by = "type")# V-ATPase

plot_gene("SalTrChr07Ag154120", dds, x_var = "trt", color_by = "type")# GSTU

plot_gene("SalTrChr08Ag177890", dds, x_var = "trt", color_by = "type")# pAPX

plot_gene("SalTrChr02Ag049280", dds, x_var = "trt", color_by = "type")# SDR1

pdf("03_output/05_DESeq2/Figure_files/All_24/other_GOI_1.pdf",
    width = 8, height = 8)

par(mfrow=c(2,2))

plotCounts(dds, gene=which(rownames(res_0_24) == "SalTrChr09Bg423040"), intgroup=c("type", "time", "trt"))

plotCounts(dds, gene=which(rownames(res_0_24) == "SalTrChr07Ag154120"), intgroup=c("type", "time", "trt"))

plotCounts(dds, gene=which(rownames(res_0_24) == "SalTrChr08Ag177890"), intgroup=c("type", "time", "trt"))

plotCounts(dds, gene=which(rownames(res_0_24) == "SalTrChr02Ag049280"), intgroup=c("type", "time", "trt"))

dev.off()
dev.off()

plot_gene("SalTrChr07Ag148570", dds, x_var = "trt", color_by = "type")# NAC

plot_gene("SalTrChr07Bg379770", dds, x_var = "trt", color_by = "type")# PRxQ

plot_gene("SalTrChr07Bg376420", dds, x_var = "trt", color_by = "type")#MnSOD

plot_gene("SalTrChr04Ag088810", dds, x_var = "trt", color_by = "type")# SOS1

pdf("03_output/05_DESeq2/Figure_files/All_24/other_GOI_2.pdf",
    width = 8, height = 8)

par(mfrow=c(2,2))

plotCounts(dds, gene=which(rownames(res_0_24) == "SalTrChr07Ag148570"), intgroup=c("type", "time", "trt"))

plotCounts(dds, gene=which(rownames(res_0_24) == "SalTrChr07Bg379770"), intgroup=c("type", "time", "trt"))

plotCounts(dds, gene=which(rownames(res_0_24) == "SalTrChr07Bg376420"), intgroup=c("type", "time", "trt"))

plotCounts(dds, gene=which(rownames(res_0_24) == "SalTrChr04Ag088810"), intgroup=c("type", "time", "trt"))

dev.off()
dev.off()

plot_gene("SalTrChr08Bg401630", dds, x_var = "trt", color_by = "type")# LEA1

plot_gene("SalTrChr02Bg244930", dds, x_var = "trt", color_by = "type")# TIP1

pdf("03_output/05_DESeq2/Figure_files/All_24/other_GOI_3.pdf",
    width = 8, height = 8)

par(mfrow=c(1,2))

plotCounts(dds, gene=which(rownames(res_0_24) == "SalTrChr08Bg401630"), intgroup=c("type", "time", "trt"))

plotCounts(dds, gene=which(rownames(res_0_24) == "SalTrChr02Bg244930"), intgroup=c("type", "time", "trt"))

dev.off()
dev.off()


##### VOLCANO PLOTS:###################

############ Enhanced volcano


pdf("03_output/05_DESeq2/Figure_files/All_24/Evolcano_All_24.pdf",
    width = 8, height = 8)

par(mfrow = c(1,1))

EnhancedVolcano(toptable = res_0_24, lab = rownames(res_0_24), 
                x = 'log2FoldChange',
                y = "pvalue",
                title = 'Differentially Expressed Genes 24 HAT',
                subtitle = 'A: Both tissue types and treatments',
                caption = bquote(~Log[2]~ "fold change cutoff = ±2; p-value cutoff: p < 0.05"),
                xlab = bquote(~Log[2]~ 'fold change'),
                pCutoff = 5e-2,
                FCcutoff = 2,
                selectLab = c(
                  "SalTrChr03Bg288430",
                  "SalTrChr05Ag119950",
                  "SalTrChr06Bg341440",
                  "SalTrChr09Bg416870",
                  "SalTrChr08Bg286720",
                  "SalTrChr03Ag060690",
                  "SalTrChr04Ag098380",
                  "SalTrChr08Bg395090"
                ),
                cutoffLineType = 'twodash',
                cutoffLineWidth = 0.8,
                cutoffLineCol = "grey",
                hline = c(10e-10),
                hlineCol = "black",
                hlineType = "twodash",
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


############## GETTING LISTS OF SIGNIFICANTLY CHANGING GENES #####################

# Select the significant subset of genes that are up-regulated
Up_24 <- subset(res_0_24, padj < 0.1 & log2FoldChange > 0.5)
Up_24 <- Up_24[order(Up_24$padj),] #order them

head(Up_24) # Check them
dim(Up_24)

# Select the significant subset of genes that are down-regulated
Down_24 <- subset(res_0_24, padj < 0.1 & log2FoldChange < -0.5)
Down_24 <- Down_24[order(Down_24$padj),]

head(Down_24)
dim(Down_24)

# Save these lists to output files:
write(rownames(Up_24), file = "03_output/05_DESeq2/Up_all_0_24.txt", sep = "\n")
write(rownames(Down_24), file = "03_output/05_DESeq2/Down_all_0_24.txt", sep = "\n")
write(rownames(res_0_24), file = "03_output/05_DESeq2/All_all_0_24.txt", sep = "\n")

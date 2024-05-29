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
library(extrafont)
#library(tidyverse) # causes some minor conflicts, but you can install this later for creating plots


font_import()
y
loadfonts()


#################################################

# Set the working directory:

setwd("~/SalTrGenome/RNA_Seq/RT_StressTrials/RNA-seq_data")

# Import the counts data file

countsData <- read.table(file = "03_output/03_feature/STAR_counts_v6.txt", header = FALSE, row.names = 1, skip = 2) # 

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
# In the instruction for DESeq2, it states: "We read in a count matrix,
# which we will name cts, and the sample information table, which we will name coldata."

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

cts_contrast <- cts[,c(1:3, 13:15, 25:27)]

# Now instead of indexing the metadata accordingly I can import the metadata file I need for this run that was generated in Unix:

metadata_contrast <- metadata_full[c(1:3, 13:15, 25:27),]

# Now I'll need to rename the column headers again for this new metadata object:

colnames(metadata_contrast) <- c("fasta1", "fasta2", "names1", "names2", "type", "trt", "time", "rep")

# Next we need to make an information called coltable. We can make this out of the metadata table.
# Reorganize the metadata table so the names2 column are now row headers

rownames(metadata_contrast)<- metadata_contrast$names1

# Adjust the following based on the contrasts included if needed:

coldata <- metadata_contrast[,c("time", "rep", "trt")]
coldata$time <- as.factor(coldata$time)
coldata$trt <- as.factor((coldata$trt))
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
                              design = ~ time)


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

dds$time <- factor(dds$time, levels = c("0","6","24"))
dds$trt <- factor(dds$trt, levels = c("Contl", "Cold"))

# PERFORM DESEQ2 analysis:

dds <- DESeq(dds)

plotDispEsts(dds)

pdf("03_output/05_DESeq2/Figure_files/Cold_shoot/dispEst.pdf", height = 6, width = 8)
par(mfrow=c(1,1))

plotDispEsts(dds)

dev.off()
dev.off()

# Other interesting QC figures:
# https://cran.r-project.org/web/packages/RNAseqQC/vignettes/introduction.html
plot_total_counts(dds)
plot_library_complexity(dds)

pdf("03_output/05_DESeq2/Figure_files/Cold_shoot/total_counts.pdf", height = 7, width = 7)
par(mfrow=c(1,1))

plot_total_counts(dds)

dev.off()
dev.off()

pdf("03_output/05_DESeq2/Figure_files/Cold_shoot/library_complx.pdf", height = 7, width = 7)
par(mfrow=c(1,1))

plot_library_complexity(dds)

dev.off()
dev.off()


dds$sizeFactor

# Both of these datasets are good supplemental tables for papers

head(counts(dds, normalized = TRUE))
head(counts(dds, normalized = FALSE))

counts_norm <- counts(dds, normalized = TRUE)
counts_nonNorm <- counts(dds, normalized = FALSE)

write.table(counts_norm, file = "03_output/05_DESeq2/03_cts_files/Cold_shoot_cts_norm.txt", sep = "\t",
      row.names = TRUE, col.names = TRUE)
write.table(counts_nonNorm, file = "03_output/05_DESeq2/03_cts_files/Cold_shoot_cts_nonNorm.txt", sep = "\t",
            row.names = TRUE, col.names = TRUE)

# make the summary file for reference:

counts_sum <- summary(counts(dds, normalized = TRUE))
write.table(counts_sum, file = "03_output/05_DESeq2/03_cts_files/Cold_shoot_cts_norm_summary.txt", sep = "\t",
            row.names = TRUE, col.names = TRUE)

#' upload both of these tables for supplimental information in the manuscript

############## PERFORM LOG FOLD CHANGE SHRINKAGE FOR VISUALIZATION #####################
#######PERFORM LOG FOLD CHANGE SHRINKAGE FOR VISUALIZATION

# An input requirement of the lfcShrink function is a coef term. This is pulled from the resultsNames of dds:
# An additional method that is appropriate in some datasets is to perform lfcshrinkage.
# This uses modeling to adjust the values of normalized read numbers.
# This helps improve plotting and prepares the data for some statsitcsal analysis. This creates "homoskedactic" data. 

#resultsNames(dds)

#res_lfc <- lfcShrink(dds, coef="type_Root_vs_Shoot", res = res)

#summary(res_lfc)

#' the above did not work for my specific case for some reason. I followed
#' this work through: https://www.biostars.org/p/448959/#484944
#' originally I did not have output from resultNames(res) so I used this work around
#' here is the edited scripting to get lfcShrink to work:

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

resLFC_0_6 <- lfcShrink(dds, coef = 2, type = "apeglm")

resLFC_0_24 <- lfcShrink(dds, coef = 3, type = "apeglm")

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

plot_pca(vsd, PC_x = 1, PC_y = 2, color_by = "rep", shape_by = "time")

pdf("03_output/05_DESeq2/Figure_files/Cold_shoot/PCA.pdf", height = 8, width = 7)
par(mfrow=c(1,1))

plot_pca(vsd, PC_x = 1, PC_y = 2, color_by = "rep", shape_by = "time")

dev.off()
dev.off() 

# create a table with the loadings for each PC

pca <- plot_pca(vsd, PC_x = 1, PC_y = 2, color_by = "rep", shape_by = "time")

loadings <- as.data.frame(pca$loadings)

#export loadings into a file
write.csv(x = loadings, file = "03_output/05_DESeq2/Figure_files/Cold_shoot/PCA_loadings.csv")

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


pdf("03_output/05_DESeq2/Figure_files/Cold_shoot/heatmap.pdf", height = 6, width = 7)
par(mfrow=c(1,1))
p

dev.off()
dev.off() 

############## DIFFERENTIAL EXPRESSION ANALYSIS #####################

# calculate the statistically significant differences
# Set your contrast expression: Contrast - a character vector with exactly three elements: 
# the name of a factor in the design formula, the name of the numerator level for the fold 
# change, and the name of the denominator level for the fold change (simplest case)
# Set your lfcThreshold to a value of your choosing. This will set whether you want to 
# exclude any genes as significant because their fold change is too low. Default is 0. 
# That means, take all statistically significant genes. In this dataset, I selected 0.5 
# because there were A LOT of significant genes and I'd like to reduce the number.

# Set your alpha value. This is your p-value cut-off for significance. The default is 0.1. 

res_0_6 <- results(dds,
                            lfcThreshold = 0.5,
                            alpha = 0.05,
                            contrast=c("time", "6", "0"))

res_0_24 <- results(dds,
                   lfcThreshold = 0.5,
                   alpha = 0.05,
                   contrast=c("time", "24", "0"))

summary(res_0_6)

summary(res_0_24)

#' these values are the default for DESeq. you can change them and check the summary to see how it was
#' affected

# If theres only 1 condition term and you want to use the default settings, this will work...
#results_table_default <- results(dds)

##### MA PLOTS:

# Plot the the default MA-plot

# edit what you want:

par(mfrow=c(2,2), mar=c(3,3,3,3))
plotMA(res_0_6, main = "0 versus 6 HAT \nunshrunken", ylim = c(-7,7), 
       ylab = "LFC (6 HAT/ 0 HAT)",
       xlab = "means of normalized counts")

plotMA(resLFC_0_6, main="0 versus 6 HAT\nshrunken", ylim = c(-7,7), 
       ylab = "LFC (6 HAT / 0 HAT)",
       xlab = "means of normalized counts")

plotMA(res_0_24, main="0 versus 24 HAT\nunshrunken", ylim = c(-7,7), 
       ylab = "LFC (24 HAT / 0 HAT)",
       xlab = "means of normalized counts")

plotMA(resLFC_0_24, main="0 versus 24 HAT\nshrunken", ylim = c(-7,7), 
       ylab = "LFC (24 HAT / 0 HAT)",
       xlab = "means of normalized counts")

# Save these plots in one panel:

pdf("03_output/05_DESeq2/Figure_files/Cold_shoot/MAplots.pdf", height = 6, width = 8)

par(mfrow=c(2,2), mar=c(3,3,3,3))

plotMA(res_0_6, main = "0 versus 6 HAT \nunshrunken", ylim = c(-7,7), 
       ylab = "LFC (6 HAT/ 0 HAT)",
       xlab = "means of normalized counts")

plotMA(resLFC_0_6, main="0 versus 6 HAT\nshrunken", ylim = c(-7,7), 
       ylab = "LFC (6 HAT / 0 HAT)",
       xlab = "means of normalized counts")

plotMA(res_0_24, main="0 versus 24 HAT\nunshrunken", ylim = c(-7,7), 
       ylab = "LFC (24 HAT / 0 HAT)",
       xlab = "means of normalized counts")

plotMA(resLFC_0_24, main="0 versus 24 HAT\nshrunken", ylim = c(-7,7), 
       ylab = "LFC (24 HAT / 0 HAT)",
       xlab = "means of normalized counts")

dev.off()
dev.off() 

##################  Exploring and exporting results ##################  

##### KNOWN GENES:

# Check known genes to make sure everything is working as predicted. Check out :WBGene00018393 # aka msra-1
# EXAM PROJECT INSTRUCTIONS - This is optional - can you plot a gene of interest? 

#plotCounts(dds, gene=which(rownames(res_EcolVBsubt) == "WBGene00018393"), intgroup="treatment")
#plotCounts(dds, gene=which(rownames(res_EcolVBsubt) == "WBGene00018393"), intgroup=c("treatment", "time"))

#' this is where I can look at specific genes that I am interested in in other salsola or 
#' mesohalophyte studies and see how they are expressed in my system.

## Enrichment: Look up these genes on the Genome Browser: http://genome.ucsc.edu/s/Erin%20Osborne/GomezOrte_191210
## --> Search for msra-1 (Which is WBGene00018393)
########## NHX1
plotCounts(dds, gene=which(rownames(res_0_6) == "SalTrChr04Ag088810"), intgroup=c("trt", "time"))

antiporter_countsplot <- plotCounts(dds, gene=which(rownames(res_0_6) == "SalTrChr04Ag088810"), intgroup=c("trt", "time"))

antiporter_counts_table <- plotCounts(dds, gene=which(rownames(res_0_6) == "SalTrChr04Ag088810"), intgroup=c("trt", "time"), returnData = TRUE)

#' this is where I can look at specific genes that I am interested in in other salsola or 
#' mesohalophyte studies and see how they are expressed in my system.

## Enrichment: Look up these genes on the Genome Browser: http://genome.ucsc.edu/s/Erin%20Osborne/GomezOrte_191210
## --> Search for msra-1 (Which is WBGene00018393)

antiporter_table <- data.frame(sample = metadata_contrast$names1,
                               treatment = c("0 HAT", "0 HAT", "0 HAT",
                                             "6 HAT", "6 HAT", "6 HAT",
                                             "24 HAT", "24 HAT", "24 HAT"),
                               count = antiporter_counts_table$count,
                               time = antiporter_counts_table$time)

ggplot(data = antiporter_table, mapping = aes(
  x = time, y = count)) + 
  geom_point() + theme_light() + 
  xlab("Hours after treatment") + ylab("Normalized gene count") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        plot.caption = element_text(hjust = 0.01),
        legend.position = 'none',
        text = element_text(family = "Times New Roman", size = 12)) +
  labs(title = "Cold treated shoot tissue", subtitle = "SalTrChr04Ag088810\nNHX1: Vacuolar Na+/H+ Exchanger",
       caption = "Tukey adjusted pairwise comparison: 0 - 24 HAT: p < 0.001; 0 - 6 HAT: p = 0.046")

# save this gene count plot

antiporter_boxplot <- ggplot(data = antiporter_table, mapping = aes(
  x = time, y = count)) + 
  geom_point() + theme_light() + 
  xlab("Hours after treatment") + ylab("Normalized gene count") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        plot.caption = element_text(hjust = 0.01),
        legend.position = 'none',
        text = element_text(family = "Times New Roman", size = 12)) +
  labs(title = "Cold treated shoot tissue", subtitle = "SalTrChr04Ag088810\nNHX1: Vacuolar Na+/H+ Exchanger",
       caption = "Tukey adjusted pairwise comparison: 0 - 24 HAT: p < 0.001; 0 - 6 HAT: p = 0.046")


pdf(file = "03_output/05_DESeq2/Figure_files/Cold_shoot/NHX1_boxplot.pdf", width = 6, height = 6)

par(mfrow=c(1,1))

antiporter_boxplot

dev.off()
dev.off()

# Pairwise comparison for each timepoint

antiporter_anova <- anova(lm(count ~ treatment, data = antiporter_table))

antiporter_emmeans <- emmeans(lm(count ~ treatment, data = antiporter_table),
                              ~ treatment)

antiporter_pairs <- pairs(
  antiporter_emmeans, adjust = "Tukey"
)

kbl(antiporter_pairs, 
    caption = "NA+/H+ Echanger ANOVA", 
    booktabs = T,
    digits = 3,
) %>%
  kable_styling(latex_options = c("striped", "hold_position"))

############### CAX1
plotCounts(dds, gene=which(rownames(res_0_6) == "SalTrChr07Ag163170"), intgroup=c("trt", "time"))

antiporter_countsplot <- plotCounts(dds, gene=which(rownames(res_0_6) == "SalTrChr07Ag163170"), intgroup=c("trt", "time"))

antiporter_counts_table <- plotCounts(dds, gene=which(rownames(res_0_6) == "SalTrChr07Ag163170"), intgroup=c("trt", "time"), returnData = TRUE)

#' this is where I can look at specific genes that I am interested in in other salsola or 
#' mesohalophyte studies and see how they are expressed in my system.

## Enrichment: Look up these genes on the Genome Browser: http://genome.ucsc.edu/s/Erin%20Osborne/GomezOrte_191210
## --> Search for msra-1 (Which is WBGene00018393)

antiporter_table <- data.frame(sample = metadata_contrast$names1,
                               treatment = c("0 HAT", "0 HAT", "0 HAT",
                                             "6 HAT", "6 HAT", "6 HAT",
                                             "24 HAT", "24 HAT", "24 HAT"),
                               count = antiporter_counts_table$count,
                               time = antiporter_counts_table$time)

ggplot(data = antiporter_table, mapping = aes(
  x = time, y = count)) + 
  geom_point() + theme_light() +
  xlab("Hours after treatment") + ylab("Normalized gene count") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        plot.caption = element_text(hjust = 0.01),
        text = element_text(family = "Times New Roman", size = 12)) +
  labs(title = "Cold treated shoot tissue", subtitle = "SalTrChr07Ag163170\nCAX1: Vacuolar Ca2+/H+ Exchanger",
       caption = "Tukey adjusted pairwise comparison: 0 - 24 HAT: p = 0.136; 0 - 6 HAT: p = 0.695")

# save this gene count plot

antiporter_boxplot <- ggplot(data = antiporter_table, mapping = aes(
  x = time, y = count)) + 
  geom_point() + theme_light() +
  xlab("Hours after treatment") + ylab("Normalized gene count") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        plot.caption = element_text(hjust = 0.01),
        text = element_text(family = "Times New Roman", size = 12)) +
  labs(title = "Cold treated shoot tissue", subtitle = "SalTrChr07Ag163170\nCAX1: Vacuolar Ca2+/H+ Exchanger",
       caption = "Tukey adjusted pairwise comparison: 0 - 24 HAT: p = 0.136; 0 - 6 HAT: p = 0.695")


pdf(file = "03_output/05_DESeq2/Figure_files/Cold_shoot/CAX1_boxplot.pdf", width = 6, height = 6)

par(mfrow=c(1,1))

antiporter_boxplot

dev.off()
dev.off()

# Pairwise comparison for each timepoint

antiporter_anova <- anova(lm(count ~ treatment, data = antiporter_table))

antiporter_emmeans <- emmeans(lm(count ~ treatment, data = antiporter_table),
                              ~ treatment)

antiporter_pairs <- pairs(
  antiporter_emmeans, adjust = "Tukey"
)

kbl(antiporter_pairs, 
    caption = "NA+/H+ Echanger ANOVA", 
    booktabs = T,
    digits = 3,
) %>%
  kable_styling(latex_options = c("striped", "hold_position"))

####################

#vignette notes from:
# https://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#exploring-and-exporting-results

# This will pull the gene with the lowest adjusted p-value instead of having to
# find the gene itself. Cool. I should do what I did above with these..

plotCounts(dds, gene=which.min(res_0_6$padj), intgroup="time")

plotCounts(dds, gene=which.min(res_0_24$padj), intgroup="time")

# A different way to plot genes from the vignette using the dds object:

plot_gene("SalTrChr04Ag088810", dds, x_var = "time", color_by = "rep")

pdf(file = "03_output/05_DESeq2/Figure_files/Cold_shoot/antiporter_dotplot.pdf", width = 6, height = 6)

par(mfrow=c(1,1))

plot_gene("SalTrChr04Ag088810", dds, x_var = "time", color_by = "rep")

dev.off()
dev.off()

############## GETTING LISTS OF SIGNIFICANTLY CHANGING GENES #####################

# Select the significant subset of genes that are up-regulated
Up_0_6 <- subset(res_0_6, padj < 0.1 & log2FoldChange > 0.5)
Up_0_6 <- Up_0_6[order(Up_0_6$padj),] #order them

head(Up_0_6) # Check them
dim(Up_0_6)

# Select the significant subset of genes that are down-regulated
Down_0_6 <- subset(res_0_6, padj < 0.1 & log2FoldChange < -0.5)
Down_0_6 <- Down_0_6[order(Down_0_6$padj),]

head(Down_0_6)
dim(Down_0_6)

# Save these lists to output files:
write(rownames(Up_0_6), file = "03_output/05_DESeq2/Up_shoot_Cold_0_6.txt", sep = "\n")
write(rownames(Down_0_6), file = "03_output/05_DESeq2/Down_shoot_Cold_0_6.txt", sep = "\n")
write(rownames(res_0_6), file = "03_output/05_DESeq2/All_shoot_Cold_0_6.txt", sep = "\n")

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
write(rownames(Up_0_24), file = "03_output/05_DESeq2/Up_shoot_Cold_0_24.txt", sep = "\n")
write(rownames(Down_0_24), file = "03_output/05_DESeq2/Down_shoot_Cold_0_24.txt", sep = "\n")
write(rownames(res_0_24), file = "03_output/05_DESeq2/All_shoot_Cold_0_24.txt", sep = "\n")

#' these are the lists of genes that you'll use for the GO websites. Paste in your
#' query list and your full list as the background list onto these websites. 

# Fun thing to do. See what KEGG pathways or GO Ontology terms are associated with your different lists of genes.
# Go to DAVID: https://david.ncifcrf.gov/
#      Click on "Start Analysis
#      Copy and paste your WBGENEID list into the "Upload" tab and "A: Paste a list" field
#      Under "Step 2: Select Identifier" select "WORMBASE_GENE_ID"
#      Under "Step 3:" select "Gene List"
#      Submit

############## MAKE HIERARCHICALLY CLUSTERED HEATMAPS OF ALL CHANGING GENES #####################

# Remember we assessed statistically significantly changing genes as...

#res_contrast_2 <- results(dds,
#                         lfcThreshold = 0.5,
#                          alpha = 0.1,
#                          contrast=c("", "", ""))

# Next, let's see what genes differ between factor levels included in the contrast:
# get these using resultsNames(dds)

res_contrast1_2 <- results(dds,
                       lfcThreshold = 0.5,
                       alpha = 0.1,
                       contrast=c("time", "0", "6"))

res_contrast3_2 <- results(dds,
                       lfcThreshold = 0.5,
                       alpha = 0.1,
                       contrast=c("time", "6", "24"))

#Subset each results table for just the differentially expressed genes:
# I'm making this a little more stringent 
#sign_contrast <- subset(res_contrast_2, padj < 0.1)
#dim(subset(res_contrast_2, padj < 0.1))

sign_contrast1 <- subset(res_contrast1_2, padj < 0.1)
sign_contrast3 <- subset(res_contrast3_2, padj < 0.1)

#Determine how many genes were captured and merge them:
changing_genes <- rbind(sign_contrast1, sign_contrast3)

dim(changing_genes)
length(unique(rownames(changing_genes)))
# Get all r-stabilized log values of all the data:
rdl_all <- assay(rlog(dds, blind=FALSE))

#Subset just the changing genes:
changing_lrt_rdl <- subset(rdl_all, rownames(rdl_all) %in% rownames(changing_genes))
dim(changing_lrt_rdl)
head(changing_lrt_rdl)

# Make sure it is in matrix form:
class(changing_lrt_rdl)

# Draw a heat map
# Scale by row, listed below as scale = "row". 
# This is really important. It sets the mean of every row to 0 and the standard deviation of every row to 1:
p <- pheatmap(changing_lrt_rdl, 
              scale="row", 
              color = colorRampPalette(c("blue", "white", "red"), space = "Lab")(100),
              cluster_rows=TRUE, 
              cluster_cols=TRUE, 
              clustering_distance_rows = "euclidean", 
              clustering_method = "complete",
              show_rownames = FALSE)

p
help(pheatmap)
# Tip: If euclidean doesn't look good, try "correlation" other clustering distances
# Tip: If complete doesn't look good, try other methods
# Save the clustered heatmap plot as a .pdf:

pdf("03_output/05_DESeq2/Figure_files/Cold_shoot/clustered_heatmap.pdf", width = 6, height = 8)

p

dev.off()
dev.off()





############## ADVANCED OPTIONAL TOPICS
#' ok I want PCA plots and I can't find it in the vignette..
#' I think I have to install the following package:

#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("ReportingTools")

#BiocManager::install("pcaExplorer")

library(pcaExplorer)

#BiocManager::install('airway')

pcaExplorer(dds = dds)

# This looks really cool but It'll take some time to figure out. Looks like I
# can feed it an annotation file so that's awesome. 



# Get versions
sessionInfo()

# Please cite R in your publications when using it
citation()


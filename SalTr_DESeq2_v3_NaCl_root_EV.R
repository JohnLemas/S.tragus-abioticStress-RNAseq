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
#library(tidyverse) # causes some minor conflicts, but you can install this later for creating plots

#################################################

# Set the working directory:

setwd("~/SalTrGenome/RNA_Seq/RT_StressTrials/RNA-seq_data")

# Import the counts data file

countsData <- read.table(file = "03_output/03_feature/STAR_counts_v4.txt", header = FALSE, row.names = 1, skip = 2) # 

# Read in the metadata
metadata_full <- read.table(file = "01_input/S.tragus_metadata_v3.txt", header = FALSE)



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

cts_contrast <- cts[,c(4:6,16:18,28:30,40:42)]

# Now instead of indexing the metadata accordingly I can import the metadata file I need for this run that was generated in Unix:

metadata_contrast <- read.table("01_input/S.tragus_metadata_v3_root_NaCl.txt", header = FALSE)

# Now I'll need to rename the column headers again for this new metadata object:

colnames(metadata_contrast) <- c("fasta1", "fasta2", "names1", "names2", "type", "trt", "time", "rep")

# Next we need to make an information called coltable. We can make this out of the metadata table.
# Reorganize the metadata table so the names2 column are now row headers

rownames(metadata_contrast)<- metadata_contrast$names1

# Adjust the following based on the contrasts included if needed:

coldata <- metadata_contrast[,c("time", "rep")]
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

dds$time <- factor(dds$time, levels = c("0","3", "8", "24"))

# PERFORM DESEQ2 analysis:

dds <- DESeq(dds)

dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)
dds <- nbinomWaldTest(dds)

relevel(dds$time, ref = "0")

# Now when I run resultNames(dds) the reference level is defined as above

resultsNames(dds)

# Now we need to re-reun the Wald test:

dds <- nbinomWaldTest(dds)
resultsNames(dds)

# ok now we can run lfcShrink(dds = dds, coef = [2, 3, or 4], type = "ageglm")
# coef = 2 means 0 vs 3; coef = 3 means 0 vs 8; coef = 4 means 0 vs 24

resLFC_0_3 <- lfcShrink(dds, coef = 2, type = "apeglm")

resLFC_0_8 <- lfcShrink(dds, coef = 3, type = "apeglm")

resLFC_0_24 <- lfcShrink(dds, coef = 4, type = "apeglm")

res_0_3 <- results(dds,
                   lfcThreshold = 0.5,
                   alpha = 0.05,
                   contrast=c("time", "3", "0"))

res_0_8 <- results(dds,
                   lfcThreshold = 0.5,
                   alpha = 0.05,
                   contrast=c("time", "8", "0"))

res_0_24 <- results(dds,
                    lfcThreshold = 0.5,
                    alpha = 0.05,
                    contrast=c("time", "24", "0"))

############## Enhanced Volcano plots:
# 3 HAT

pdf("03_output/05_DESeq2/Figure_files/NaCl_root/Evolcano_0_3.pdf", 
    width = 8, height = 8)

par(mfrow = c(1,1))

EnhancedVolcano(toptable = res_0_3, lab = rownames(res_0_3), 
                x = 'log2FoldChange',
                y = "pvalue",
                title = 'Differentially Expressed Genes',
                subtitle = 'D. Salt treated root tissue 3 HAT',
                caption = bquote(~Log[2]~ "fold change cutoff = ±2; p-value cutoff = p < 0.05"),
                xlab = bquote(~Log[2]~ 'fold change'),
                selectLab = c('SalTrChr07Bg363570',
                              'SalTrChr09Ag196440',
                              'SalTrChr06Ag138210',
                              'SalTrChr04Bg308350',
                              'SalTrChr07Bg364500'),
                pCutoff = 5e-2,
                FCcutoff = 2,
                cutoffLineType = 'twodash',
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

# 8 HAT

pdf("03_output/05_DESeq2/Figure_files/NaCl_root/Evolcano_0_8.pdf",
    width = 8, height = 8)

par(mfrow = c(1,1))

EnhancedVolcano(toptable = res_0_8, lab = rownames(res_0_8), 
                x = 'log2FoldChange',
                y = "pvalue",
                title = 'Differentially Expressed Genes',
                subtitle = 'E. Salt treated root tissue 8 HAT',
                caption = bquote(~Log[2]~ "fold change cutoff = ±2; p-value cutoff = p < 0.05"),
                xlab = bquote(~Log[2]~ 'fold change'),
                selectLab = c('SalTrChr04Bg313440',
                              'SalTrChr04Ag097900',
                              'SalTrChr09Bg431380',
                              'SalTrChr06Bg349590',
                              'SalTrChr09Ag215460',
                              'SalTrChr06Ag141140',
                              'SalTrChr08Ag184120',
                              'SalTrChr07Ag146460',
                              'SalTrChr06Bg354650',
                              'SalTrChr05Bg336460',
                              'SalTrChr06Ag129910',
                              'SalTrChr06Bg350180',
                              'SalTrChr05Ag121400'),
                pCutoff = 5e-2,
                FCcutoff = 2,
                cutoffLineType = 'twodash',
                cutoffLineCol = 'grey',
                hline = 10e-30,
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

# 24 HAT

pdf("03_output/05_DESeq2/Figure_files/NaCl_root/Evolcano_0_24.pdf",
    width = 8, height = 8)

par(mfrow = c(1,1))

EnhancedVolcano(toptable = res_0_24, lab = rownames(res_0_24), 
                x = 'log2FoldChange',
                y = "pvalue",
                title = 'Differentially Expressed Genes',
                subtitle = 'F. Salt treated root tissue 24 HAT',
                caption = bquote(~Log[2]~ "fold change cutoff = ±2; p-value cutoff = p < 0.05"),
                xlab = bquote(~Log[2]~ 'fold change'),
                selectLab = c('SalTrChr07Bg363570',
                              'SalTrChr01Bg222620',
                              'SalTrChr07Ag165940',
                              'SalTrChr03Bg286720',
                              'SalTrChr02Ag032720',
                              'SalTrChr05Ag118930',
                              'SalTrChr04Ag079170',
                              'SalTrChr01Ag010940',
                              'SalTrChr05Ag118100'),
                pCutoff = 5e-2,
                FCcutoff = 2,
                cutoffLineType = 'twodash',
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

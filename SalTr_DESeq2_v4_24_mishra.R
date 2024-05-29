#######################################

#__Date:17 May 2024
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
library(patchwork)
library(extrafont)

# load fonts

font_import()
y
loadfonts()


##################################

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

countColNum <- dim(countsData)[2]

cts <- as.matrix(countsData[,6:countColNum])

cts_contrast <- cts[,c(1:6,25:36)]

metadata_contrast <- metadata_full[c(1:6,25:36),]

rownames(metadata_contrast) <- metadata_contrast$names1

coldata <- metadata_contrast[,c("type", "trt", "time", "rep")]
coldata$type <- as.factor(coldata$type)
coldata$trt <- as.factor(coldata$trt)
coldata$time <- as.factor(coldata$time)
coldata$rep <- as.factor(coldata$rep)
str(coldata)

all(rownames(coldata) == colnames(cts_contrast))

dds <- DESeqDataSetFromMatrix(countData = cts_contrast,
                              colData = coldata,
                              design = ~ time + type)

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

dds$type <- factor(dds$type, levels = c("Shoot","Root"))
dds$time <- factor(dds$time, levels = c("0", "24"))
dds$trt <- factor(dds$trt, levels = c("Contl", "NaCl", "Cold"))

dds <- DESeq(dds)

dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)
dds <- nbinomWaldTest(dds)

relevel(dds$time, ref = "0")
relevel(dds$trt, ref = "Contl")
dds <- nbinomWaldTest(dds)
resultsNames(dds)
resLFC_0_24 <- lfcShrink(dds, coef = 2, type = "apeglm")
res_0_24 <- results(dds,
                    lfcThreshold = 0.5,
                    alpha = 0.05,
                    contrast=c("time", "24", "0"))

#####################################################################################################
######### Best genes of interest from Mishra et al to plot in one pane:

plot_gene("SalTrChr09Bg423040", dds, x_var = "trt", color_by = "type")# V-ATPase

ATPase_counts_table <- plotCounts(dds, gene=which(rownames(res_0_24) == "SalTrChr09Bg423040"), intgroup=c("type", "time", "trt"), returnData = TRUE)

ATPase_table <- data.frame(sample = metadata_contrast$names1,
                           treatment = c("0 HAT", "0 HAT", "0 HAT",
                           "24 HAT", "24 HAT", "24 HAT"),
                            count =ATPase_counts_table$count,
                            time = ATPase_counts_table$time,
                            type = ATPase_counts_table$type,
                            trt = ATPase_counts_table$trt)

###############

plot_gene("SalTrChr07Ag154120", dds, x_var = "trt", color_by = "type")# GSTU

GSTU_counts_table <- plotCounts(dds, gene=which(rownames(res_0_24) == "SalTrChr07Ag154120"), intgroup=c("type", "time", "trt"), returnData = TRUE)

GSTU_table <- data.frame(sample = metadata_contrast$names1,
                         treatment = c("0 HAT", "0 HAT", "0 HAT",
                         "24 HAT", "24 HAT", "24 HAT"),
                         count =GSTU_counts_table$count,
                         time = GSTU_counts_table$time,
                         type = GSTU_counts_table$type,
                         trt = GSTU_counts_table$trt)

#############

plot_gene("SalTrChr08Ag177890", dds, x_var = "trt", color_by = "type")# pAPX

pAPX_counts_table <- plotCounts(dds, gene=which(rownames(res_0_24) == "SalTrChr08Ag177890"), intgroup=c("type", "time", "trt"), returnData = TRUE)

pAPX_table <- data.frame(sample = metadata_contrast$names1,
                         treatment = c("0 HAT", "0 HAT", "0 HAT",
                                       "24 HAT", "24 HAT", "24 HAT"),
                         count =pAPX_counts_table$count,
                         time = pAPX_counts_table$time,
                         type = pAPX_counts_table$type,
                         trt = pAPX_counts_table$trt)

################

plot_gene("SalTrChr02Ag049280", dds, x_var = "trt", color_by = "type")# SDR1

SDR1_counts_table <- plotCounts(dds, gene=which(rownames(res_0_24) == "SalTrChr02Ag049280"), intgroup=c("type", "time", "trt"), returnData = TRUE)

SDR1_table <- data.frame(sample = metadata_contrast$names1,
                         treatment = c("0 HAT", "0 HAT", "0 HAT",
                                       "24 HAT", "24 HAT", "24 HAT"),
                         count =SDR1_counts_table$count,
                         time = SDR1_counts_table$time,
                         type = SDR1_counts_table$type,
                         trt = SDR1_counts_table$trt)

##############

plot_gene("SalTrChr07Ag148570", dds, x_var = "trt", color_by = "type")# NAC

NAC_counts_table <- plotCounts(dds, gene=which(rownames(res_0_24) == "SalTrChr07Ag148570"), intgroup=c("type", "time", "trt"), returnData = TRUE)

NAC_table <- data.frame(sample = metadata_contrast$names1,
                        treatment = c("0 HAT", "0 HAT", "0 HAT",
                                      "24 HAT", "24 HAT", "24 HAT"),
                        count =NAC_counts_table$count,
                        time = NAC_counts_table$time,
                        type = NAC_counts_table$type,
                        trt = NAC_counts_table$trt)
    
############

plot_gene("SalTrChr07Bg379770", dds, x_var = "trt", color_by = "type")# PRxQ

PRxQ_counts_table <- plotCounts(dds, gene=which(rownames(res_0_24) == "SalTrChr07Bg379770"), intgroup=c("type", "time", "trt"), returnData = TRUE)

PRxQ_table <- data.frame(sample = metadata_contrast$names1,
                         treatment = c("0 HAT", "0 HAT", "0 HAT",
                                       "24 HAT", "24 HAT", "24 HAT"),
                         count =PRxQ_counts_table$count,
                         time = PRxQ_counts_table$time,
                         type = PRxQ_counts_table$type,
                         trt = PRxQ_counts_table$trt)
 
 ############
 
plot_gene("SalTrChr07Bg376420", dds, x_var = "trt", color_by = "type")#MnSOD

MnSOD_counts_table <- plotCounts(dds, gene=which(rownames(res_0_24) == "SalTrChr07Bg376420"), intgroup=c("type", "time", "trt"), returnData = TRUE)

MnSOD_table <- data.frame(sample = metadata_contrast$names1,
                          treatment = c("0 HAT", "0 HAT", "0 HAT",
                                        "24 HAT", "24 HAT", "24 HAT"),
                          count =MnSOD_counts_table$count,
                          time = MnSOD_counts_table$time,
                          type = MnSOD_counts_table$type,
                          trt = MnSOD_counts_table$trt)

###########
        
plot_gene("SalTrChr04Ag088810", dds, x_var = "trt", color_by = "type")# SOS1

SOS1_counts_table <- plotCounts(dds, gene=which(rownames(res_0_24) == "SalTrChr04Ag088810"), intgroup=c("type", "time", "trt"), returnData = TRUE)

SOS1_table <- data.frame(sample = metadata_contrast$names1,
                                 treatment = c("0 HAT", "0 HAT", "0 HAT",
                                               "24 HAT", "24 HAT", "24 HAT"),
                                 count =SOS1_counts_table$count,
                                 time = SOS1_counts_table$time,
                                 type = SOS1_counts_table$type,
                                 trt = SOS1_counts_table$trt)

###########

plot_gene("SalTrChr08Bg401630", dds, x_var = "trt", color_by = "type")# LEA1

LEA1_counts_table <- plotCounts(dds, gene=which(rownames(res_0_24) == "SalTrChr08Bg401630"), intgroup=c("type", "time", "trt"), returnData = TRUE)

LEA1_table <- data.frame(sample = metadata_contrast$names1,
                         treatment = c("0 HAT", "0 HAT", "0 HAT",
                                       "24 HAT", "24 HAT", "24 HAT"),
                         count =LEA1_counts_table$count,
                         time = LEA1_counts_table$time,
                         type = LEA1_counts_table$type,
                         trt = LEA1_counts_table$trt)

###########

plot_gene("SalTrChr02Bg244930", dds, x_var = "trt", color_by = "type")# TIP1

TIP1_counts_table <- plotCounts(dds, gene=which(rownames(res_0_24) == "SalTrChr02Bg244930"), intgroup=c("type", "time", "trt"), returnData = TRUE)

TIP1_table <- data.frame(sample = metadata_contrast$names1,
                         treatment = c("0 HAT", "0 HAT", "0 HAT",
                                       "24 HAT", "24 HAT", "24 HAT"),
                         count =TIP1_counts_table$count,
                         time = TIP1_counts_table$time,
                         type = TIP1_counts_table$type,
                         trt = TIP1_counts_table$trt)

############# build ggplots

p1 <- ggplot(data = ATPase_table) +
  geom_point(mapping = aes(
    x = trt, y = count, color = type)) +
  theme_light() +
  xlab("") + ylab("Normalized gene count") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        plot.caption = element_text(hjust = 0.01),
        text = element_text(family = "Times New Roman", size = 12)) +
  labs(subtitle = "Vacuolar ATPase\nSalTrChr09Bg423040", color = "Tissue type") +
  scale_color_brewer(palette = "Dark2") +
  guides(color = "none")

p2 <- ggplot(data = GSTU_table) +
  geom_point(mapping = aes(
    x = trt, y = count, color = type)) +
  theme_light() +
  xlab("") + ylab("") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        plot.caption = element_text(hjust = 0.01),
        text = element_text(family = "Times New Roman", size = 12)) +
  labs(subtitle = "Tau class glutithione transferases\nSalTrChr07Ag154120", color = "Tissue type") +
  scale_color_brewer(palette = "Dark2") +
  guides(color = "none")
 
p3 <- ggplot(data = pAPX_table) +
  geom_point(mapping = aes(
    x = trt, y = count, color = type)) +
  theme_light() +
  xlab("") + ylab("Normalized gene count") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        plot.caption = element_text(hjust = 0.01),
        text = element_text(family = "Times New Roman", size = 12)) +
  labs(subtitle = "Stroma ascorbate peroxidase\nSalTrChr08Ag177890", color = "Tissue type") +
  scale_color_brewer(palette = "Dark2") +
  guides(color = "none")

p4 <- ggplot(data = SDR1_table) +
  geom_point(mapping = aes(
    x = trt, y = count, color = type)) +
  theme_light() +
  xlab("Treatment") + ylab("") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        plot.caption = element_text(hjust = 0.01),
        text = element_text(family = "Times New Roman", size = 12)) +
  labs(subtitle = "Salt and drought responsive gene\nSalTrChr02Ag049280", color = "Tissue type") +
  scale_color_brewer(palette = "Dark2") +
  guides(color = "none")

p5 <- ggplot(data = NAC_table) +
  geom_point(mapping = aes(
    x = trt, y = count, color = type)) +
  theme_light() +
  xlab("") + ylab("") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        plot.caption = element_text(hjust = 0.01),
        text = element_text(family = "Times New Roman", size = 12)) +
  labs(subtitle = "NAC transcription factor family\nSalTrChr07Ag148570", color = "Tissue type") +
  scale_color_brewer(palette = "Dark2") +
  guides(color = "none")

p6 <- ggplot(data = PRxQ_table) +
  geom_point(mapping = aes(
    x = trt, y = count, color = type)) +
  theme_light() +
  xlab("") + ylab("") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        plot.caption = element_text(hjust = 0.01),
        text = element_text(family = "Times New Roman", size = 12)) +
  labs(subtitle = "Chhloroplast-located peroxiredoxin\nSalTrChr07Bg379770", color = "Tissue type") +
  scale_color_brewer(palette = "Dark2")


p7 <- ggplot(data = MnSOD_table) +
  geom_point(mapping = aes(
    x = trt, y = count, color = type)) +
  theme_light() +
  xlab("") + ylab("") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        plot.caption = element_text(hjust = 0.01),
        text = element_text(family = "Times New Roman", size = 12)) +
  labs(subtitle = "Manganses superoxide dismutase\nSalTrChr07Bg376420", color = "Tissue type") +
  scale_color_brewer(palette = "Dark2") +
  guides(color = "none")

p8 <- ggplot(data = SOS1_table) +
  geom_point(mapping = aes(
    x = trt, y = count, color = type)) +
  theme_light() +
  xlab("") + ylab("") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        plot.caption = element_text(hjust = 0.01),
        text = element_text(family = "Times New Roman", size = 12)) +
  labs(subtitle = "Salt overly sensitive\nSalTrChr04Ag088810", color = "Tissue type") +
  scale_color_brewer(palette = "Dark2") +
  guides(color = "none")

p9 <- ggplot(data = LEA1_table) +
  geom_point(mapping = aes(
    x = trt, y = count, color = type)) +
  theme_light() +
  xlab("") + ylab("") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        plot.caption = element_text(hjust = 0.01),
        text = element_text(family = "Times New Roman", size = 12)) +
  labs(subtitle = "Late embryogenesis abundant\nSalTrChr08Bg401630", color = "Tissue type") +
  scale_color_brewer(palette = "Dark2")

p10 <- ggplot(data = TIP1_table) +
  geom_point(mapping = aes(
    x = trt, y = count, color = type)) +
  theme_light() +
  xlab("") + ylab("") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        plot.caption = element_text(hjust = 0.01),
        text = element_text(family = "Times New Roman", size = 12)) +
  labs(subtitle = "Tonoplast AQP gene\nSalTrChr02Bg244930", color = "Tissue type") +
  scale_color_brewer(palette = "Dark2") +
  guides(color = "none")

########################################################################################

pdf(file = "03_output/05_DESeq2/Figure_files/All_24/Other_mishra_plots_1.pdf", width = 8, height = 8)

(p1 | p2) / (p3 | p10)

dev.off()
dev.off()

pdf(file = "03_output/05_DESeq2/Figure_files/All_24/Other_mishra_plots_2.pdf", width = 8, height = 8)

(p5 | p6) / (p7 | p9)

dev.off()
dev.off()

pdf(file = "03_output/05_DESeq2/Figure_files/All_24/Other_mishra_plots_3.pdf", width = 4, height = 8)

p8 / p4

dev.off()
dev.off()

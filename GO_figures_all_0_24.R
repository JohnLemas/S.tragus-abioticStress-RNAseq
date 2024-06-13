#############
#' Title: Gene Ontology ggplots
#' Author: John Lemas
#' Date: April 28 2024
##############

##############
# packages:
library(ggplot2)
library(tidyr)
library(RColorBrewer)
library(extrafont)

# load fonts

font_import()
y
loadfonts()


# Set directory

setwd("~/SalTrGenome/RNA_Seq/RT_StressTrials/RNA-seq_data/03_output/05_DESeq2/GO_enrichment/")

####### All upregulated DEGs

# Read in the data
Up_all <- read.delim("Up_all_0_24_cleaned.txt", header = T, sep = "\t")

# Clean it up
Up_all <- Up_all[,c(2:4,6,11)]

colnames(Up_all) <- c("Term", "Response", "Count", "p-value", "Fold_Enrichment")

# P-value threshold
Up_all <- Up_all[Up_all$`p-value`<0.05,]

# Plot the data
#ggplot(data = go_table_0_3, 
#       mapping = aes(x = Count, y = reorder(Response, Count), fill = Fold_Enrichment)) +
#  geom_col() + theme_light()

ggplot(data = Up_all, aes(x = -log10(`p-value`), y = reorder(Response, -log10(`p-value`)), 
                                color = Fold_Enrichment, size = Count)) + 
  geom_point() +
  scale_color_gradient(low = "blue", high = "red") +
  theme_light() + 
  ylab("") + 
  xlab("-log(p-value)") + 
  ggtitle("GO enrichment analysis\nA. Both treatments and tissue types 24 HAT") +
  geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "grey") +
  theme(text = element_text(family = "Times New Roman", size = 12))

# Save the plot
pdf("GO_figures/Up_all_dotplot.pdf", width = 8, height = 8)

par(mfrow=c(1,1))

ggplot(data = Up_all, aes(x = -log10(`p-value`), y = reorder(Response, -log10(`p-value`)), 
                          color = Fold_Enrichment, size = Count)) + 
  geom_point() +
  scale_color_gradient(low = "blue", high = "red") +
  theme_light() + 
  ylab("") + 
  xlab("-log(p-value)") + 
  ggtitle("GO enrichment analysis\nA. Both treatments and tissue types 24 HAT") +
  geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "grey") +
  theme(text = element_text(family = "Times New Roman", size = 12))

dev.off()
dev.off()

####### Shoot tissue Upregulated DEGs

# Read in the next table
Up_shoot <- read.delim(file = "Up_shoot_0_24_v2_cleaned.txt", header = T, sep = "\t")

# Clean it up
Up_shoot <- Up_shoot[,c(2:4,6,11)]

colnames(Up_shoot) <- c("Term", "Response", "Count", "p-value", "Fold_Enrichment")

# P-value threshold
Up_shoot <- Up_shoot[Up_shoot$`p-value`<0.05,]

# Plot the data
#ggplot(data = go_table_0_6, 
#       mapping = aes(x = Count, y = reorder(Response, Count), fill = Fold_Enrichment)) +
#  geom_col() + theme_light()

ggplot(data = Up_shoot, aes(x = -log10(`p-value`), y = reorder(Response, -log10(`p-value`)), 
                                color = Fold_Enrichment, size = Count)) + 
  geom_point() +
  scale_color_gradient(low = "blue", high = "red") +
  theme_light() + 
  ylab("") + 
  xlab("-log10(p-value)") + 
  ggtitle("GO enrichment analysis\nB. Cold and salt treated shoot tissues 24 HAT") +
  geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "grey") +
  theme(text = element_text(family = "Times New Roman", size = 12))

# Save the plot
pdf("GO_figures/Up_shoot_24_dotplot.pdf", width = 8, height = 8)

par(mfrow=c(1,1))

ggplot(data = Up_shoot, aes(x = -log10(`p-value`), y = reorder(Response, -log10(`p-value`)), 
                            color = Fold_Enrichment, size = Count)) + 
  geom_point() +
  scale_color_gradient(low = "blue", high = "red") +
  theme_light() + 
  ylab("") + 
  xlab("-log10(p-value)") + 
  ggtitle("GO enrichment analysis\nB. Cold and salt treated shoot tissues 24 HAT") +
  geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "grey") +
  theme(text = element_text(family = "Times New Roman", size = 12))

dev.off()
dev.off()

####### Root tissues upregulated DEGs

# Read in the table
Up_root <- read.delim(file = "Up_root_0_24_cleaned.txt", header = TRUE, sep = "\t")

# Clean it up
Up_root <- Up_root[,c(2:4,6,11)]

colnames(Up_root) <- c("Term", "Response", "Count", "p-value", "Fold_Enrichment")

# P-value threshold
Up_root <- Up_root[Up_root$`p-value`<0.05,]

# Plot
#ggplot(data = go_table_0_8, 
#       mapping = aes(x = Count, y = reorder(Response, Count), fill = Fold_Enrichment)) +
#  geom_col() + theme_light()

ggplot(data = Up_root, aes(x = -log10(`p-value`), y = reorder(Response, -log10(`p-value`)), 
                                color = Fold_Enrichment, size = Count)) + 
  geom_point() +
  scale_color_gradient(low = "blue", high = "red") +
  theme_light() + 
  ylab("") + 
  xlab("-log(p-value)") + 
  ggtitle("GO enrichment analysis\nC. Cold and salt treated root tissues 24 HAT") +
  geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "grey") +
  theme(text = element_text(family = "Times New Roman", size = 12))

# Save the plot:
pdf("GO_figures/Up_root_0_24_dotplot.pdf", width = 8, height = 8)

par(mfrow=c(1,1))

ggplot(data = Up_root, aes(x = -log10(`p-value`), y = reorder(Response, -log10(`p-value`)), 
                           color = Fold_Enrichment, size = Count)) + 
  geom_point() +
  scale_color_gradient(low = "blue", high = "red") +
  theme_light() + 
  ylab("") + 
  xlab("-log(p-value)") + 
  ggtitle("GO enrichment analysis\nC. Cold and salt treated root tissues 24 HAT") +
  geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "grey") +
  theme(text = element_text(family = "Times New Roman", size = 12))

dev.off()
dev.off()

####### Downregulated DEGS
# Shoot tissues

# Read in the table
Down_shoot <- read.delim("Down_shoot_0_24_cleaned.txt", header = T, sep = "\t")

# Clean it up
Down_shoot <- Down_shoot[,c(2:4,6,11)]

colnames(Down_shoot) <- c("Term", "Response", "Count", "p-value", "Fold_Enrichment")

# P-value threshold
Down_shoot <- Down_shoot[Down_shoot$`p-value`<0.05,]

#ggplot(data = go_table_0_24_cold, 
#       mapping = aes(x = Count, y = reorder(Response, Count), fill = Fold_Enrichment)) +
#  geom_col() + theme_light()

# Plot it
ggplot(data = Down_shoot, aes(x = -log10(`p-value`), y = reorder(Response, -log10(`p-value`)), 
                                color = Fold_Enrichment, size = Count)) + 
  geom_point() +
  scale_color_gradient(low = "blue", high = "red") +
  theme_light() + 
  ylab("") + 
  xlab("-log(p-value)") + 
  ggtitle("GO enrichment analysis\nShoot tissues 24 HAT") +
  geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "grey")

# Save the plot:
pdf("GO_figures/Down_shoot_0_24_dotplot.pdf", width = 8, height = 8)

par(mfrow=c(1,1))

ggplot(data = Down_shoot, aes(x = -log10(`p-value`), y = reorder(Response, -log10(`p-value`)), 
                              color = Fold_Enrichment, size = Count)) + 
  geom_point() +
  scale_color_gradient(low = "blue", high = "red") +
  theme_light() + 
  ylab("") + 
  xlab("-log(p-value)") + 
  ggtitle("GO enrichment analysis\nShoot tissues 24 HAT") +
  geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "grey")

dev.off()
dev.off()

####### ROot tissues

# Read in the data
Down_root <- read.delim("Down_root_0_24_cleaned.txt", header = T, sep = "\t")

# Clean it up
Down_root <- Down_root[,c(2:4,6,11)]

colnames(Down_root) <- c("Term", "Response", "Count", "p-value", "Fold_Enrichment")

# P-value threshold
Down_root <- Down_root[Down_root$`p-value`<0.05,]

# Plot
#ggplot(data = go_table_0_24, 
#       mapping = aes(x = Count, y = reorder(Response, Count), fill = Fold_Enrichment)) +
#  geom_col() + theme_light()

ggplot(data = Down_root, aes(x = -log10(`p-value`), y = reorder(Response, -log10(`p-value`)), 
                                      color = Fold_Enrichment, size = Count)) + 
  geom_point() +
  scale_color_gradient(low = "blue", high = "red") +
  theme_light() + 
  ylab("") + 
  xlab("-log(p-value)") + 
  ggtitle("GO enrichment analysis\nRoot tissues 24 HAT") +
  geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "grey")

# Save the plot:
pdf("GO_figures/Down_root_0_24_dotplot.pdf", width = 8, height = 8)

par(mfrow=c(1,1))

ggplot(data = Down_root, aes(x = -log10(`p-value`), y = reorder(Response, -log10(`p-value`)), 
                             color = Fold_Enrichment, size = Count)) + 
  geom_point() +
  scale_color_gradient(low = "blue", high = "red") +
  theme_light() + 
  ylab("") + 
  xlab("-log(p-value)") + 
  ggtitle("GO enrichment analysis\nRoot tissues 24 HAT") +
  geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "grey")

dev.off()
dev.off()



                       
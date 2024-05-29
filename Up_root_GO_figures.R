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

# Set directory

setwd("~/SalTrGenome/RNA_Seq/RT_StressTrials/RNA-seq_data/03_output/05_DESeq2/GO_enrichment/")

########################################################
####### Salt 3 HAT

# Read in the table
go_table_0_3 <- read.table("Up_root_NaCl_0_3_cleaned.txt", header = T, sep = "\t")

# Clean it up
go_table_0_3 <- go_table_0_3[,c(2:4,6,11)]

colnames(go_table_0_3) <- c("Term", "Response", "Count", "p-value", "Fold_Enrichment")

# P-value threshold
go_table_0_3 <- go_table_0_3[go_table_0_3$`p-value`<0.05,]

# Plot it
#ggplot(data = go_table_0_3, 
#       mapping = aes(x = Count, y = reorder(Response, Count), fill = Fold_Enrichment)) +
#  geom_col() + theme_light()

ggplot(data = go_table_0_3, aes(x = -log10(`p-value`), y = reorder(Response, -log10(`p-value`)), 
                                color = Count, size = Fold_Enrichment)) + 
  geom_point() +
  scale_color_gradient(low = "blue", high = "red") +
  theme_light() + 
  ylab("") + 
  xlab("-log(p-value)") + 
  ggtitle("GO enrichment analysis\nA. Salt treated root tissue 3 HAT") +
  geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "grey") +
  theme(text = element_text(family = "Times New Roman", size = 12))

# Save the plot:
pdf("GO_figures/Up_root_NaCl_0_3_dotplot.pdf", width = 8, height = 8)

par(mfrow=c(1,1))

ggplot(data = go_table_0_3, aes(x = -log10(`p-value`), y = reorder(Response, -log10(`p-value`)), 
                                color = Count, size = Fold_Enrichment)) + 
  geom_point() +
  scale_color_gradient(low = "blue", high = "red") +
  theme_light() + 
  ylab("") + 
  xlab("-log(p-value)") + 
  ggtitle("GO enrichment analysis\nA. Salt treated root tissue 3 HAT") +
  geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "grey") +
  theme(text = element_text(family = "Times New Roman", size = 12))

dev.off()
dev.off()

########################################################
####### Cold 6 HAT

# Read in the data
go_table_0_6 <- read.delim(file = "Up_root_Cold_0_6_cleaned.txt", header = T, sep = "\t")

# Clean it up
go_table_0_6 <- go_table_0_6[,c(2:4,6,11)]

colnames(go_table_0_6) <- c("Term", "Response", "Count", "p-value", "Fold_Enrichment")

# P-value threshold
go_table_0_6 <- go_table_0_6[go_table_0_6$`p-value`<0.05,]

# Plot it
#ggplot(data = go_table_0_6, 
#       mapping = aes(x = Count, y = reorder(Response, Count), fill = Fold_Enrichment)) +
#  geom_col() + theme_light()

ggplot(data = go_table_0_6, aes(x = -log10(`p-value`), y = reorder(Response, -log10(`p-value`)), 
                                color = Fold_Enrichment, size = Count)) + 
  geom_point() +
  scale_color_gradient(low = "blue", high = "red") +
  theme_light() + 
  ylab("") + 
  xlab("-log(p-value)") + 
  ggtitle("GO enrichment analysis\nA. Cold treated root tissue 6 HAT") +
  geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "grey") +
  theme(text = element_text(family = "Times New Roman", size = 12))

# Save the plot:
pdf("GO_figures/Up_root_Cold_0_6_dotplot.pdf", width = 8, height = 8)

par(mfrow=c(1,1))

ggplot(data = go_table_0_6, aes(x = -log10(`p-value`), y = reorder(Response, -log10(`p-value`)), 
                                color = Fold_Enrichment, size = Count)) + 
  geom_point() +
  scale_color_gradient(low = "blue", high = "red") +
  theme_light() + 
  ylab("") + 
  xlab("-log(p-value)") + 
  ggtitle("GO enrichment analysis\nA. Cold treated root tissue 6 HAT") +
  geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "grey") +
  theme(text = element_text(family = "Times New Roman", size = 12))

dev.off()
dev.off()

########################################################
####### Salt 8 HAT

# Read in the table
go_table_0_8 <- read.delim(file = "Up_root_NaCl_0_8_cleaned.txt", header = TRUE, sep = "\t")

# Clean it up
go_table_0_8 <- go_table_0_8[,c(2:4,6,11)]

colnames(go_table_0_8) <- c("Term", "Response", "Count", "p-value", "Fold_Enrichment")

# P-value threshold
go_table_0_8 <- go_table_0_8[go_table_0_8$`p-value`<0.05,]

# Plot it
#ggplot(data = go_table_0_8, 
#       mapping = aes(x = Count, y = reorder(Response, Count), fill = Fold_Enrichment)) +
#  geom_col() + theme_light()

ggplot(data = go_table_0_8, aes(x = -log10(`p-value`), y = reorder(Response, -log10(`p-value`)), 
                                color = Fold_Enrichment, size = Count)) + 
  geom_point() +
  scale_color_gradient(low = "blue", high = "red") +
  theme_light() + 
  ylab("") + 
  xlab("-log(p-value)") + 
  ggtitle("GO enrichment analysis\nB. Salt treated root tissue 8 HAT") +
  geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "grey") +
  theme(text = element_text(family = "Times New Roman", size = 12))

# Save the plot:
pdf("GO_figures/Up_root_NaCl_0_8_dotplot.pdf", width = 8, height = 8)

par(mfrow=c(1,1))

ggplot(data = go_table_0_8, aes(x = -log10(`p-value`), y = reorder(Response, -log10(`p-value`)), 
                                color = Fold_Enrichment, size = Count)) + 
  geom_point() +
  scale_color_gradient(low = "blue", high = "red") +
  theme_light() + 
  ylab("") + 
  xlab("-log(p-value)") + 
  ggtitle("GO enrichment analysis\nB. Salt treated root tissue 8 HAT") +
  geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "grey") +
  theme(text = element_text(family = "Times New Roman", size = 12))

dev.off()
dev.off()

########################################################
####### Cold 24 HAT

#Read in the table
go_table_0_24_cold <- read.delim("Up_root_Cold_0_24_cleaned.txt", header = T, sep = "\t")

# Clean it up
go_table_0_24_cold <- go_table_0_24_cold[,c(2:4,6,11)]

colnames(go_table_0_24_cold) <- c("Term", "Response", "Count", "p-value", "Fold_Enrichment")

# P-value threshold
go_table_0_24_cold <- go_table_0_24_cold[go_table_0_24_cold$`p-value`<0.05,]

# Plot it
#ggplot(data = go_table_0_24_cold, 
#       mapping = aes(x = Count, y = reorder(Response, Count), fill = Fold_Enrichment)) +
#  geom_col() + theme_light()

ggplot(data = go_table_0_24_cold, aes(x = -log10(`p-value`), y = reorder(Response, -log10(`p-value`)), 
                                      color = Fold_Enrichment, size = Count)) + 
  geom_point() +
  scale_color_gradient(low = "blue", high = "red") +
  theme_light() + 
  ylab("") + 
  xlab("-log(p-value)") + 
  ggtitle("GO enrichment analysis\nB. Cold treated root tissue 24 HAT") +
  geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "grey") +
  theme(text = element_text(family = "Times New Roman", size = 12))

# Save the plot:
pdf("GO_figures/Up_root_Cold_0_24_dotplot.pdf", width = 8, height = 8)

par(mfrow=c(1,1))

ggplot(data = go_table_0_24_cold, aes(x = -log10(`p-value`), y = reorder(Response, -log10(`p-value`)), 
                                      color = Fold_Enrichment, size = Count)) + 
  geom_point() +
  scale_color_gradient(low = "blue", high = "red") +
  theme_light() + 
  ylab("") + 
  xlab("-log(p-value)") + 
  ggtitle("GO enrichment analysis\nB. Cold treated root tissue 24 HAT") +
  geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "grey") +
  theme(text = element_text(family = "Times New Roman", size = 12))

dev.off()
dev.off()

########################################################
####### Salt shoot 24 HAT

# Read in the table
go_table_0_24 <- read.delim("Up_root_NaCl_0_24_cleaned.txt", header = T, sep = "\t")

# Clean it up
go_table_0_24 <- go_table_0_24[,c(2:4,6,11)]

colnames(go_table_0_24) <- c("Term", "Response", "Count", "p-value", "Fold_Enrichment")

# P-value threshold
go_table_0_24 <- go_table_0_24[go_table_0_24$`p-value`<0.05,]

# Plot it
#ggplot(data = go_table_0_24, 
#       mapping = aes(x = Count, y = reorder(Response, Count), fill = Fold_Enrichment)) +
#  geom_col() + theme_light()

ggplot(data = go_table_0_24, aes(x = -log10(`p-value`), y = reorder(Response, -log10(`p-value`)), 
                                 color = Fold_Enrichment, size = Count)) + 
  geom_point() +
  scale_color_gradient(low = "blue", high = "red") +
  theme_light() + 
  ylab("") + 
  xlab("-log(p-value)") + 
  ggtitle("GO enrichment analysis\nC. Salt treated root tissue 24 HAT") +
  geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "grey") +
  theme(text = element_text(family = "Times New Roman", size = 12))

# Save the plot:
pdf("GO_figures/Up_root_NaCl_0_24_dotplot.pdf", width = 8, height = 8)

par(mfrow=c(1,1))

ggplot(data = go_table_0_24, aes(x = -log10(`p-value`), y = reorder(Response, -log10(`p-value`)), 
                                 color = Fold_Enrichment, size = Count)) + 
  geom_point() +
  scale_color_gradient(low = "blue", high = "red") +
  theme_light() + 
  ylab("") + 
  xlab("-log(p-value)") + 
  ggtitle("GO enrichment analysis\nC. Salt treated root tissue 24 HAT") +
  geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "grey") +
  theme(text = element_text(family = "Times New Roman", size = 12))

dev.off()
dev.off()



##*****************************
##Thomas Tekle - 1079669
##
##University of Guelph
##March 11, 2025
##
##Bioinformatics for Genomics - Assignment 2 
##
##******************************


library(edgeR)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(glue)
library(ggvenn)
library(gt)


# Load count matrix (raw counts from featureCounts or similar tool)
counts <- read.delim("//wsl.localhost/Ubuntu/home/hb/Genomics_A2/count_output.txt", comment.char = "#")

# Remove unnecessary columns (Start, End, Length, Strand) to retain only gene counts
Seq_counts <- counts %>% 
  select(-c(Start,End,Length,Strand)) 

# Remove duplicated gene entries and rows with missing values 
Seq_counts.with.chr <- Seq_counts[!duplicated(counts$Geneid),] %>% 
  na.omit()

# Rename columns to reflect sample groups and biological replicates
col_names <- c("Gene_ID","Early_1","Early_2","Early_3","Thin_1","Thin_2","Thin_3","Mature_1", "Mature_2", "Mature_3")

# Remove chromosome column and assign new column names
Seq_counts.no.chr <- Seq_counts.with.chr %>% 
  select(-c(Chr))


colnames(Seq_counts.no.chr) <- col_names

# Create DGEList object for edgeR analysis, specifying counts and gene IDs
Seq_counts.DGE <- DGEList(counts=Seq_counts.no.chr[,2:10], genes=Seq_counts.no.chr[,1])

# Define treatment groups
groups <- as.factor(c("Early","Early","Early","Thin","Thin","Thin","Mature","Mature","Mature"))


Seq_counts.DGE$samples$group <- groups

# Filter out lowly expressed genes to retain only those expressed above threshold
Filtered.Seq_counts <- Seq_counts.DGE[filterByExpr(Seq_counts.DGE, group = Seq_counts.DGE$samples$group),,keep.lib.size = F]

# Normalize counts using TMM normalization
normal.GeneExpr <- calcNormFactors(Filtered.Seq_counts)

# Create design matrix for the linear model (no intercept to compare all groups directly)
design <- model.matrix(~0 + normal.GeneExpr$samples$group) %>% 
  `colnames<-`(levels(groups))

# Fit a quasi-likelihood negative binomial GLM
fit <- glmQLFit(normal.GeneExpr, design = design)

# Define contrasts for pairwise comparisons:
# 1. Early vs Mature
contrast.earlyVsMature <- makeContrasts(EarlyVsMature = Early - Mature, levels = design)

# 2. Early vs average of Thin and Mature
contrast.earlyVsMature_Thin <- makeContrasts(EarlyVsThin_Mature = Early - (Thin + Mature)/2, levels = design)

# Perform likelihood ratio test (LRT) for both contrasts
lrtFit.E_M <- glmLRT(fit, contrast = contrast.earlyVsMature)
lrtFit.E_TM <- glmLRT(fit, contrast = contrast.earlyVsMature_Thin)

# Function to classify DEGs based on logFC and FDR thresholds
diffexpress.checker <- function(data){
  data$table$differential_expression <- "Insignificant"
  data$table$differential_expression[data$table$logFC > 1 & data$table$FDR < 0.05]  <- "Up Regulated"
  data$table$differential_expression[data$table$logFC < -1 & data$table$FDR < 0.05]  <- "Down Regulated"
  
  return(data$table)
}

# Extract and classify DEGs for both contrasts
topTags_EM <- topTags(lrtFit.E_M, n = Inf) %>% 
  diffexpress.checker()

topTags_ETM <- topTags(lrtFit.E_TM, n = Inf) %>% 
  diffexpress.checker()

#Volcano Plot for both DEG lists

#Function to create volcano plots for DEGs
Volcano_Plotter <- function(data, treatment, show_legend = TRUE){
  plot <- ggplot(data = data, aes(x = logFC, y = -log10(FDR), colour = differential_expression))+
  geom_point()+
  geom_vline(xintercept = c(-1, 1), col = "gray", linetype = 'dashed')+
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
  geom_point(size = 2) +
  scale_color_manual(values = c("blue", "grey", "red"),labels = c("Downregulated", "Not significant", "Upregulated"))+
  labs(title = glue("Early vs {treatment}"),
       x = "Log2 Fold Change",
       y = "-Log10 FDR") 

  if (!show_legend) {
    plot <- plot + theme(legend.position = "none")
  }
  return(plot)
}

# Generate volcano plots for both contrasts
EM_volcano <- Volcano_Plotter(topTags_EM, "Mature", show_legend = FALSE)
ETM_volcano <- Volcano_Plotter(topTags_ETM, "Thin and Mature")

# Arrange volcano plots side by side
ggarrange(EM_volcano, ETM_volcano)

#Top DEGs between different treatment comparisons
# Function to create bar plot of top N DEGs

topCharter <- function(data, top_n, comparison){
  
top.genes <- head(topTags_EM[order(data$FDR),], top_n)

graph <- ggplot(top.genes, aes(x = reorder(genes, logFC), y = logFC, fill = logFC > 0)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_manual(values = c("blue", "red")) +
  labs(title = glue("{comparison}"),
       x = "Gene",
       y = "Log2 Fold Change") +
  theme_minimal()+
  theme(legend.position = "none")

return(graph)
}

# Create bar plots for top 20 DEGs from both contrasts
top.genes.EM <- topCharter(topTags_EM, 20, "Early vs Mature")
top.genes.ETM <- topCharter (topTags_ETM, 20, "Early vs Thin and Mature")


ggarrange(top.genes.EM, top.genes.ETM) 


#Venn Diagram to discern the number of shared DEGs between the two groups

# Extract significantly up- or down-regulated genes for Venn diagram
top.EM.DEG <- topTags_EM %>% 
  filter(differential_expression == c("Up Regulated", "Down Regulated"))

top.ETM.DEG <- topTags_ETM %>% 
  filter(differential_expression == c("Up Regulated", "Down Regulated"))

#Create list for Venn diagram
venn_list <- list(
  `Early vs Mature` = top.EM.DEG$genes,
  `Early vs Thin and Mature` = top.ETM.DEG$genes
)

#Plot Venn diagram to visualize shared and unique DEGs
ggvenn(venn_list, 
       fill_color = c("blue", "yellow"), 
       stroke_size = 0.5, 
       set_name_size = 5,)

#Making table to show the top genes for each list of DEGs 
head(top.EM.DEG[order(top.EM.DEG$FDR),],10) %>% 
  gt()

head(top.ETM.DEG[order(top.ETM.DEG$FDR),],10) %>% 
  gt()


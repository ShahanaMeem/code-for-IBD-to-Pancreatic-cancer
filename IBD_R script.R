#install packages
install.packages("plotly")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("EnhancedVolcano")
install.packages("ggplot2")

#load library
library("readr")
library("readxl")
library("ggplot2")
library("DESeq2")
library("ggrepel")
library("plotly")
library("dplyr")
library("tidyverse")
library("EnhancedVolcano")

#load_data
sampleTable <- read_excel("D:/Pancreatic cancer Vs IBD/Github Files/IBD_GSM2.xlsx")
rownames(sampleTable) <- sampleTable$GSM
count_data <- read_excel("D:/Pancreatic cancer Vs IBD/Github Files/IBD_FPKM_files.xlsx")
rownames(count_data) <- count_data[[1]]
count_data <- count_data[,-1] 
colnames(count_data)
dim(count_data)
dim(sampleTable)
duplicated_genes <- duplicated(count_data$Gene_ID)
if (any(duplicated_genes)) {
  print("Duplicate Gene IDs found:")
  print(count_data[duplicated_genes, ])
}
print(duplicated_genes)
colnames(count_data) <- sampleTable$GSM
count_data <- round(count_data)
dds <- DESeqDataSetFromMatrix(countData = count_data,colData = sampleTable, design = ~ Conditions)
dds_ibd <- DESeq(dds)
res_1 <- results(dds_ibd, contrast = c("Conditions", "Wild", "mCD"))
res_2 <- results(dds_ibd, contrast = c("Conditions", "Wild", "mUC"))
res_1
res_2

plotMA(res_1, ylim=c(-2,2))
plotMA(res_2, ylim=c(-2,2))

#up_down_regulation_res_1
filter_res_1 <- (res_1[ which(res_1$padj <= .01 & abs(res_1$log2FoldChange) > 1.5), ])
filter_res_1 <- as.data.frame(filter_res_1)
filter_res_1 <- filter_res_1 %>% mutate( de_status = case_when(
  log2FoldChange > 1.5 ~ "up",
  log2FoldChange < -1.5 ~ "down"
))
count_data1 <- read_excel("D:/Pancreatic cancer Vs IBD/Github Files/IBD_FPKM_file.xlsx")
res_1_df <- as.data.frame(res_1)
gene_ids <- count_data1$Gene_ID
res_1_df <- cbind(Gene_ID = gene_ids, res_1_df)
res_1_df
# View the final filtered and formatted dataframe
View(res_1_df)
# Load the required library
library(openxlsx)

# Specify the file path and name
output_file <- "D:/Pancreatic cancer Vs IBD/Github Files/Filtered_DEGs_1.xlsx"

# Write the dataframe to an Excel file
write.xlsx(res_1_df, file = output_file)

# Confirm the file has been written
cat("File successfully saved at:", output_file)

#up_down_regulation_res_2
filter_res_2 <- (res_2[ which(res_2$padj <= .01 & abs(res_2$log2FoldChange) > 1.5), ])
filter_res_2 <- as.data.frame(filter_res_2)
filter_res_2 <- filter_res_2 %>% mutate( de_status = case_when(
  log2FoldChange > 1.5 ~ "up",
  log2FoldChange < -1.5 ~ "down"
))
count_data2 <- read_excel("D:/Pancreatic cancer Vs IBD/Github Files/IBD_FPKM_file.xlsx")
res_2_df <- as.data.frame(res_2)
gene_ids <- count_data2$Gene_ID
res_2_df <- cbind(Gene_ID = gene_ids, res_2_df)
res_2_df
# View the final filtered and formatted dataframe
View(res_2_df)
# Load the required library
library(openxlsx)

# Specify the file path and name
output_file <- "D:/Pancreatic cancer Vs IBD/Github Files/Filtered_DEGs_2.xlsx"

# Write the dataframe to an Excel file
write.xlsx(res_2_df, file = output_file)

# Confirm the file has been written
cat("File successfully saved at:", output_file)

#Res_1_volcano_plot
res_df_1 <- as.data.frame(res_1)
res_df_1$significant <- (res_df_1$padj < 0.01)  # Good practice to define significance threshold
res_df_1$color <- "Not Significant"  # Default color for non-significant
res_df_1$color[res_df_1$padj < 0.01 & res_df_1$log2FoldChange > 1.5] <- "Upregulated"   # Change color for upregulated
res_df_1$color[res_df_1$padj < 0.01 & res_df_1$log2FoldChange < -1.5] <- "Downregulated"  # Change color for downregulated
ggplot(res_df_1, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = color), alpha = 0.5) +
  geom_hline(yintercept = -log10(0.01), linetype = "dashed", color = "blue") +  # Line for significance level
  geom_vline(xintercept = c(-1.5, 1.5), linetype = "dashed", color = "blue") +  # Lines for fold change thresholds
  scale_color_manual(values = c("Not Significant" = "grey", "Upregulated" = "red", "Downregulated" = "darkgreen")) +
  labs(title = "Healthy vs CD",  # Adding the new title here
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted P-value") +
  theme_minimal() +
  theme(legend.position = "bottom")

#Res_2_volcano_plot
res_df_2 <- as.data.frame(res_2)
res_df_2$significant <- (res_df_2$padj < 0.01)  # Good practice to define significance threshold
res_df_2$color <- "Not Significant"  # Default color for non-significant
res_df_2$color[res_df_2$padj < 0.01 & res_df_2$log2FoldChange > 1.5] <- "Upregulated"   # Change color for upregulated
res_df_2$color[res_df_2$padj < 0.01 & res_df_2$log2FoldChange < -1.5] <- "Downregulated"  # Change color for downregulated
ggplot(res_df_2, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = color), alpha = 0.5) +
  geom_hline(yintercept = -log10(0.01), linetype = "dashed", color = "blue") +  # Line for significance level
  geom_vline(xintercept = c(-1.5, 1.5), linetype = "dashed", color = "blue") +  # Lines for fold change thresholds
  scale_color_manual(values = c("Not Significant" = "grey", "Upregulated" = "red", "Downregulated" = "darkgreen")) +
  labs(title = "Healthy vs UC",  # Adding the new title here
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted P-value") +
  theme_minimal() +
  theme(legend.position = "bottom")


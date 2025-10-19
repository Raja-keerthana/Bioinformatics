library(DESeq2)


counts <- read.csv("/path_to/Counts.csv", row.names = 1,header=T, check.names=F )


counts <- round(counts)


metadata <- read.csv("/path_to/metadata.csv", row.names = 1,header=T, check.names=F)


all(rownames(metadata) %in% colnames(counts))
all(colnames(counts) == rownames(metadata))



dds <- DESeqDataSetFromMatrix(countData = counts, colData = metadata, design = ~ Condition)
dds

# Run DESeq
dds <- DESeq(dds)
res <- results(dds, contrast = c("Condition", "WT", "MUT"))
res
summary(res)



res$Gene <- rownames(res)
res <- res[, c("Gene", setdiff(colnames(res), "Gene"))]


write.table(res, "/path_to/results/deseq2/deseq_results.txt", sep = "\t", row.names = FALSE, col.names = TRUE)


upregulated_genes <- res[!is.na(res$padj) & res$log2FoldChange > 0 & res$padj < 0.05, "Gene"]


downregulated_genes <- res[!is.na(res$padj) & res$log2FoldChange < 0 & res$padj < 0.05, "Gene"]


upregulated_genes
downregulated_genes



write.table(upregulated_genes, "/path_to/results/deseq2/upregulated_genes.txt", sep = "\t", row.names = FALSE, col.names = FALSE)
write.table(downregulated_genes, "/path_to/results/deseq2/downregulated_genes.txt", sep = "\t", row.names = FALSE, col.names = FALSE)


#####################################################################################################
                                            Volcano plot
#####################################################################################################

library(ggplot2)
library(ggrepel)
library(dplyr)   


results_df <- as.data.frame(res)
head(results_df)


results_df$significance <- "Not Significant"
results_df$significance[results_df$padj < 0.05 & results_df$log2FoldChange > 0] <- "Upregulated"
results_df$significance[results_df$padj < 0.05 & results_df$log2FoldChange < 0] <- "Downregulated"

volcano_plot <- ggplot(results_df, aes(x = log2FoldChange, y = -log10(padj), color = significance)) +
  geom_point(alpha = 0.6, size = 2) +
  theme_minimal() +
  labs(title = "Volcano Plot",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted P-value") +
  scale_color_manual(values = c("red", "gray", "blue")) +
  geom_text_repel(data = subset(results_df, significance != "Not Significant"),
                  aes(label = Gene), 
                  vjust = 1.5, 
                  size = 3.5)

print(volcano_plot)

ggsave(filename = "/path_to/results/deseq2/volcano_plot.png", plot = volcano_plot, width = 8, height = 6, dpi = 300)

################################################################################################
                                            Heatmap
################################################################################################
library(DESeq2)
library(pheatmap)


upregulated_genes <- res[!is.na(res$padj) & res$log2FoldChange > 0 & res$padj < 0.05, "Gene"]
downregulated_genes <- res[!is.na(res$padj) & res$log2FoldChange < 0 & res$padj < 0.05, "Gene"]
significant_genes <- unique(c(upregulated_genes, downregulated_genes))


normalized_counts <- assay(vst(dds, blind = FALSE))[significant_genes, ]

meta.data <- as.data.frame(colData(dds))


meta.data <- meta.data[order(meta.data$Condition), ]


heat_map <- pheatmap(normalized_counts,
         scale = "row",                  
         show_rownames = TRUE,            
         show_colnames = TRUE,            
         main = "Heatmap of Significant Genes",
         cluster_rows = FALSE,             
         cluster_cols = FALSE,            
         annotation_col = meta.data[, "Condition", drop = FALSE], 
         annotation_colors = list(Condition = c("WT" = "blue", "MUT" = "red")))  
print(heat_map)

graphics.off()
png("/path_to/results/deseq2/heatmap.png", width = 800, height = 600, res = 300)
print(heat_map)
dev.off()



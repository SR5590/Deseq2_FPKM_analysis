# install and load the necessary packages
#install.packages("DESeq2")
library(DESeq2)
library(tidyverse)
library(edgeR)
library(ggplot2)
library(factoextra)
library(gplots)
library(RColorBrewer)
library(scatterplot3d)
library(plotly)  #new library

# set the path to the count data file
count_data_file = "FeatureCount_GRCm39.tsv"

# read the count data into R
count_data = read.table(count_data_file, header = TRUE, row.names = 1)

# set the path to the metadata file
metadata_file = "sample_info.csv"

# read the metadata into R
metadata = read.csv(metadata_file, header = TRUE, row.names = 1)
# read the metadata into R
#metadata = read.csv("sample_info.csv", header = TRUE, sep = ",")
#metadata$Type = as.factor(metadata$Type)
# create the design formula for the analysis
# (this will depend on the specific data and analysis you are doing)

# create the DESeqDataSet object using the count data and metadata
dds = DESeqDataSetFromMatrix(countData = count_data, colData = metadata, design = ~Type)

# run the differential expression analysis
dds = DESeq(dds)
dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)

# fit the negative binomial GLM with increased maximum iterations
dds <- nbinomWaldTest(dds)

# extract the results of the analysis
res = results(dds)

# view the results
#head(results)

# Creation of normalized counts matrix
normCounts <- counts(dds, normalized = TRUE)
normCounts_tibble <-
  normCounts %>% as_tibble(rownames = "Gene_ID")

# Creation of log2 normalized counts matrix
lognormcounts <- log2(normCounts + 1)
lognormcounts_tibble <-
  lognormcounts %>% as_tibble(rownames = "Gene_ID")

# Sort DESEQ results by p-value as tibble
res_tibble <-
  res[order(res$padj),] %>% as_tibble(rownames = "Gene_ID")

# Combine DESEQ results and log2 normalized counts as tibble
resdata_tibble <-
  full_join(res_tibble, lognormcounts_tibble, by = "Gene_ID")

# Regularized log transformation
rld <- rlog(dds, blind = FALSE)

# Creation of a list of TOP50 genes (lowest padj values) from DESEQ results
res_tibble_top_genes50 <-
  res_tibble %>% arrange(padj) %>% head(50) %>% arrange(desc(log2FoldChange))

# Creation of a list of TOP250 genes (lowest padj values) from DESEQ results
res_tibble_top_genes250 <-
  res_tibble %>% arrange(padj) %>% head(250) %>% arrange(desc(log2FoldChange))

# Creation of a list of genes 2fold up and 2fold down from DESEQ results
# !! all p-values will be included
res_tibble_top_genes_2fold_up_and_down <-
  res_tibble %>% filter(log2FoldChange >= 2 |
                          log2FoldChange <= -2) %>% arrange(desc(log2FoldChange))

################################################################################
#### DATA EXPORT ####
dir.create("02_deseq", showWarnings = FALSE)

write.table(
  normCounts_tibble,
  file = "02_deseq/01_normalized_counts_german.tsv",
  sep = "\t",
  dec = ",",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)
write.table(
  normCounts_tibble,
  file = "02_deseq/01_normalized_counts_english.tsv",
  sep = "\t",
  dec = ".",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)
write.table(
  res_tibble,
  file = "02_deseq/02_DEseq_result_german.tsv",
  sep = "\t",
  dec = ",",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)
write.table(
  res_tibble,
  file = "02_deseq/02_DEseq_result_english.tsv",
  sep = "\t",
  dec = ".",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)
write.table(
  resdata_tibble,
  file = "02_deseq/03_log_normalized_counts_german.tsv",
  sep = "\t",
  dec = ",",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)
write.table(
  resdata_tibble,
  file = "02_deseq/03_log_normalized_counts_english.tsv",
  sep = "\t",
  dec = ".",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)
write.table(
  res_tibble_top_genes250,
  file = "02_deseq/04_TOP250_genes_german_pvalue.tsv",
  append = FALSE,
  sep = "\t",
  dec = ",",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)
write.table(
  res_tibble_top_genes250,
  file = "02_deseq/04_TOP250_genes_english_pvalue.tsv",
  append = FALSE,
  sep = "\t",
  dec = ".",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)
write.table(
  res_tibble_top_genes_2fold_up_and_down,
  file = "02_deseq/05_TOP_genes_german_2fold.tsv",
  append = FALSE,
  sep = "\t",
  dec = ",",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)
write.table(
  res_tibble_top_genes_2fold_up_and_down,
  file = "02_deseq/05_TOP_genes_english_2fold.tsv",
  append = FALSE,
  sep = "\t",
  dec = ".",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

################################################################################
#### GRAPHICAL DATA EXPORT ####

pdf("02_deseq/Graphical_Output.pdf", paper = "A4")

# Dispersions plot
# Gene dispersion is the key determinant of the read count bias in RNA-seq data analysis
# Dispersions plot shows the per-gene dispersion estimates together with the fitted mean-dispersion relationship.
plotDispEsts(
  dds,
  main = "Dispersion plot",
  legend = TRUE,
  genecol = "black",
  fitcol = "red",
  finalcol = "dodgerblue",
  cex = 0.25,
  log = "xy"
)

# Histogram of log2-fold changes
vn_genes_2fold_up <- res_tibble %>% filter(log2FoldChange >= 2) %>% nrow()
vn_genes_2fold_up_p_value <- res_tibble %>% filter(log2FoldChange >= 2, padj <= 0.05) %>% nrow()
vn_genes_2fold_down <- res_tibble %>% filter(log2FoldChange <= -2) %>% nrow()
vn_genes_2fold_down_p_value <- res_tibble %>% filter(log2FoldChange <= -2, padj <= 0.05) %>% nrow()


# Save differentially expressed genes to a TSV file
differentially_expressed_genes <- res_tibble %>% filter(abs(log2FoldChange) >= 2, padj <= 0.05)
write.table(differentially_expressed_genes, file = "differentially_expressed_genes.tsv", sep = "\t", quote = FALSE, row.names = FALSE)


ggplot(data = res_tibble, aes(x = log2FoldChange)) +
  geom_histogram(aes(y = ..density..), colour = "black", fill = "white", bins = 30) +
  geom_vline(aes(xintercept=0),
             color = "blue", linetype = "dashed", size = 0.5) +
  geom_vline(aes(xintercept=2),
             color = "red", linetype = "dashed", size = 0.5) +
  geom_vline(aes(xintercept=-2),
             color = "red", linetype ="dashed", size = 0.5) +
  geom_density(alpha=0.2, fill="#FF6666") +
  annotate("text", x = 4.5, y = 1.5, label = paste0(vn_genes_2fold_up, "\ngenes\n\n\n",vn_genes_2fold_up_p_value, "\nsign. genes\n(padj<0.05)")) +
  annotate("text", x = -4.5, y = 1.5, label = paste0(vn_genes_2fold_down, "\ngenes\n\n\n",vn_genes_2fold_down_p_value, "\nsign. genes\n(padj<0.05)")) +
  ggtitle("Distribution of differentially expressed genes") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black")) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 2.5)) +
  scale_x_continuous(limits = c(-6, 6)) +
  xlab("log2 fold changes") +
  ylab("density")


# Density Histogram of p-Values
ggplot(data = res_tibble, aes(x = pvalue)) +
  geom_histogram(aes(y = ..density..), colour = "#00000000", fill = "#00000000") +
  geom_vline(aes(xintercept=0.05),
             color = "red", linetype ="dashed", size = 0.5) +
  geom_density(alpha=0.2, fill="#FF6666") +
  ggtitle("Distribution of p-values") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black")) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, NA)) +
  xlab("p-values") +
  ylab("density")


# PCAAnalysis-1
select <- order(rowMeans(lognormcounts), decreasing = TRUE)[1:40000]
hc <- lognormcounts[select,]
pca <- t(hc)
pca <- prcomp(pca)
group <- as.factor(metadata$Type)
fviz_pca_ind(
  pca,
  col.ind = group,
  palette = colorRampPalette(c("#00AFBB", "#FC4E07"))(7),  #colors <-  c("#00AFBB", "#FC4E07")
  addEllipses = TRUE,
  legend.title = "Group",
  repel = TRUE
)



# PCAAnalysis-2
select <- order(rowMeans(lognormcounts), decreasing = TRUE)[1:40000]
hc <- lognormcounts[select,]
pca <- t(hc)
pca <- prcomp(pca)
pca_df <- data.frame(PCA1 = pca$x[,1], PCA2 = pca$x[,2], PCA3 = pca$x[,3])
pca_df$group <- metadata$Type
plot_ly(pca_df, x = pca_df$PCA1, y = pca_df$PCA2, z = pca_df$PCA3, type = "scatter3d", mode = "markers") %>%
  add_markers(color = pca_df$group, colorscale = 'Plasma', text = pca_df$group) %>%
  layout(scene = list(xaxis = list(title = "PC1", showgrid = T, zeroline = F, showline = T, mirror = T),
                      yaxis = list(title = "PC2", showgrid = T, zeroline = F, showline = T, mirror = T),
                      zaxis = list(title = "PC3", showgrid = T, zeroline = F, showline = T, mirror = T)),
         paper_bgcolor = 'rgba(0,0,0,0)', plot_bgcolor = 'rgba(0,0,0,0)')
# PCAAnalysis-3
plotPCA(rld, intgroup = "Type") +
  scale_color_brewer(palette = "Set1")
se <- SummarizedExperiment(log2(counts(dds, normalized = TRUE) + 1), colData =
                             colData(dds))

plotPCA(DESeqTransform(se), intgroup = "Type") +
  scale_color_brewer(palette = "Set1")





# PCAAnalysis 5
mycols <- brewer.pal(8, "Set2")[1:length(unique(group))]
rld_pca <-
  function (rld,
            intgroup = "Type",
            ntop = 500,
            colors = NULL,
            main = "PCA Biplot",
            legendpos = "bottomleft",
            textcx = 0.75,
            ...) {
    require(genefilter)
    require(calibrate)
    require(RColorBrewer)
    rv = rowVars(assay(rld))
    select = order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
    pca = prcomp(t(assay(rld)[select,]))
    fac = factor(apply(as.data.frame(colData(rld)[, intgroup, drop = FALSE]), 1, paste, collapse = " : "))
    if (is.null(colors)) {
      if (nlevels(fac) >= 3) {
        colors = brewer.pal(nlevels(fac), "Paired")
      }   else {
        colors = c("black", "red")
      }
    }
    pc1var <- round(summary(pca)$importance[2, 1] * 100, digits = 1)
    pc2var <- round(summary(pca)$importance[2, 2] * 100, digits = 1)
    pc1lab <- paste0("PC1 (", as.character(pc1var), "%)")
    pc2lab <- paste0("PC1 (", as.character(pc2var), "%)")
    plot(
      PC2 ~ PC1,
      data = as.data.frame(pca$x),
      bg = colors[fac],
      pch = 21,
      xlab = pc1lab,
      ylab = pc2lab,
      main = main,
      ...
    )
    with(as.data.frame(pca$x),
         textxy(PC1, PC2, labs = rownames(as.data.frame(pca$x)), cex = textcx))
    legend(legendpos,
           legend = levels(fac),
           col = colors,
           pch = 20)
  }
rld_pca(rld,
        colors = mycols,
        intgroup = "Batch",
        xlim = c(-100, 100))

# Sample distance heatmap 1
sampleDists <- assay(rld) %>% t() %>% dist()
sampleDistMatrix <- as.matrix(sampleDists)
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
Heatmap(
  name = "Heatmap of sample distances",
  matrix = sampleDistMatrix,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  col = colors,
  show_heatmap_legend = TRUE,
  heatmap_legend_param = list(title = "Matrix")
)

# Sample distance heatmap 2
heatmap.2(
  as.matrix(sampleDistMatrix),
  main = "Heatmap of sample distances",
  key = FALSE,
  trace = "none",
  col = colors,
  ColSideColors = mycols[group],
  RowSideColors = mycols[group],
  margin = c(10, 10),
)


# Volcano Plot 2
#EnhancedVolcano(res,
#lab = rownames(res),
#x = 'log2FoldChange',
#y = 'pvalue')


# Heatmap Genes padj < 0.05
select_genes <- rownames(subset(res, padj < 0.05))
lognormcounts_heatmap <- lognormcounts
heatmap.2(
  as.matrix(lognormcounts_heatmap)[select_genes, ],
  #col = colorRampPalette(c("red", "green")),
  main = "Heatmap genes with p-value < 0.05",
  scale = "row",
  Rowv = TRUE,
  Colv = TRUE,
  cexRow = 0.5,
  cexCol = 1,
  srtCol = 45,
  labRow = "",
  trace = "none",
  dendrogram = "both",
  density.info = "histogram",
  margins = c(8, 8),
  lmat=rbind(c(0, 3), c(2, 1), c(0, 4)),
  lhei=c(1.5, 5, 2)
)
colnames(lognormcounts_heatmap) <- group
heatmap.2(
  as.matrix(lognormcounts_heatmap)[select_genes, ],
  #col = colorRampPalette(c("red", "green")),
  main = "Heatmap genes with p-value < 0.05",
  scale = "row",
  Rowv = TRUE,
  Colv = TRUE,
  cexRow = 0.5,
  cexCol = 1,
  srtCol = 45,
  trace = "none",
  dendrogram = "both",
  density.info = "histogram",
  margins = c(8, 8),
  lmat=rbind(c(0, 3), c(2, 1), c(0, 4)),
  lhei=c(1.5, 5, 2)
)


#####
# Create a color palette from blue to red
palette <- colorRampPalette(c("blue", "white", "red"))

# Use the palette as the color for the heatmap
heatmap.2(
  as.matrix(lognormcounts_heatmap)[select_genes, ],
  col = palette,
  main = "Heatmap genes with p-value < 0.05",
  scale = "row",
  Rowv = TRUE,
  Colv = TRUE,
  cexRow = 0.5,
  cexCol = 1,
  srtCol = 45,
  labRow = "",
  trace = "none",
  dendrogram = "both",
  density.info = "histogram",
  margins = c(8, 8),
  lmat=rbind(c(0, 3), c(2, 1), c(0, 4)),
  lhei=c(1.5, 5, 2)
)

colnames(lognormcounts_heatmap) <- group

# Use the palette as the color for the heatmap
heatmap.2(
  as.matrix(lognormcounts_heatmap)[select_genes, ],
  col = palette,
  main = "Heatmap genes with p-value < 0.05",
  scale = "row",
  Rowv = TRUE,
  Colv = TRUE,
  cexRow = 0.5,
  cexCol = 1,
  srtCol = 45,
  trace = "none",
  dendrogram = "both",
  density.info = "histogram",
  margins = c(8, 8),
  lmat=rbind(c(0, 3), c(2, 1), c(0, 4)),
  lhei=c(1.5, 5, 2)
)

# Remove gene names from the heatmap
heatmap.2(
  as.matrix(lognormcounts_heatmap)[select_genes, ],
  col = palette,
  main = "Heatmap genes with p-value < 0.05",
  scale = "row",
  Rowv = TRUE,
  Colv = TRUE,
  cexRow = 0.5,
  cexCol = 1,
  srtCol = 45,
  labRow = "",  # Set labRow to an empty string
  trace = "none",
  dendrogram = "both",
  density.info = "histogram",
  margins = c(8, 8),
  lmat=rbind(c(0, 3), c(2, 1), c(0, 4)),
  lhei=c(1.5, 5, 2)
)


# Heatmap Top50-genes with lowest p-value
mat  <- assay(rld)[res_tibble_top_genes50$Gene_ID, ]
mat  <- mat - rowMeans(mat)
heatmap.2(
  mat,
  main = "Heatmap Top50-genes with lowest p-value",
  trace = "none",
  scale = "row",
  dendrogram = "both",
  density.info = "histogram",
  cexRow = 0.5,
  cexCol = 1,
  srtCol = 45,
  margins = c(5, 8),
  lmat=rbind(c(0, 3), c(2, 1), c(0, 4)),
  lhei=c(1.5, 5, 2)
)

#Heatmap(mat)
colnames(mat) <- group
heatmap.2(
  mat,
  main = "Heatmap Top50-genes with lowest p-value",
  trace = "none",
  scale = "row",
  dendrogram = "both",
  density.info = "histogram",
  cexRow = 0.5,
  cexCol = 1,
  srtCol = 45,
  margins = c(8, 8),
  lmat=rbind(c(0, 3), c(2, 1), c(0, 4)),
  lhei=c(1.5, 5, 2)
)

dev.off()

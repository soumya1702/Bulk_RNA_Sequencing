countdata <- read.table("~/Desktop/counts.txt", header=TRUE, row.names=1)

# Remove first five columns (chr, start, end, strand, length)
countdata <- countdata[ ,6:ncol(countdata)]

head(countdata)
# Remove .bam or .sam from filenames
colnames(countdata) <- gsub("\\Aligned.sortedByCoord.out.bam", "", colnames(countdata))

# Convert to matrix
countdata <- as.matrix(countdata)
head(countdata)
(celltype <- factor(c(rep("KID90",7),rep("KID21",7),"KID19","KID18",rep("KID06",7),rep("KID04",6),rep("KID19",6),rep("KID18",6),"KID04")))

#####################
#metadata#
library(DESeq2)
library(readxl)

metadata <- read_excel("~/Desktop/cancer.xlsx")


# Add sampleID's to the mapping file
metadata$sampleid <- row.names(metadata)

# Reorder sampleID's to match featureCounts column order. 
metadata <- metadata[match(colnames(countdata), metadata$Run), ]

# Make sure ID's are correct
head(metadata)


############################################################
# - countData : count dataframe
# - colData : sample metadata in the dataframe with row names as sampleID's
# - design : The design of the comparisons to use. 
#            Use (~) before the name of the column variable to compare
ddsMat <- DESeqDataSetFromMatrix(countData = countdata,
                                 colData = metadata,
                                 design = ~celltype)


# Find differential expressed genes
ddsMat <- DESeq(ddsMat)

###########################################################
# Get results from testing with FDR adjust pvalues
results <- results(ddsMat, pAdjustMethod = "fdr", alpha = 0.05)

# Generate summary of testing. 
summary(results)

####################################################################
mcols(results, use.names = T)

###########################################################
# Mouse genome database (Select the correct one)
library(org.Hs.eg.db) 

# Add gene full name
results$description <- mapIds(x = org.Hs.eg.db,
                              keys = row.names(results),
                              column = "GENENAME",
                              keytype = "SYMBOL",
                              multiVals = "first")

# Add gene symbol
results$symbol <- row.names(results)

# Add ENTREZ ID
results$entrez <- mapIds(x = org.Hs.eg.db,
                         keys = row.names(results),
                         column = "ENTREZID",
                         keytype = "SYMBOL",
                         multiVals = "first")

# Add ENSEMBL
results$ensembl <- mapIds(x = org.Hs.eg.db,
                          keys = row.names(results),
                          column = "ENSEMBL",
                          keytype = "SYMBOL",
                          multiVals = "first")

# Subset for only significant genes (q < 0.05)
results_sig <- subset(results, padj < 0.05)
head(results_sig)


##################################################################################################
# Convert all samples to rlog
ddsMat_rlog <- rlog(ddsMat, blind = FALSE)
library(ggplot2)
# Plot PCA by column variable
plotPCA(ddsMat_rlog, intgroup = "celltype", ntop = 500) +
  theme_bw() + # remove default ggplot2 theme
  geom_point(size = 5) + # Increase point size
  scale_y_continuous(limits = c(-20, 20)) + # change limits to fix figure dimensions
  ggtitle(label = "Principal Component Analysis (PCA)", 
          subtitle = "Top 500 most variable genes") 

#########################################################################################
# Convert all samples to rlog
ddsMat_rlog <- rlog(ddsMat, blind = FALSE)

# Gather 30 significant genes and make matrix
mat <- assay(ddsMat_rlog[row.names(results_sig)])[1:40, ]

# Choose which column variables you want to annotate the columns by.
annotation_col = data.frame(
  Group = factor(colData(ddsMat_rlog)$celltype), 
  row.names = colData(ddsMat_rlog)$Run
)

# Specify colors you want to annotate the columns by.
ann_colors = list(
  Group = c(CD4_T = "lightblue", CD8_T = "forestgreen", gd_T = "orange")
)


# Make Heatmap with pheatmap function.
## See more in documentation for customization
pheatmap(mat = mat, 
         color = colorRampPalette(brewer.pal(9, "YlOrBr"))(255), 
         scale = "row",
         annotation_col = annotation_col, # Add multiple annotations to the samples
         annotation_colors = ann_colors,# Scale genes to Z-score (how many standard deviations)
         # Change the default colors of the annotations
         fontsize = 6.5, # Make fonts smaller
         cellwidth = 5, # Make the cells wider
         show_colnames = T)

#########################################
# Gather Log-fold change and FDR-corrected pvalues from DESeq2 results
## - Change pvalues to -log10 (1.3 = 0.05)
data <- data.frame(gene = row.names(results),
                   pval = -log10(results$padj), 
                   lfc = results$log2FoldChange)

# Remove any rows that have NA as an entry
data <- na.omit(data)
library(dplyr)
# Color the points which are up or down
## If fold-change > 0 and pvalue > 1.3 (Increased significant)
## If fold-change < 0 and pvalue > 1.3 (Decreased significant)
data <- mutate(data, color = case_when(data$lfc > 0 & data$pval > 1.3 ~ "Increased",
                                       data$lfc < 0 & data$pval > 1.3 ~ "Decreased",
                                       data$pval < 1.3 ~ "nonsignificant"))

# Make a basic ggplot2 object with x-y values
vol <- ggplot(data, aes(x = lfc, y = pval, color = color))

# Add ggplot2 layers
vol +   
  ggtitle(label = "Volcano Plot", subtitle = "Colored by fold-change direction") +
  geom_point(size = 2.5, alpha = 0.8, na.rm = T) +
  scale_color_manual(name = "Directionality",
                     values = c(Increased = "#008B00", Decreased = "#CD4F39", nonsignificant = "darkgray")) +
  theme_bw(base_size = 5) + # change overall theme
  theme(legend.position = "right") + # change the legend
  xlab(expression(log[2]("LoGlu" / "HiGlu"))) + # Change X-Axis label
  ylab(expression(-log[10]("adjusted p-value"))) + # Change Y-Axis label
  geom_hline(yintercept = 1.3, colour = "darkgrey") + # Add p-adj value cutoff line
  scale_y_continuous(trans = "log1p") # Scale yaxis due to large p-values

#######################################################################################
plotMA(results, ylim = c(-5, 5))

###############################################################
plotDispEsts(ddsMat)

########################################
# Convert all samples to rlog
ddsMat_rlog <- rlog(ddsMat, blind = FALSE)

# Get gene with highest expression
top_gene <- rownames(results)[which.min(results$log2FoldChange)]

# Plot single gene
plotCounts(dds = ddsMat, 
           gene = top_gene, 
           intgroup = "celltype", 
           normalized = T, 
           transform = T)

##################################################################

# Remove any genes that do not have any entrez identifiers
#results_sig_entrez <- subset(results_sig, is.na(entrez) == FALSE)

# Create a matrix of gene log2 fold changes
gene_matrix <- results_sig$log2FoldChange

# Add the entrezID's as names for each logFC entry
names(gene_matrix) <- results_sig$entrez

# View the format of the gene matrix
##- Names = ENTREZ ID
##- Values = Log2 Fold changes
head(gene_matrix)
################################################################
kegg_enrich <- enrichKEGG(gene = names(gene_matrix),
                          organism = 'human',
                          pvalueCutoff = 0.05, 
                          qvalueCutoff = 0.10)

# Plot results
barplot(kegg_enrich, 
        drop = TRUE, 
        showCategory = 10, 
        title = "KEGG Enrichment Pathways",
        font.size = 8)

#####################################################
go_enrich <- enrichGO(gene = names(gene_matrix),
                      OrgDb = 'org.Hs.eg.db', 
                      readable = T,
                      ont = "BP",
                      pvalueCutoff = 0.05, 
                      qvalueCutoff = 0.10)

# Plot results
barplot(go_enrich, 
        drop = TRUE, 
        showCategory = 10, 
        title = "GO Biological Pathways",
        font.size = 8)
### PACKAGE INSTALLATION & LOADING
# This block ensures all packages are compatible with R 4.4.0

if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")

# Force update BiocManager for R 4.4.0
if (BiocManager::version() != "3.19") {
  BiocManager::install(version = "3.19", ask = FALSE, update = TRUE)
}

# List of required packages
packages <- c("affy", "limma", "ggplot2", "hgu133plus2.db", "AnnotationDbi", 
              "clusterProfiler", "enrichplot", "pheatmap", "plotly", "factoextra", "ggrepel")

# Install missing packages
install_if_missing <- function(p) {
  if (!require(p, character.only = TRUE)) {
    BiocManager::install(p, ask = FALSE)
    library(p, character.only = TRUE)
  }
}

invisible(lapply(packages, install_if_missing))

# Standard Loading
library(affy)
library(limma)
library(ggplot2)
library(hgu133plus2.db)
library(AnnotationDbi)
library(clusterProfiler)
library(enrichplot)
library(pheatmap)
library(plotly)
library(factoextra)
library(ggrepel)

###LOADING TARGETS FILE 
#loading our target file 
targets <- readTargets("./data/data.txt")

###BACKGROUND CORRECTION AND NORMALISATION.
#Loading the data and Normalisation 

eset <- justRMA(filenames = targets$Sample_ID, celfile.path = "./data/input/")

#Setting the coloum names to Title

colnames(eset) <- targets$Title

#ANNOTATION
#Annotation 

Alias <- rownames(exprs(eset))
fData(eset) <- AnnotationDbi::select(hgu133plus2.db,Alias,c("SYMBOL","ENTREZID","GENENAME"),keytype="PROBEID")

#Identify duplicated probe IDs in the feature data
duplicated_probes <- fData(eset)$PROBEID[duplicated(fData(eset)$PROBEID)]


# filtering out duplicated probes 
fData(eset) <- fData(eset)[!duplicated(fData(eset)$PROBEID),]

# Removing NA using symbol 
HasSymbol <- !is.na(fData(eset)$SYMBOL)
eset <- eset[HasSymbol, ]

# DIFFERENTIAL EXPRESSION 

#Generating a model matrix specifying your design 
Group <- factor(targets$Condition, levels = c("Healthy", "Arthritis"))
design <- model.matrix(~0 + Group)
colnames(design) = levels(factor(Group))

#Fit the model to the design 
fit <- lmFit(eset, design)

#Making contrast
contrast <- makeContrasts('Arthritis-Healthy', levels = design)

#Empirical Bayes model
fit2 <- contrasts.fit(fit, contrast)
fit2 <- eBayes(fit2)

#Extracting our differential expressed results
Res <- topTable(fit2, number=Inf, adjust.method = 'fdr')

# Summary of our results
Results <- decideTests(fit2)

#creating a new coloum to classify our genes as upregulated or downregulated.

Res$deff.express.gene <- "NO"
# if log2FC > 0 and adjusted p-value < 0.08, set as "UP" 
Res$deff.express.gene[Res$logFC > 0 & Res$adj.P.Val < 0.08] <- "UP"
# if log2FC < 0 and adjusted p-value < 0.08, set as "DOWN"
Res$deff.express.gene[Res$logFC < 0 & Res$adj.P.Val < 0.08] <- "DOWN"

Res$delabel <- NA

Res$delabel[Res$deff.express.gene != "NO"] <- Res$SYMBOL [Res$deff.express.gene != "NO"]

## GENE ENRICHMENT ANALYSIS
#GO Gene Enrichment Analysis for the differential expressed gene

#Define our universe gene symbols (total gene symbols differential expressed)
Universe <- Res$SYMBOL

#Define genes symbols that are upregulated or downregulated
Up_genes <- Res$SYMBOL[Res$deff.express.gene == "UP"]
Down_genes <- Res$SYMBOL[Res$deff.express.gene == "DOWN"]

#Enrichment analysis of upregulated genes

enrich_up <- enrichGO(
  Up_genes,
  hgu133plus2.db,
  keyType = "SYMBOL",
  ont = "ALL",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.2,
  Universe,
  minGSSize = 10,
  maxGSSize = 1000,
  readable = FALSE,
  pool = FALSE
)

#Enrichment analysis of Downregulated genes
enrich_Down <- enrichGO(
  Down_genes,
  hgu133plus2.db,
  keyType = "SYMBOL",
  ont = "ALL",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.2,
  Universe,
  minGSSize = 10,
  maxGSSize = 1000,
  readable = FALSE,
  pool = FALSE
)

# RESULTS - VISUALIZATION

#CHECKING OUR EXPRESSION VALUES AND NORMALIZATION RESULTS
#showing our target file 
targets
#viewing our expression data 
head(exprs(eset))
#Boxplot of of our expression data after normalisation 
boxplot(exprs(eset), main = "Expression data after normalisation", col = c("blue", "red"))

### PCA PLOT
# Prepare data for PCA (transpose so samples are rows)
pca_data <- t(exprs(eset))
pca_res <- prcomp(pca_data, scale. = TRUE)

# Visualize PCA
fviz_pca_ind(pca_res,
             col.ind = Group, # Color by condition
             palette = c("#00AFBB", "#FC4E07"),
             addEllipses = TRUE, # Concentration ellipses
             legend.title = "Condition",
             repel = TRUE,
             title = "PCA - Sample Clustering")

###RESULTS FOR ANNOTATION
#list our keytypes
keytypes(hgu133plus2.db)
#viewing our annotated data 
head(fData(eset))
#viewing our Duplicates probes
print(duplicated_probes)
#showing summary of our logic report
table(HasSymbol)
#checking the dimension for the eset 
dim(eset)

###RESULTS FOR DIFFERENTIAL EXPRESSION
#Viewing the differential expression results
head(Res)
#showing the summary of the results. 
summary(Results)
# Volcano plot Differential expressed gene (Static) with top genes labeled
ggplot(Res, aes(y = -log10(adj.P.Val), x = logFC, label = delabel)) +
  geom_point(aes(colour = deff.express.gene), alpha = 0.6, size = 1.2) +
  geom_text(data = subset(Res, logFC < -1 | logFC > 1), aes(label = delabel), nudge_y = 0.3, size = 2) +
  scale_colour_manual(values = c("UP" = "red", "NO" = "lightgrey", "DOWN" = "blue")) +
  geom_vline(xintercept = c(-1, 1), colour = "black", linetype = "dashed") +
  geom_hline(yintercept = -log10(0.08), colour = "black", linetype = "dashed") +
  labs(y = "-log10(adj. p-value)", x = "log2(Fold Change)", colour = "Condition") +
  ggtitle("Volcano Plot for DEGs (Top labels)") +
  theme_minimal()

# Interactive Volcano plot
p <- ggplot(Res, aes(y = -log10(adj.P.Val), x = logFC, text = paste("Gene:", SYMBOL, "<br>Adj.P:", round(adj.P.Val, 4)))) +
  geom_point(aes(colour= deff.express.gene), alpha = 0.5) + 
  scale_colour_manual(values = c("UP" = "red", "NO" = "lightgrey", "DOWN" = "blue")) +
  geom_vline(xintercept = c(-1, 1), colour = "black", linetype = "dashed") +
  geom_hline(yintercept = -log10(0.08), colour = "black", linetype = "dashed")+
  labs(y = "-log10(adj. p-value)", x = "log2(Fold Change)") +
  theme_minimal()

ggplotly(p, tooltip = "text")

### HEATMAP OF TOP UPREGULATED AND DOWNREGULATED GENES
# Filter for UP and DOWN genes
up_genes <- subset(Res, deff.express.gene == "UP")
down_genes <- subset(Res, deff.express.gene == "DOWN")

# Select top 25 upregulated and top 25 downregulated genes sorted by absolute logFC
up_heat <- head(up_genes[order(-up_genes$logFC), ], 25)
down_heat <- head(down_genes[order(down_genes$logFC), ], 25)

# Combine into a single dataframe
top_heat <- rbind(up_heat, down_heat)

gene_matrix <- exprs(eset)[rownames(top_heat), ]
# Replace probe IDs with Gene Symbols
rownames(gene_matrix) <- top_heat$SYMBOL

# Create heatmap
pheatmap(gene_matrix, 
         scale = "row", 
         clustering_distance_rows = "correlation",
         annotation_col = data.frame(Condition = targets$Condition, row.names = targets$Title),
         main = "Top Upregulated and Downregulated Genes",
         color = colorRampPalette(c("blue", "white", "red"))(100))


###RESULTS FOR GENE ENRICHMENT ANALYSIS.
#Enrichment Annalysis

#For upregulated genes
p2 <- dotplot(enrich_up, showCategory=10) + ggtitle("Upregulated Gene ontolgy Enrichment analysis")
p2
#For downregulated genes 
p3 <- dotplot(enrich_Down, showCategory=10) + ggtitle("Downregulated Gene ontolgy Enrichment analysis")
p3

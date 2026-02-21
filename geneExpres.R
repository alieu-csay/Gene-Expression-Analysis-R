###DOWNLOADING AND LOADING REQUIRED PACKAGES 

#Downloading and Loading Required Packages

if (!require("BiocManager", quietly = TRUE))
   install.packages("BiocManager")

#From the Bioconductor
#BiocManager::install("affy")
#BiocManager::install("limma")
#BiocManager::install("hgu133plus2.db")
#BiocManager::install("AnnotationDbi")
#BiocManager::install("clusterProfiler", version = "3.18")
#BiocManager::install("pathview")
#BiocManager::install("enrichplot")

#CRAN packages
#install.packages("ggplot2")
#install.packages("pheatmap")

#Loading the librarys 
library(affy)
library(limma)
library(ggplot2)
library(hgu133plus2.db)
library(AnnotationDbi)
library(clusterProfiler)
library(enrichplot)
library(pheatmap)

###LOADING TARGETS FILE 
#loading our target file 
targets <- readTargets("/Users/ALIEU CSAY/Desktop/Masters Bioinformatics/Principles of Bioinformatics/E-portfolio/E_portfolio_ass/data.txt")

###BACKGROUND CORRECTION AND NORMALISATION.
#Loading the data and Normalisation 

eset <- justRMA(filenames = targets$Sample_ID, celfile.path = "/Users/ALIEU CSAY/Desktop/Masters Bioinformatics/Principles of Bioinformatics/E-portfolio/E_portfolio_ass/input/")

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
# if log2FC > 0 and adj.pvalue < 0.05, set as "UP" 
Res$deff.express.gene[Res$logFC > 0 & Res$adj.P.Val < 0.08] <- "UP"
# if log2FC < 0 and adj.pvalue < 0.05, set as "DOWN"
Res$deff.express.gene[Res$logFC < 0 & Res$adj.P.Val < 0.08] <- "DOWN"

Res$delabel <- NA

Res$delabel[Res$deff.express.gene != "NO"] <- Res$SYMBOL [Res$deff.express.gene != "NO"]

## GENE ENRICHMENT ANALYSIS
#Differential expression 

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
# if log2FC > 0 and adj.pvalue < 0.05, set as "UP" 
Res$deff.express.gene[Res$logFC > 0 & Res$adj.P.Val < 0.08] <- "UP"
# if log2FC < 0 and adj.pvalue < 0.05, set as "DOWN"
Res$deff.express.gene[Res$logFC < 0 & Res$adj.P.Val < 0.08] <- "DOWN"

Res$delabel <- NA

Res$delabel[Res$deff.express.gene != "NO"] <- Res$SYMBOL [Res$deff.express.gene != "NO"]

```

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
#Volcano plot Differential expressed gene
ggplot(Res, aes(y = -log10(P.Value), x = logFC))+
geom_text(data = subset(Res, logFC < -1 | logFC >  1),aes(label = delabel),nudge_y = 0.3, size = 2) +
geom_point() + labs(y = "-log10(p-value)", x = "log2(Fold Change)") +
ggtitle("Volcano Plot for DEGs") +
geom_point(aes(colour= Res$deff.express.gene)) +
scale_colour_manual(values = c("red","yellow","blue")) +
geom_vline(xintercept = c(-1, 1), colour = "red" ) +
geom_hline(yintercept = -log10(0.08), colour = "green")+
theme_minimal() 


###RESULTS FOR GENE ENRICHMENT ANALYSIS.
#Enrichment Annalysis

#For upregulated genes
p2 <- dotplot(enrich_up, showCategory=10) + ggtitle("Upregulated Gene ontolgy Enrichment analysis")
p2
#For downregulated genes 
p3 <- dotplot(enrich_Down, showCategory=10) + ggtitle("Downregulated Gene ontolgy Enrichment analysis")
p3



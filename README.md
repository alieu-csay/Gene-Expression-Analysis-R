R Code for Differential Gene Expression Analysis (Microarray Data)
---------------------------------------------------------------------------------------------------------------

__Project Title__

Differential Expression of Genes in Osteoarthritis (OA) and Normal Synovial Tissue

Overview
-----------------------------------------------------------------------------------------------------------------
This repository contains R code for performing microarray-based differential gene expression analysis comparing Osteoarthritis (OA) synovial tissue with Normal (Healthy) synovial tissue.

The pipeline includes:

1. Loading and preprocessing Affymetrix CEL files

2. Background correction and normalization (RMA)

3. Probe annotation

4. Differential expression analysis using limma

5. Volcano plot visualization

6. Gene Ontology (GO) enrichment analysis

Note: The workflow follows standard Bioconductor microarray analysis practices.

Data Source
--------------------------------------------------------------------------------------------------------------------
The dataset was obtained from:

1. Gene Expression Omnibus (GEO)

2. Platform: Affymetrix HG-U133_Plus_2 Array

3. Raw data format: .CEL files


Required R Packages
-----------------------------------------------------------------------------------------------------------------------

__Bioconductor Packages__

1. affy
  
2. limma
  
3. hgu133plus2.db
  
4. AnnotationDbi
  
5. clusterProfiler
  
6. enrichplot

__CRAN Packages__

1. ggplot2

2. pheatmap




R Code for Differential Gene Expression Analysis (Microarray Data)
---------------------------------------------------------------------------------------------------------------

__Project Title__

Differential Expression of Genes in Osteoarthritis (OA) and Normal Synovial Tissue

Check the Full project:
[https://alieu-csay.github.io/Gene-Expression-Analysis-R/](https://alieu-csay.github.io/Gene-Expression-Analysis-R/)

Overview
-----------------------------------------------------------------------------------------------------------------
This repository contains R code for performing microarray-based differential gene expression analysis comparing Osteoarthritis (OA) synovial tissue with Normal (Healthy) synovial tissue.

The pipeline includes:

1. Loading and preprocessing Affymetrix CEL files

2. QC and normalization (RMA)
  
3. EDA 

5. Probe annotation

6. Differential expression analysis using limma

7. Volcano plot and Heatmap visualization

8. Gene Ontology (GO) enrichment analysis

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




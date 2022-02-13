# GenomicSuperSignature
**Interpretation of RNA-seq experiments through robust, efficient comparison to public databases**

### PURPOSE
Thousands of RNA sequencing profiles have been deposited in public archives, yet 
remain unused for the interpretation of most newly performed experiments. Methods 
for leveraging these public resources have focused on the interpretation of existing 
data, or analysis of new datasets independently, but do not facilitate direct comparison 
of new to existing experiments. The interpretability of common unsupervised analysis 
methods such as Principal Component Analysis would be enhanced by efficient comparison 
of the results to previously published datasets.

### METHODS
To help identify replicable and interpretable axes of variation in any given gene 
expression dataset, we performed principal component analysis (PCA) on 536 studies 
comprising 44,890 RNA sequencing profiles. Sufficiently similar loading vectors, 
when compared across studies, were combined through simple averaging. We annotated 
the collection of resulting average loading vectors, which we call Replicable Axes 
of Variation (RAV), with details from the originating studies and gene set enrichment 
analysis. Functions to match PCA of new datasets to RAVs from existing studies, 
extract interpretable annotations, and provide intuitive visualization, are implemented 
as the GenomicSuperSignature R package, to be submitted to Bioconductor. 

### Results
Preprint is available [**HERE**](https://www.biorxiv.org/content/10.1101/2021.05.26.445900v1). 
Usecases and benchmarking examples included in the preprint are further documented at [GenomicSuperSignaturePaper](https://shbrief.github.io/GenomicSuperSignaturePaper/) page.




## Installation
You can install GenomicSuperSignature in Bioconductor. This can be done using 
`BiocManager`:
```
if (!require("BiocManager"))
    install.packages("BiocManager")

library(BiocManager)
install("GenomicSuperSignature")
```

RAVmodel can be directly downloaded from Google bucket with no cost. The sizes of 
RAVmodels`RAVmodel_C2.rds` and `RAVmodel_PLIERpriors.rds` are 476.1MB and 475.1MB, 
respectively. You can use `wget` or `GenomicSuperSignature::getModel` function.

```
## Download RAVmodel with wget
wget https://storage.googleapis.com/genomic_super_signature/RAVmodel_C2.rds
wget https://storage.googleapis.com/genomic_super_signature/RAVmodel_PLIERpriors.rds

## Download RAVmodel with getModel function
getModel("C2")
getModel("PLIERpriors")
```

## Schematic
#### Overview of GenomicSuperSignature 
Schematic illustration of RAVmodel construction and GenomicSuperSignature application.  Building the RAVmodel (components in grey) is performed once on a time scale of hours on a high-memory, high-storage server. Users can apply RAVmodel on their data (component in red) using the GenomicSuperSignature R/Bioconductor package (components in blue), which operates on a time scale of seconds for exploratory data analyses (components in orange) on a typical laptop computer.   

<img src="https://raw.githubusercontent.com/shbrief/GenomicSuperSignature/master/vignettes/GSig_overview.png" width="90%" height="90%"/>

<br>

#### User's perspective
The GenomicSuperSignature package allows users to access a RAVmodel (Z matrix, blue) and annotation information on each RAV. From a gene expression matrix (Y matrix, grey), users can calculate dataset-level validation score or sample score matrix (B matrix, red). Through the RAV of your interest, additional information such as related studies, GSEA, and MeSH terms can be easily extracted. 

<img src="https://raw.githubusercontent.com/shbrief/GenomicSuperSignature/master/vignettes/GSig_model_usage_diagram.png" width="90%" height="90%"/>

#### Information assembled by GenomicSuperSignature
GenomicSuperSignature connects different public databases and prior information through RAVindex, creating the knowledge graph illustrated here. Users can instantly access data and metadata resources from multiple entry points, such as gene expression profiles, MeSH terms, gene sets, and keywords. 

<img src="https://raw.githubusercontent.com/shbrief/GenomicSuperSignaturePaper/master/inst/images/GSig_knowledge_graph.png" width="90%" height="90%"/>



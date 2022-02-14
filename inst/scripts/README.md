# GenomciSuperSignature Command Line Interface (CLI) Tool

## How to use
The minimum inputs you need to provide is,    
1) count matrix, 2) your choice of RAVmodel, and 3) the output directory. 

```
## Example
Rscript gss.R --input bcellViperExpr_10C.tsv --model microRAVmodel_v2.rds --outDir .
```

## What it does
Connect new gene expression profile with the relevant information from the 
existing databases, such as previous publications, MeSH terms, and gene sets.


## Inputs
### Count Files
GenomicSuperSignature takes a count matrix as an input. Input file should have
row names in gene symbols and column names in sample ID. For the best result,
we recommend a data transformation (e.g. log2) for the input to follow a 
normal distribution, while scaling is NOT recommended. Currently for 
validation, inputs need at least eight samples.

Example of input format:

|      |GSM44075  |GSM44078 |GSM44080  |GSM44081|
|------|----------|---------|----------|--------|
|ADA   |9.571369 |10.599436 |8.740659 |10.104469|
|CDH2  |6.175890  |5.312704 |5.651928  |4.462205|
|MED6  |9.671113  |8.773383 |9.190276  |9.526235|
|NR2E3 |9.847733  |9.582061 |9.628792  |8.820422|

### RAVmodel
*R*eplicable *A*xes of *V*ariation (RAV) consists of principal components 
repeatedly observed in an independent analysis of multiple published datasets.
RAVs connect different databases that are both linked to the originated study 
or associated with the RAV itself through the gene rankings of it. RAVmodel 
contains the collection of RAVs (RAVindex), metadata from model building 
process and the additional annotations. Currently, two RAVmodels are available
based on the gene sets used for annotation.

1) C2 : RAVmodel annotated with Molecular Signatures Database (MSigDB) 
curated gene sets (version 7.1)    
2) PLIERpriors : RAVmodel annotated with the three gene sets provided in the 
[PLIER package](https://github.com/wgmao/PLIER) - bloodCellMarkersIRISDMAP, 
svmMarkers, and canonicalPathways


## Outputs
There are four categories of outputs from this tool, which is one html file and
three csv tabular files. The actual number of csv files will vary depending on
the parameter, *--numOut*, and the validated RAVs.

#### validate.csv

| Column |Description  |
|------|----------|
|score | the maximum pearson correlation coefficient between the top 8 PCs of the input and RAVs|
|PC | one of the top 8 PCs of the input, which gives the highest *score* |
|sw | the average silhouette width of the RAV |
|cl_size | the size of each RAV |
|cl_num | the RAV number |

#### Genesets 
This is the enriched gene sets for the target RAV, calculated from the ranked 
gene list. Gene sets with the adjusted p-value < 0.05 are included.

| Column |Description  |
|------|----------|
|Description | name of the gene sets |
|NES | normalized enrichment score (ES) |
|pvalue | statistical significance |
|qvalues | p-value adjusted for the FDR |

#### Literatures  

| Column  |Description            |
|---------|-----------------------|
|studyName|study accession        |
|title    |the title of the study |

#### report.html      
A html file with the summary of the main analyses by GenomicSuperSignature. 
It includes MeSH terms in word cloud and an interactive plot overviewing the 
validated RAVs, in addition to the previews of the tabular output files.  


## Citations
Oh, S., Geistlinger, L., Ramos, M., Taroni, J.N., Carey, V.J., Greene, C.S., 
Waldron, L., & Davis, S.R. (2021). GenomicSuperSignature: interpretation of 
RNA-seq experiments through robust, efficient comparison to public databases. 
bioRxiv.

## References
GenomicSuperSignature package: [webpage](https://shbrief.github.io/GenomicSuperSignature/)       
GenomicSuperSignature usecases: [webpage](https://shbrief.github.io/GenomicSuperSignaturePaper/)

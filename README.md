[![Travis-CI Build Status](https://travis-ci.org/waldronlab/ProjectAsPackage.svg?branch=master)](https://travis-ci.org/waldronlab/ProjectAsPackage)
[![Coverage Status](https://codecov.io/github/waldronlab/ProjectAsPackage/coverage.svg?branch=master)](https://codecov.io/github/waldronlab/ProjectAsPackage?branch=master)


# PCAGenomicSignatures

## Software connecting new dataset to the public expression database

## Installation
```
devtools::install_github("shbrief/PCAGenomicSignatures")
```

## Schematic
Here is the overview of using PCAGenomicSignatures.
<img src="vignettes/GSig_model_usage_diagram.png"/>

PCAGenomicSignatures allows you to connect your gene expression data to the existing 
database through the expression profile itself. You only need to provide your gene
expression matrix.
<img src="vignettes/GSig_knowledge_network.png"/>


## `refinebioData.csv`
Information on refine.bio RNAseq-compendium data downloaded on 04.10.20.    
- studyName: Name of the study   
- downloaded: the number of samples downloaded for each study   
- metadata: the number of samples for each study based on the metadata   
- imported: whether `tximport` succeedded or not   
- PCAmodel_677: whether the study is a part of 677 studies for PCAmodel (test ver.)   
- PCAmodel_1399: whether the study is a part of 1,399 studies for PCAmodel (production ver.)   
- PCAmodel_536: whether the study is a part of 536 studies for PCAmodel (production ver.2),
include studies with >50 samples and succesfully imported

[Temp.] This table combines `numSample.csv` and `rseq_meta.csv`


## `MeSH_terms_1399refinebio.rds`
MeSH terms assigned to 1,399 studies from refine.bio. 

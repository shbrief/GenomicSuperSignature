### `studyMeta.csv`
Information on refine.bio RNAseq-compendium data downloaded on 04.10.20. You can
find out how this table is created from [`available_samples.Rmd`]().

- **studyName**: Name of the study   
- **downloaded**: The number of samples downloaded for each study   
- **metadata**: The number of samples for each study based on the metadata   
- **imported**: Whether `tximport::tximport` succeeded or not   
- **PCAmodel_677**: Whether the study is a part of 677 studies for PCAmodel (test ver.1)   
- **PCAmodel_1399**: Whether the study is a part of 1,399 studies for PCAmodel (test ver.2)   
- **PCAmodel_536**: Whether the study is a part of 536 studies for PCAmodel (production ver.),
include studies with >50 samples and succesfully imported

### `MeSH_terms_1399refinebio.rds`
MeSH terms assigned to 1,399 studies, which is `TRUE` in `PCAmodel_1399` column
of `studyMeta.csv` file. 

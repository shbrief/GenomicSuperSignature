### `studyMeta.csv`
Information on refine.bio RNAseq-compendium data downloaded on 04.10.20. You can
find out how this table is created from [here](https://github.com/shbrief/GenomicSuperSignaturePaper/blob/master/Methods/training_datasets/available_samples.Rmd).

- **studyName**: Name of the study   
- **downloaded**: The number of samples downloaded for each study   
- **metadata**: The number of samples for each study based on the metadata   
- **imported**: Whether `tximport::tximport` succeeded or not   
- **RAVmodel_677**: Whether the study is a part of 677 studies for RAVmodel (test ver.1)   
- **RAVmodel_1399**: Whether the study is a part of 1,399 studies for RAVmodel (test ver.2)   
- **RAVmodel_536**: Whether the study is a part of 536 studies for RAVmodel (production ver.),
include studies with >50 samples and successfully imported



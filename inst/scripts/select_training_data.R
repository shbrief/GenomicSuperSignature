##### How to select training datasets ##########################################

# Exclude studies with > 1,000 samples
x <- read.table("/inst/extdata/refinebioData.csv")
ind_1000 <- which(x$metadata > 1000)
studies_1000 <- x$studyName[ind_1000]   # 38 studies

# Exclude studies assigned with "Single-Cell Analysis" MeSH term
y <- readRDS("/inst/extdata/MeSH_terms_1399refinebio.rds")
ind_sca <- which(y$name == "Single-Cell Analysis")
studies_sca <- y$identifier[ind_sca]   # 30 studies

all_sc_studies <- unique(c(studies_1000, studies_sca))   # 67 studies

# Studies with > 50 downloaded samples and succesfully imported
sum(x$downloaded > 50 & x$imported == TRUE & !x$studyName %in% all_sc_studies)   # 536 studies
ind <- which(x$downloaded > 50 & x$imported == TRUE & !x$studyName %in% all_sc_studies)

x$PCAmodel_536 <- FALSE
x$PCAmodel_536[ind] <- TRUE

write.table(x, "~/data2/PCAGenomicSignatures/inst/extdata/refinebioData.csv")

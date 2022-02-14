## CLI test dataset
## Expression matrix of `dset` from bcellViper package

library(bcellViper)
data(bcellViper)
write.table(exprs(dset),
            file="inst/scripts/bcellViperExpr.tsv", row.names=TRUE, sep="\t")

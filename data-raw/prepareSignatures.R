library(devtools)
library(roxygen2)

Signatures <- RJSONIO::fromJSON("signatures")
usethis::use_data(Signatures, overwrite = TRUE)

Mapgene  = read.table("geneSymbol_to_geneID.txt", header = TRUE, sep = "\t")
Hgenes  = read.table("HumanGeneSymbol_to_geneID.txt", header = TRUE, sep = "\t")
Mgenes  = read.table("MouseGeneSymbol_to_geneID.txt", header = TRUE, sep = "\t")
usethis::use_data(Mapgene, overwrite = TRUE)
usethis::use_data(Hgenes, overwrite = TRUE)
usethis::use_data(Mgenes, overwrite = TRUE)

example_input = read.table("example_input", header=TRUE, row.names=1)
usethis::use_data(example_input, overwrite = TRUE)


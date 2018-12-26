setwd("data-raw")

Signatures <- RJSONIO::fromJSON("signatures")
usethis::use_data(Signatures)

Mapgene  = read.table("geneSymbol_to_geneID.txt", header = TRUE, sep = "\t")
Hgenes  = read.table("HumanGeneSymbol_to_geneID.txt", header = TRUE, sep = "\t")
Mgenes  = read.table("MouseGeneSymbol_to_geneID.txt", header = TRUE, sep = "\t")
usethis::use_data(Mapgene)
usethis::use_data(Hgenes)
usethis::use_data(Mgenes)

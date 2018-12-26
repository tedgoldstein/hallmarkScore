setwd("data-raw")
Signatures <- RJSONIO::fromJSON("signatures")
usethis::use_data(Signatures)


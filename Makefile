

all:
	R -e "library(devtools);library(roxygen2);document()"
	(cd data-raw; Rscript prepareSignatures.R)

push:
	git push 

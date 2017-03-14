# Download the GEOQuery package. Using bioconductor to download. 
# I had to update a tonne of packages to do this. 

## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("GEOquery")
library(GEOquery)
library(purrr)

data_meth <- getGEO("GSE85566", getGPL = FALSE)
yoo <- getGEOfile("GSE85568")

load(raw_data/GSE85568)

gds <- getGEO(filename="GSE85568.soft.gz")

gds_eset <- GDS2eSet(gds@gsms[[1]]@dataTable)

names <- purrr:map(gds, )

p <- phenoData(gds@gsms[[1]]@dataTable)


# OK so the above command does not work on my computer. I keep getting and error. 
# Error in download.file(myurl, destfile, mode = mode, quiet = TRUE, method = getOption("download.file.method.GEOquery")) : 
# cannot open URL 'http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?targ=self&acc=GPL11154&form=text&view=full'

counttable <- read.table("~/Downloads/GSE85567_RNASeq_normalizedcounts.txt")
dontknow <- read.table("~/Downloads/GSE85566_processed_meth_data_covarates.txt", header=TRUE, sep=",", row.names=1)
counttable <- t(counttable)

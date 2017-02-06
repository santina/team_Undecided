# Download the GEOQuery package. Using bioconductor to download. 
# I had to update a tonne of packages to do this. 

## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("GEOquery")
library(GEOquery)

data <- getGEO("GSE85568")
yoo <- getGEOfile("GSE85568")

# OK so the above command does not work on my computer. I keep getting and error. 
# Error in download.file(myurl, destfile, mode = mode, quiet = TRUE, method = getOption("download.file.method.GEOquery")) : 
# cannot open URL 'http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?targ=self&acc=GPL11154&form=text&view=full'
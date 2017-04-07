Differential Expression 3 Group
================
Arjun Baghela
3/29/2017

Our last step was [k-means clustering](https://github.com/STAT540-UBC/team_Undecided/blob/master/src/2_kmeans_clustering/Cluster.md), where we assigned each asthma patient to a group (Th2-high or Th2-low) depending on their expression in three genes.

In this document, we will perform differential expression analysis of the RNA-seq data between the three groups (control vs high, control vs low, and high vs low). The resulting three lists will be combined, and we'll filter out about 500 genes (deemed "interesting"), which will be passed to the next stage for constructing the correlation network and performing differential expression analysis.
Load the necessary packages.

``` r
library(edgeR)
```

    ## Loading required package: limma

``` r
library(limma)
library(tidyverse)
```

    ## Loading tidyverse: ggplot2
    ## Loading tidyverse: tibble
    ## Loading tidyverse: tidyr
    ## Loading tidyverse: readr
    ## Loading tidyverse: purrr
    ## Loading tidyverse: dplyr

    ## Conflicts with tidy packages ----------------------------------------------

    ## filter(): dplyr, stats
    ## lag():    dplyr, stats

``` r
library(magrittr)
```

    ## 
    ## Attaching package: 'magrittr'

    ## The following object is masked from 'package:purrr':
    ## 
    ##     set_names

    ## The following object is masked from 'package:tidyr':
    ## 
    ##     extract

Load in the normalized RNA-seq count data.

``` r
countdata <- read.table(file= "../../data/raw_data/rna_seq_data/GSE85567_RNASeq_normalizedcounts.txt", check.names = FALSE)
metadata <- read.csv(file= "../../data/processed_data/metaCluster.csv", row.names = 1)

metadata %<>% filter(ID %in% colnames(countdata)) # Remove metadata rows that do not have

metadata %>% group_by(Status) %>% tally() # See how many patients there are in each group
```

    ## # A tibble: 2 × 2
    ##    Status     n
    ##    <fctr> <int>
    ## 1  Asthma    57
    ## 2 Control    28

``` r
metadata %<>% arrange(cluster) # arrange by classes
metadata %>% group_by(cluster) %>% tally() # how many patients are there in each cluster- Control, Th2 high (1), Th2 low (2). 
```

    ## # A tibble: 3 × 2
    ##   cluster     n
    ##     <int> <int>
    ## 1       0    28
    ## 2       1    28
    ## 3       2    29

``` r
countdata <- countdata[,as.character(metadata$ID)] # Maintain the same order in both metadata and count data
colnames(countdata) == metadata$ID
```

    ##  [1] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
    ## [15] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
    ## [29] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
    ## [43] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
    ## [57] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
    ## [71] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
    ## [85] TRUE

``` r
metadata$cluster <- as.factor(metadata$cluster)
str(metadata)
```

    ## 'data.frame':    85 obs. of  9 variables:
    ##  $ ID                : Factor w/ 85 levels "0447_47e4","0610_40f0",..: 57 62 28 43 18 72 21 32 12 22 ...
    ##  $ Status            : Factor w/ 2 levels "Asthma","Control": 2 2 2 2 2 2 2 2 2 2 ...
    ##  $ Gender            : Factor w/ 2 levels "Female","Male": 1 2 2 1 2 1 1 1 1 1 ...
    ##  $ Age               : num  37 48 35 36 45 45 42 35 53 42 ...
    ##  $ Ethnicity         : Factor w/ 3 levels "AA","EA","Other": 1 1 1 1 1 1 1 1 1 1 ...
    ##  $ current_smoker    : Factor w/ 3 levels "","N","Y": 2 3 2 2 2 3 2 2 3 2 ...
    ##  $ Smoke_Ever        : Factor w/ 3 levels "","N","Y": 2 3 2 2 2 3 2 3 3 2 ...
    ##  $ smoke_pack_years_1: Factor w/ 16 levels "",".","0","0.25",..: 3 8 3 3 3 7 3 5 11 3 ...
    ##  $ cluster           : Factor w/ 3 levels "0","1","2": 1 1 1 1 1 1 1 1 1 1 ...

``` r
remove <- c("ENSG00000016490", "ENSG00000197632", "ENSG00000133110") # Remove the three genes we clustered on, to remove any biases. 
countdata <- countdata[!(rownames(countdata) %in% remove),]
```

We'll process the data some more.

``` r
DGElist <- DGEList(counts= countdata, group= metadata$cluster) # create DGEList to store data in
DGElist$samples$lib.size %>% min() # Find lib size
```

    ## [1] 18164788

``` r
DGElist$counts %>% nrow()
```

    ## [1] 16532

``` r
keep <- rowSums(cpm(DGElist)>0.275) >= 28 # Filter lowly expressed genes. Genes with less than 5 counts as determined by library with smallest size. 

DGElistFilt <- DGElist[keep, , keep.lib.sizes=FALSE] 
DGElistFilt$counts %>% nrow()
```

    ## [1] 15219

``` r
DGElistFiltNorm<- calcNormFactors(DGElistFilt) # calculate Norm factors
```

Now, let's perform the differential expression analysis, using edgeR.

``` r
design <- model.matrix(~0+group, data=DGElistFiltNorm$samples) #Create model matrix
design %>% head()
```

    ##           group0 group1 group2
    ## a7d1_4fec      1      0      0
    ## b3ae_412f      1      0      0
    ## 3d3b_40f5      1      0      0
    ## 78d1_4946      1      0      0
    ## 2800_4691      1      0      0
    ## dd18_4326      1      0      0

``` r
DGElistFiltNormDisp <- estimateDisp(DGElistFiltNorm, design) # Calculate dispersion 
fit <- glmFit(DGElistFiltNormDisp, design) #GLM fit

lrt <- glmLRT(fit, contrast=c(-1,1,0)) 
filtHighCon <- (topTags(lrt, n= Inf))$table %>% rownames_to_column(var="Gene") %>% filter(FDR <= .01)

lrt <- glmLRT(fit, contrast=c(-1,0,1))
filtLowCon <- (topTags(lrt, n= Inf))$table %>% rownames_to_column(var="Gene") %>% filter(FDR <= .25)

lrt <- glmLRT(fit, contrast=c(0,1,-1))
filtHighLow <- (topTags(lrt, n= Inf))$table %>% rownames_to_column(var="Gene") %>% filter(FDR <= .01)

union <- union(filtHighCon$Gene, filtLowCon$Gene) %>% union(filtHighLow$Gene) # Union of Important genes. 
union %>% length() # List is 571 genes to do differential co-expression analysis on now. 
```

    ## [1] 571

``` r
write.table(union, "../../data/processed_data/DEGene_ForDiffCoexpAnalysis.txt") # Write it to a table. 
```

Now, we can move onto the next stage, [differential methylation](https://github.com/STAT540-UBC/team_Undecided/blob/master/src/1_data_inspect_%26_4_diff%20met/Cleaning_methylation_data.md#rna-seq-confirming-that-pcs-do-not-correlate-with-covariates), which is necessary for weighting the correlation network we'll generate later.

##### DONT LOOK AT THIS.

PRODUCES THE SAME RESULTS AS ABOVE.

``` r
# model.matrix(~cluster*Gender*Age*current_smoker, metadata) %>% colnames()

design <- model.matrix(~group, data=DGElistFiltNorm$samples)

colnames(design) <- c("(Intercept)", "Th2 High", "Th2 Low")

DGElistFiltNormDisp <- estimateDisp(DGElistFiltNorm, design)
plotBCV(DGElistFiltNormDisp)
```

``` r
fit <- glmFit(DGElistFiltNormDisp, design)
lrtHigh_Con <- glmLRT(fit, coef = 2)
lrtLow_Con <- glmLRT(fit, coef = 3)

DEHigh_Con <- (topTags(lrtHigh_Con, n= Inf))$table
DELow_Con <- (topTags(lrtLow_Con, n= Inf))$table

DEHigh_ConFilt <- DEHigh_Con %>% rownames_to_column(var="Gene") %>% filter(FDR <= .25)
DELow_ConFilt <- DELow_Con %>% rownames_to_column(var="Gene") %>% filter(FDR <= .25)

DEHigh_ConFilt %>% nrow()
DELow_ConFilt %>% nrow()

intersect(DEHigh_ConFilt$Gene,DELow_ConFilt$Gene) %>% length()
finalGeneList <- union(DEHigh_ConFilt$Gene,DELow_ConFilt$Gene)
```

``` r
DGElistFiltNorm2 <- DGElistFiltNorm
DGElistFiltNorm2$samples$group <- relevel(DGElistFiltNorm$samples$group, ref="2")

design2 <- model.matrix(~group, data=DGElistFiltNorm2$samples)

DGElistFiltNorm2Disp <- estimateDisp(DGElistFiltNorm2, design2)

fit2 <- glmFit(DGElistFiltNorm2Disp, design2)
lrt2 <-  glmLRT(fit2, coef=3)

DEHigh_Low <- (topTags(lrt2, n= Inf))$table

DEHigh_LowFilt <- DEHigh_Low %>% rownames_to_column(var="Gene") %>% filter(FDR <= .25)
DEHigh_LowFilt %>% nrow()

intersect(finalGeneList, DEHigh_LowFilt$Gene) %>% length()
```

Preliminary dina analysis
================
Eric Chu
2017-03-16

In this document, I have done some preliminary/exploratory in differential network analysis.

1.  I find genes that are differentially expressed (between control and asthma patients).
2.  I assess how these genes have differential connectivity between the correlation networks of two groups, using the dna package.
3.  I also sampled matching number of random genes to see if these are less likely to have differential connectivity.
4.  I ran dna on the two asthma endotype groups.... just to see if we observe any genes with differential connectivity.

Key discussions from this preliminary analysis 
* We're confident that differential network analysis will be informative of the disease mechanism. 
* We need to formalize the way we narrow down the number of genes to look at. We will do this by: 
+ looking at biological pathways associated with the highly differentially expressed genes 
+ finding other genes that are associated with these identified pathways; conduct literature review for validation 
+ these genes will be taken as genes important for asthma 
* In addition to dna, we will also look at differential coexpression as previously suggested by Amrit - also a technique of differential network analysis (not demonstrated in this markdown doc)

``` r
library(tidyverse)
```

    ## Warning: package 'tidyverse' was built under R version 3.3.2

    ## Loading tidyverse: ggplot2
    ## Loading tidyverse: tibble
    ## Loading tidyverse: tidyr
    ## Loading tidyverse: readr
    ## Loading tidyverse: purrr
    ## Loading tidyverse: dplyr

    ## Warning: package 'ggplot2' was built under R version 3.3.2

    ## Warning: package 'tidyr' was built under R version 3.3.2

    ## Conflicts with tidy packages ----------------------------------------------

    ## filter(): dplyr, stats
    ## lag():    dplyr, stats

``` r
library(dplyr)
# packages
library(dna)
```

    ## dna 1.1-1 loaded

``` r
library(limma)
```

    ## Warning: package 'limma' was built under R version 3.3.3

``` r
# ----- functions -----

firstRowAsName <- function(tibble) {
  if (nrow(tibble) >= 2) {
    names(tibble) <- tibble %>% slice(1) %>% unlist(use.names = FALSE) %>% as.character()
    tibble[2:nrow(tibble),]
  } else {
    tibble
  }
}

transposeTibble <- function(tibble) {
  tibble %>% as.data.frame(stringsAsFactors = FALSE) %>% t() %>% 
    as.data.frame(stringsAsFactors = FALSE) %>% 
    firstRowAsName() %>% rownames_to_column() %>% as_tibble()
}

prepareTransposedCounts <- function(tibble) {
  tibbleT <- transposeTibble(tibble) %>% select(sample = rowname, everything())
  colNames <- tibbleT %>% colnames()
  geneNames <- colNames[2:length(colNames)]
  for (currGeneName in geneNames) {
    tibbleT[[currGeneName]] <- as.numeric(tibbleT[[currGeneName]])
  }
  return(tibbleT)
}
```

First, load the raw data. Also construct a transposed count matrix for easy maniplation. (dna requires genes to be columns)

``` r
# setwd("/Users/ericchu/ws/team_Undecided/")

# raw data
rawMetadata <- read_csv("../Raw_Data/GSE85566_metadata.txt")
```

    ## Warning: Missing column names filled in: 'X1' [1]

    ## Parsed with column specification:
    ## cols(
    ##   X1 = col_integer(),
    ##   ID = col_character(),
    ##   Status = col_character(),
    ##   Gender = col_character(),
    ##   Age = col_double(),
    ##   Ethnicity = col_character(),
    ##   current_smoker = col_character(),
    ##   Smoke_Ever = col_character(),
    ##   smoke_pack_years_1 = col_character()
    ## )

``` r
rawCounts <- read.table("../Raw_Data/GSE85567_RNASeq_normalizedcounts.txt", check.names = FALSE) %>% 
  rownames_to_column() %>% as_tibble() %>% 
  select(geneId = rowname, everything())

rawCountsT <- rawCounts %>% prepareTransposedCounts()
```

Differential expression analysis using Limma
--------------------------------------------

Below, we find genes that are differentially expressed between the controls and asthmatics. The hypothesis is that differentially expressed genes would play some role in the disease etiology and therefore are "important" genes to look at. We hope that this analysis will help us identify genes to look at in the subsequent differential network analyses.

``` r
# find differentially expressed genes using limma
sampleMetadata <- rawMetadata %>% filter(ID %in% colnames(rawCounts))
sampleMetadata$Status <- sampleMetadata$Status %>% as.factor() %>% relevel(ref = "Control")
sampleMetadata <- sampleMetadata %>% arrange(Status)
design <- model.matrix(~Status, sampleMetadata)

# prepare expression matrix
expressionMatrix <- rawCounts[sampleMetadata$ID] %>% as.data.frame()
rownames(expressionMatrix) <- rawCounts$geneId

fit <- lmFit(expressionMatrix, design)
fit <- eBayes(fit)

diffGenes <- fit %>% topTable(coef = "StatusAsthma", 
                              number = Inf, adjust.method = "fdr",
                              p.value = 0.1,
                              sort.by = "p")


# only 6 genes are differentially expressed with FDR set at 0.1. Though this is obviously a very stringent requirement!
diffGenes
```

    ##                      logFC   AveExpr         t      P.Value   adj.P.Val
    ## ENSG00000109956  -42.76629  59.00000 -5.483558 4.340549e-07 0.007177099
    ## ENSG00000133121   41.02130 137.29412  4.898817 4.655215e-06 0.038486988
    ## ENSG00000197629 -140.12845 272.35294 -4.736757 8.784643e-06 0.048418024
    ## ENSG00000138829   34.55952 102.28235  4.415791 2.990608e-05 0.098056934
    ## ENSG00000172348  -24.38033  48.32941 -4.403895 3.126776e-05 0.098056934
    ## ENSG00000172765 -139.09712 815.25882 -4.369267 3.558159e-05 0.098056934
    ##                          B
    ## ENSG00000109956 -0.8678741
    ## ENSG00000133121 -1.5092984
    ## ENSG00000197629 -1.6844595
    ## ENSG00000138829 -2.0262673
    ## ENSG00000172348 -2.0387830
    ## ENSG00000172765 -2.0751468

Below I run the dna package to genes that have differential connectivity between healthy controls and asthma patients.

Differentially expressed genes are must more likely to also exhibit differential connectivity!

``` r
# losen the fdr requirement to get more genes; 51 genes
diffGenes <- fit %>% topTable(coef = "StatusAsthma", 
                              number = Inf, adjust.method = "fdr",
                              p.value = 0.2,
                              sort.by = "p")

diffGenes %>% nrow()
```

    ## [1] 51

``` r
# separate into 2 groups - control & asthma
asthmaSamplesMetadata <- rawMetadata %>% filter(Status == "Asthma")
controlSamplesMetadata <- rawMetadata %>% filter(Status == "Control")


# try running dna with differentially expressed genes

# find genes with differential connectivity! - individual genes!

asthmaExpression <- rawCountsT %>% filter(sample %in% asthmaSamplesMetadata$ID)
asthmaExpression <- asthmaExpression[rownames(diffGenes)] %>% as.data.frame()
controlExpression <- rawCountsT %>% filter(sample %in% controlSamplesMetadata$ID)
controlExpression <- controlExpression[rownames(diffGenes)] %>% as.data.frame()
invisible(
  geneLevelResult <- test.individual.genes(controlExpression, asthmaExpression,
                                           scores = "cor", distance = "abs",
                                           num.permutations = 1000)
)
```

``` r
geneLevelResult %>% summary()
```

    ## Tests for differential connectivity of individual genes
    ## 
    ## 2 genes are significant at level 0.001
    ## 6 genes are significant at level 0.005
    ## 8 genes are significant at level 0.01
    ## 18 genes are significant at level 0.05

``` r
# now try running it with random set of genes and see if we observe differences!

randomGenes <- rawCounts$geneId %>% sample(51)
asthmaExpressionRandom <- rawCountsT %>% filter(sample %in% asthmaSamplesMetadata$ID)
asthmaExpressionRandom <- asthmaExpressionRandom[randomGenes] %>% as.data.frame()
controlExpressionRandom <- rawCountsT %>% filter(sample %in% controlSamplesMetadata$ID)
controlExpressionRandom <- controlExpressionRandom[randomGenes] %>% as.data.frame()
invisible(
  geneLevelResultRandom <- test.individual.genes(controlExpressionRandom, asthmaExpressionRandom,
                                                 scores = "cor", distance = "abs",
                                                 num.permutations = 1000)
)
```

``` r
geneLevelResultRandom %>% summary()
```

    ## Tests for differential connectivity of individual genes
    ## 
    ## 0 genes are significant at level 0.001
    ## 0 genes are significant at level 0.005
    ## 0 genes are significant at level 0.01
    ## 3 genes are significant at level 0.05

Next, I try running it on the patient clusters we have.

``` r
# ----- try running dna on different samples using differentially expressed genes -----


thHighSamples <- read_table("data/th2_high_samples.txt")
thLowSamples <- read_table("data/th2_low_samples.txt")
thHighExpression <- rawCountsT %>% filter(sample %in% thHighSamples$sample)
thHighExpression <- thHighExpression[rownames(diffGenes)] %>% as.data.frame()
thLowExpression <- rawCountsT %>% filter(sample %in% thLowSamples$sample)
thLowExpression <- thLowExpression[rownames(diffGenes)] %>% as.data.frame()
thResultGene <- test.individual.genes(thLowExpression, thHighExpression,
                                      scores = "cor", distance = "abs",
                                      num.permutations = 1000)
```

``` r
thResultGene %>% summary()
```

    ## Tests for differential connectivity of individual genes
    ## 
    ## 0 genes are significant at level 0.001
    ## 0 genes are significant at level 0.005
    ## 1 genes are significant at level 0.01
    ## 4 genes are significant at level 0.05

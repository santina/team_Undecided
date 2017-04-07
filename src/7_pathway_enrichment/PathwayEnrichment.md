Pathway Enrichment
================
Arjun Baghela
4/3/2017

Load Packages. I will be using SIGORA. It looks at pairs of genes that, as a combination, are specific to a single pathway.

``` r
library(tidyverse)
```

    ## Warning: package 'tidyverse' was built under R version 3.2.5

    ## Loading tidyverse: ggplot2
    ## Loading tidyverse: tibble
    ## Loading tidyverse: tidyr
    ## Loading tidyverse: readr
    ## Loading tidyverse: purrr
    ## Loading tidyverse: dplyr

    ## Warning: package 'ggplot2' was built under R version 3.2.5

    ## Warning: package 'tibble' was built under R version 3.2.5

    ## Warning: package 'tidyr' was built under R version 3.2.5

    ## Warning: package 'readr' was built under R version 3.2.5

    ## Warning: package 'purrr' was built under R version 3.2.5

    ## Warning: package 'dplyr' was built under R version 3.2.5

    ## Conflicts with tidy packages ----------------------------------------------

    ## filter(): dplyr, stats
    ## lag():    dplyr, stats

``` r
library(sigora)
```

    ## Warning: package 'sigora' was built under R version 3.2.5

Read in the data. That is, the gene pair list that came out of the gene pair analysis.

``` r
conHighRaw <- readRDS("../../data/processed_data/controlHighResults.rds")
conLowRaw <- readRDS("../../data/processed_data/controlLowResults.rds")
highLowRaw <- readRDS("../../data/processed_data/highLowResults.rds")

# Here, i am just getting gene pairs with an FDR of 0, But this is also arbtirary. 

conHigh <- conHighRaw[conHighRaw$fdr==0,]$gene_pair %>% strsplit(".", fixed= TRUE) %>% unlist() %>% unique()

conLow <- conLowRaw[conLowRaw$fdr==0,]$gene_pair %>% strsplit(".", fixed= TRUE) %>% unlist() %>% unique()

highLow <- highLowRaw[highLowRaw$fdr==0,]$gene_pair %>% strsplit(".", fixed= TRUE) %>% unlist() %>% unique()

conHigh %>% length() #532 genes from the pairs are unique with a fdr of 0. This is interesting, considering there were only 571 genes put into the diff network. 
```

    ## [1] 532

``` r
conLow %>% length()
```

    ## [1] 391

``` r
highLow %>% length()
```

    ## [1] 383

``` r
# I can also see which genes have the highest connectivity given a threshold 
conHighRaw[conHighRaw$fdr==0,]$gene_pair %>% strsplit(".", fixed= TRUE) %>% unlist() %>% table() %>% sort() %>% tail(10)
```

    ## .
    ## ENSG00000107821 ENSG00000124507 ENSG00000129467 ENSG00000204149 
    ##              40              40              42              44 
    ## ENSG00000007129 ENSG00000168329 ENSG00000171097 ENSG00000146205 
    ##              45              45              50              52 
    ## ENSG00000196296 ENSG00000143127 
    ##              57              59

``` r
conLowRaw[conLowRaw$fdr==0,]$gene_pair %>% strsplit(".", fixed= TRUE) %>% unlist() %>% table() %>% sort() %>% tail(10)
```

    ## .
    ## ENSG00000168542 ENSG00000110619 ENSG00000120457 ENSG00000166391 
    ##              13              15              15              17 
    ## ENSG00000096060 ENSG00000174370 ENSG00000196296 ENSG00000243056 
    ##              18              18              19              22 
    ## ENSG00000179593 ENSG00000141744 
    ##              25              29

``` r
highLowRaw[highLowRaw$fdr==0,]$gene_pair %>% strsplit(".", fixed= TRUE) %>% unlist() %>% table() %>% sort() %>% tail(10)
```

    ## .
    ## ENSG00000121335 ENSG00000167779 ENSG00000174059 ENSG00000012223 
    ##              13              13              13              17 
    ## ENSG00000204577 ENSG00000170373 ENSG00000007171 ENSG00000104044 
    ##              17              20              23              25 
    ## ENSG00000006606 ENSG00000156510 
    ##              27              34

The actual pathway enrichment is easy.

``` r
data("kegH")
sigoraResconHigh <- sigora(GPSrepo= kegH, level= 2, queryList = conHigh)
```

    ##    pathwy.id                               description   pvalues
    ## 38  hsa04380                Osteoclast differentiation 1.163e-38
    ## 24  hsa04070     Phosphatidylinositol signaling system 2.738e-17
    ## 44  hsa04611                       Platelet activation 5.903e-13
    ## 98  hsa05414                    Dilated cardiomyopathy 9.553e-10
    ## 64  hsa04970                        Salivary secretion 1.562e-09
    ## 80  hsa05146                                Amoebiasis 2.465e-09
    ## 65  hsa04972                      Pancreatic secretion 1.986e-06
    ## 16  hsa04012                    ErbB signaling pathway 4.840e-06
    ## 48  hsa04650 Natural killer cell mediated cytotoxicity 6.834e-06
    ## 23  hsa04062               Chemokine signaling pathway 7.887e-06
    ## 71  hsa05031                     Amphetamine addiction 2.946e-05
    ##    Bonferroni successes PathwaySize      N sample.size
    ## 38  3.256e-36     30.33      417.58 350891      596.45
    ## 24  7.666e-15     26.56     1623.80 350891      596.45
    ## 44  1.653e-10     36.32     4911.93 350891      596.45
    ## 98  2.675e-07     27.00     3807.25 350891      596.45
    ## 64  4.374e-07      9.00      274.01 350891      596.45
    ## 80  6.902e-07     16.81     1344.41 350891      596.45
    ## 65  5.561e-04     19.72     3056.27 350891      596.45
    ## 16  1.355e-03     14.46     1843.03 350891      596.45
    ## 48  1.914e-03     19.09     3336.26 350891      596.45
    ## 23  2.208e-03      7.64      408.85 350891      596.45
    ## 71  8.249e-03     21.85     4384.48 350891      596.45

``` r
sigoraResconLow <- sigora(GPSrepo= kegH, level= 2, queryList = conLow)
```

    ##    pathwy.id                description   pvalues Bonferroni successes
    ## 47  hsa04970         Salivary secretion 2.929e-12  8.201e-10      9.00
    ## 48  hsa04972       Pancreatic secretion 1.371e-09  3.839e-07     17.08
    ## 58  hsa05146                 Amoebiasis 2.272e-09  6.362e-07     12.67
    ## 29  hsa04380 Osteoclast differentiation 3.360e-09  9.408e-07      8.56
    ## 13  hsa04012     ErbB signaling pathway 4.361e-06  1.221e-03     10.62
    ## 33  hsa04611        Platelet activation 4.978e-06  1.394e-03     16.77
    ## 60  hsa05152               Tuberculosis 3.093e-05  8.660e-03      5.28
    ##    PathwaySize      N sample.size
    ## 47      274.01 350891      292.14
    ## 48     3056.27 350891      292.14
    ## 58     1344.41 350891      292.14
    ## 29      417.58 350891      292.14
    ## 13     1843.03 350891      292.14
    ## 33     4911.93 350891      292.14
    ## 60      420.28 350891      292.14

``` r
sigoraReshighLow <- sigora(GPSrepo= kegH, level= 2, queryList = highLow)
```

    ##    pathwy.id                           description   pvalues Bonferroni
    ## 23  hsa04380            Osteoclast differentiation 2.172e-32  6.082e-30
    ## 43  hsa04970                    Salivary secretion 1.794e-11  5.023e-09
    ## 70  hsa05414                Dilated cardiomyopathy 4.293e-09  1.202e-06
    ## 28  hsa04611                   Platelet activation 1.238e-08  3.466e-06
    ## 55  hsa05146                            Amoebiasis 2.173e-08  6.084e-06
    ## 44  hsa04972                  Pancreatic secretion 8.691e-07  2.433e-04
    ## 13  hsa04070 Phosphatidylinositol signaling system 1.247e-06  3.492e-04
    ## 8   hsa04012                ErbB signaling pathway 4.134e-06  1.158e-03
    ## 57  hsa05152                          Tuberculosis 5.577e-06  1.562e-03
    ##    successes PathwaySize      N sample.size
    ## 23     23.76      417.58 350891      358.38
    ## 43      9.00      274.01 350891      358.38
    ## 70     20.25     3807.25 350891      358.38
    ## 28     22.42     4911.93 350891      358.38
    ## 55     12.25     1344.41 350891      358.38
    ## 44     15.53     3056.27 350891      358.38
    ## 13     11.73     1623.80 350891      358.38
    ## 8      11.69     1843.03 350891      358.38
    ## 57      6.36      420.28 350891      358.38

``` r
sigoraResconHigh$summary_results %>% head() # This is an example of the output
```

    ##    pathwy.id                           description   pvalues Bonferroni
    ## 38  hsa04380            Osteoclast differentiation 1.163e-38  3.256e-36
    ## 24  hsa04070 Phosphatidylinositol signaling system 2.738e-17  7.666e-15
    ## 44  hsa04611                   Platelet activation 5.903e-13  1.653e-10
    ## 98  hsa05414                Dilated cardiomyopathy 9.553e-10  2.675e-07
    ## 64  hsa04970                    Salivary secretion 1.562e-09  4.374e-07
    ## 80  hsa05146                            Amoebiasis 2.465e-09  6.902e-07
    ##    successes PathwaySize      N sample.size
    ## 38     30.33      417.58 350891      596.45
    ## 24     26.56     1623.80 350891      596.45
    ## 44     36.32     4911.93 350891      596.45
    ## 98     27.00     3807.25 350891      596.45
    ## 64      9.00      274.01 350891      596.45
    ## 80     16.81     1344.41 350891      596.45

``` r
sigoraResconLow$summary_results %>% head()
```

    ##    pathwy.id                description   pvalues Bonferroni successes
    ## 47  hsa04970         Salivary secretion 2.929e-12  8.201e-10      9.00
    ## 48  hsa04972       Pancreatic secretion 1.371e-09  3.839e-07     17.08
    ## 58  hsa05146                 Amoebiasis 2.272e-09  6.362e-07     12.67
    ## 29  hsa04380 Osteoclast differentiation 3.360e-09  9.408e-07      8.56
    ## 13  hsa04012     ErbB signaling pathway 4.361e-06  1.221e-03     10.62
    ## 33  hsa04611        Platelet activation 4.978e-06  1.394e-03     16.77
    ##    PathwaySize      N sample.size
    ## 47      274.01 350891      292.14
    ## 48     3056.27 350891      292.14
    ## 58     1344.41 350891      292.14
    ## 29      417.58 350891      292.14
    ## 13     1843.03 350891      292.14
    ## 33     4911.93 350891      292.14

``` r
sigoraReshighLow$summary_results %>% head()
```

    ##    pathwy.id                description   pvalues Bonferroni successes
    ## 23  hsa04380 Osteoclast differentiation 2.172e-32  6.082e-30     23.76
    ## 43  hsa04970         Salivary secretion 1.794e-11  5.023e-09      9.00
    ## 70  hsa05414     Dilated cardiomyopathy 4.293e-09  1.202e-06     20.25
    ## 28  hsa04611        Platelet activation 1.238e-08  3.466e-06     22.42
    ## 55  hsa05146                 Amoebiasis 2.173e-08  6.084e-06     12.25
    ## 44  hsa04972       Pancreatic secretion 8.691e-07  2.433e-04     15.53
    ##    PathwaySize      N sample.size
    ## 23      417.58 350891      358.38
    ## 43      274.01 350891      358.38
    ## 70     3807.25 350891      358.38
    ## 28     4911.93 350891      358.38
    ## 55     1344.41 350891      358.38
    ## 44     3056.27 350891      358.38

``` r
setdiff(sigoraResconHigh$summary_results$description,sigoraResconLow$summary_results$description) # This is a way to see differences in pathway enrichment. 
```

    ##  [1] "Phosphatidylinositol signaling system"                   
    ##  [2] "Amphetamine addiction"                                   
    ##  [3] "Aldosterone synthesis and secretion"                     
    ##  [4] "TGF-beta signaling pathway"                              
    ##  [5] "Nicotinate and nicotinamide metabolism"                  
    ##  [6] "Amino sugar and nucleotide sugar metabolism"             
    ##  [7] "cAMP signaling pathway"                                  
    ##  [8] "Steroid biosynthesis"                                    
    ##  [9] "Folate biosynthesis"                                     
    ## [10] "Ras signaling pathway"                                   
    ## [11] "Rap1 signaling pathway"                                  
    ## [12] "Protein processing in endoplasmic reticulum"             
    ## [13] "ECM-receptor interaction"                                
    ## [14] "Signaling pathways regulating pluripotency of stem cells"
    ## [15] "Hematopoietic cell lineage"                              
    ## [16] "TNF signaling pathway"                                   
    ## [17] "Synaptic vesicle cycle"                                  
    ## [18] "Glutamatergic synapse"                                   
    ## [19] "Inflammatory mediator regulation of TRP channels"        
    ## [20] "Alcoholism"                                              
    ## [21] "Bacterial invasion of epithelial cells"                  
    ## [22] "Shigellosis"                                             
    ## [23] "Salmonella infection"                                    
    ## [24] "Malaria"                                                 
    ## [25] "Epstein-Barr virus infection"                            
    ## [26] "Viral carcinogenesis"                                    
    ## [27] "Rheumatoid arthritis"

We planned not do this, because we were not sure if it would be appropriate, because we are kind of losing the idea of looking at things in gene pairs. But it is definitely something we can explore more.

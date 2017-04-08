Pathway Enrichment
================
Arjun Baghela
4/3/2017

Load Packages. I will be using SIGORA. It looks at pairs of genes that, as a combination, are specific to a single pathway.

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
library(sigora)
```

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
sigoraResconHigh$summary_results %>% head(30) # This is an example of the output
```

    ##    pathwy.id                                         description   pvalues
    ## 38  hsa04380                          Osteoclast differentiation 1.163e-38
    ## 24  hsa04070               Phosphatidylinositol signaling system 2.738e-17
    ## 44  hsa04611                                 Platelet activation 5.903e-13
    ## 98  hsa05414                              Dilated cardiomyopathy 9.553e-10
    ## 64  hsa04970                                  Salivary secretion 1.562e-09
    ## 80  hsa05146                                          Amoebiasis 2.465e-09
    ## 65  hsa04972                                Pancreatic secretion 1.986e-06
    ## 16  hsa04012                              ErbB signaling pathway 4.840e-06
    ## 48  hsa04650           Natural killer cell mediated cytotoxicity 6.834e-06
    ## 23  hsa04062                         Chemokine signaling pathway 7.887e-06
    ## 71  hsa05031                               Amphetamine addiction 2.946e-05
    ## 32  hsa04145                                           Phagosome 9.336e-05
    ## 82  hsa05152                                        Tuberculosis 9.513e-05
    ## 22  hsa04060              Cytokine-cytokine receptor interaction 9.869e-05
    ## 90  hsa05205                             Proteoglycans in cancer 2.131e-04
    ## 92  hsa05214                                              Glioma 4.913e-04
    ## 51  hsa04666                    Fc gamma R-mediated phagocytosis 7.016e-04
    ## 7   hsa00350                                 Tyrosine metabolism 8.035e-04
    ## 49  hsa04660                   T cell receptor signaling pathway 2.365e-03
    ## 77  hsa05133                                           Pertussis 2.837e-03
    ## 3   hsa00130 Ubiquinone and other terpenoid-quinone biosynthesis 3.085e-03
    ## 93  hsa05215                                     Prostate cancer 3.356e-03
    ## 58  hsa04910                           Insulin signaling pathway 5.542e-03
    ## 45  hsa04621                 NOD-like receptor signaling pathway 5.713e-03
    ## 81  hsa05150                     Staphylococcus aureus infection 7.941e-03
    ## 79  hsa05145                                       Toxoplasmosis 8.287e-03
    ## 39  hsa04510                                      Focal adhesion 1.420e-02
    ## 83  hsa05164                                         Influenza A 3.585e-02
    ## 19  hsa04020                           Calcium signaling pathway 4.264e-02
    ## 74  hsa05130               Pathogenic Escherichia coli infection 4.841e-02
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
    ## 32  2.614e-02      8.45      815.48 350891      596.45
    ## 82  2.664e-02      6.36      420.28 350891      596.45
    ## 22  2.763e-02     18.66     3732.87 350891      596.45
    ## 90  5.967e-02     23.24     5815.32 350891      596.45
    ## 92  1.376e-01     10.00     1596.57 350891      596.45
    ## 51  1.964e-01      9.57     1380.67 350891      596.45
    ## 7   2.250e-01      3.87      105.13 350891      596.45
    ## 49  6.622e-01     15.34     3762.93 350891      596.45
    ## 77  7.944e-01      2.57       45.78 350891      596.45
    ## 3   8.638e-01      2.24       47.61 350891      596.45
    ## 93  9.397e-01     14.86     3525.17 350891      596.45
    ## 58  1.000e+00     22.49     7065.67 350891      596.45
    ## 45  1.000e+00      8.65     1557.34 350891      596.45
    ## 81  1.000e+00      2.56       78.32 350891      596.45
    ## 79  1.000e+00     16.54     4744.94 350891      596.45
    ## 39  1.000e+00      3.02      294.03 350891      596.45
    ## 83  1.000e+00     17.89     6134.15 350891      596.45
    ## 19  1.000e+00      3.50      451.55 350891      596.45
    ## 74  1.000e+00      6.19     1529.91 350891      596.45

``` r
sigoraResconLow$summary_results %>% head(30)
```

    ##    pathwy.id                                     description   pvalues
    ## 47  hsa04970                              Salivary secretion 2.929e-12
    ## 48  hsa04972                            Pancreatic secretion 1.371e-09
    ## 58  hsa05146                                      Amoebiasis 2.272e-09
    ## 29  hsa04380                      Osteoclast differentiation 3.360e-09
    ## 13  hsa04012                          ErbB signaling pathway 4.361e-06
    ## 33  hsa04611                             Platelet activation 4.978e-06
    ## 60  hsa05152                                    Tuberculosis 3.093e-05
    ## 6   hsa00350                             Tyrosine metabolism 1.004e-04
    ## 17  hsa04062                     Chemokine signaling pathway 4.142e-04
    ## 39  hsa04666                Fc gamma R-mediated phagocytosis 1.160e-03
    ## 37  hsa04660               T cell receptor signaling pathway 1.392e-03
    ## 55  hsa05130           Pathogenic Escherichia coli infection 1.938e-03
    ## 61  hsa05164                                     Influenza A 2.162e-03
    ## 57  hsa05145                                   Toxoplasmosis 2.370e-03
    ## 68  hsa05214                                          Glioma 2.396e-03
    ## 69  hsa05215                                 Prostate cancer 3.096e-03
    ## 73  hsa05414                          Dilated cardiomyopathy 5.072e-03
    ## 44  hsa04915                      Estrogen signaling pathway 9.000e-03
    ## 50  hsa04976                                  Bile secretion 2.270e-02
    ## 62  hsa05166                                HTLV-I infection 3.124e-02
    ## 56  hsa05133                                       Pertussis 3.757e-02
    ## 11  hsa00980    Metabolism of xenobiotics by cytochrome P450 4.263e-02
    ## 10  hsa00604 Glycosphingolipid biosynthesis - ganglio series 4.476e-02
    ## 59  hsa05150                 Staphylococcus aureus infection 6.288e-02
    ## 16  hsa04060          Cytokine-cytokine receptor interaction 9.385e-02
    ## 27  hsa04340                      Hedgehog signaling pathway 9.885e-02
    ## 42  hsa04910                       Insulin signaling pathway 1.384e-01
    ## 34  hsa04621             NOD-like receptor signaling pathway 1.415e-01
    ## 36  hsa04650       Natural killer cell mediated cytotoxicity 1.476e-01
    ## 66  hsa05205                         Proteoglycans in cancer 2.132e-01
    ##    Bonferroni successes PathwaySize      N sample.size
    ## 47  8.201e-10      9.00      274.01 350891      292.14
    ## 48  3.839e-07     17.08     3056.27 350891      292.14
    ## 58  6.362e-07     12.67     1344.41 350891      292.14
    ## 29  9.408e-07      8.56      417.58 350891      292.14
    ## 13  1.221e-03     10.62     1843.03 350891      292.14
    ## 33  1.394e-03     16.77     4911.93 350891      292.14
    ## 60  8.660e-03      5.28      420.28 350891      292.14
    ## 6   2.811e-02      3.87      105.13 350891      292.14
    ## 17  1.160e-01      4.02      408.85 350891      292.14
    ## 39  3.248e-01      6.09     1380.67 350891      292.14
    ## 37  3.898e-01     10.03     3762.93 350891      292.14
    ## 55  5.426e-01      6.29     1529.91 350891      292.14
    ## 61  6.054e-01     13.23     6134.15 350891      292.14
    ## 57  6.636e-01     11.86     4744.94 350891      292.14
    ## 68  6.709e-01      6.00     1596.57 350891      292.14
    ## 69  8.669e-01      9.98     3525.17 350891      292.14
    ## 73  1.000e+00      9.75     3807.25 350891      292.14
    ## 44  1.000e+00      3.18      506.11 350891      292.14
    ## 50  1.000e+00      2.34      276.62 350891      292.14
    ## 62  1.000e+00     11.01     6895.64 350891      292.14
    ## 56  1.000e+00      1.16       45.78 350891      292.14
    ## 11  1.000e+00      4.52     1562.19 350891      292.14
    ## 10  1.000e+00      1.00       55.04 350891      292.14
    ## 59  1.000e+00      1.28       78.32 350891      292.14
    ## 16  1.000e+00      6.00     3732.87 350891      292.14
    ## 27  1.000e+00      1.00      125.22 350891      292.14
    ## 42  1.000e+00      9.20     7065.67 350891      292.14
    ## 34  1.000e+00      3.50     1557.34 350891      292.14
    ## 36  1.000e+00      5.40     3336.26 350891      292.14
    ## 66  1.000e+00      7.86     5815.32 350891      292.14

``` r
sigoraReshighLow$summary_results %>% head(30)
```

    ##    pathwy.id                                            description
    ## 23  hsa04380                             Osteoclast differentiation
    ## 43  hsa04970                                     Salivary secretion
    ## 70  hsa05414                                 Dilated cardiomyopathy
    ## 28  hsa04611                                    Platelet activation
    ## 55  hsa05146                                             Amoebiasis
    ## 44  hsa04972                                   Pancreatic secretion
    ## 13  hsa04070                  Phosphatidylinositol signaling system
    ## 8   hsa04012                                 ErbB signaling pathway
    ## 57  hsa05152                                           Tuberculosis
    ## 49  hsa05031                                  Amphetamine addiction
    ## 31  hsa04660                      T cell receptor signaling pathway
    ## 67  hsa05215                                        Prostate cancer
    ## 56  hsa05150                        Staphylococcus aureus infection
    ## 18  hsa04144                                            Endocytosis
    ## 64  hsa05205                                Proteoglycans in cancer
    ## 30  hsa04650              Natural killer cell mediated cytotoxicity
    ## 27  hsa04610                    Complement and coagulation cascades
    ## 59  hsa05166                                       HTLV-I infection
    ## 24  hsa04510                                         Focal adhesion
    ## 52  hsa05133                                              Pertussis
    ## 11  hsa04060                 Cytokine-cytokine receptor interaction
    ## 3   hsa00350                                    Tyrosine metabolism
    ## 69  hsa05412 Arrhythmogenic right ventricular cardiomyopathy (ARVC)
    ## 1   hsa00010                           Glycolysis / Gluconeogenesis
    ## 25  hsa04512                               ECM-receptor interaction
    ## 38  hsa04910                              Insulin signaling pathway
    ## 21  hsa04260                             Cardiac muscle contraction
    ## 51  hsa05130                  Pathogenic Escherichia coli infection
    ## 29  hsa04621                    NOD-like receptor signaling pathway
    ## 54  hsa05145                                          Toxoplasmosis
    ##      pvalues Bonferroni successes PathwaySize      N sample.size
    ## 23 2.172e-32  6.082e-30     23.76      417.58 350891      358.38
    ## 43 1.794e-11  5.023e-09      9.00      274.01 350891      358.38
    ## 70 4.293e-09  1.202e-06     20.25     3807.25 350891      358.38
    ## 28 1.238e-08  3.466e-06     22.42     4911.93 350891      358.38
    ## 55 2.173e-08  6.084e-06     12.25     1344.41 350891      358.38
    ## 44 8.691e-07  2.433e-04     15.53     3056.27 350891      358.38
    ## 13 1.247e-06  3.492e-04     11.73     1623.80 350891      358.38
    ## 8  4.134e-06  1.158e-03     11.69     1843.03 350891      358.38
    ## 57 5.577e-06  1.562e-03      6.36      420.28 350891      358.38
    ## 49 5.876e-05  1.645e-02     15.65     4384.48 350891      358.38
    ## 31 5.876e-04  1.645e-01     12.09     3762.93 350891      358.38
    ## 67 1.167e-03  3.268e-01     11.51     3525.17 350891      358.38
    ## 56 2.961e-03  8.291e-01      2.56       78.32 350891      358.38
    ## 18 3.151e-03  8.823e-01     31.66    17885.99 350891      358.38
    ## 64 7.542e-03  1.000e+00     13.59     5815.32 350891      358.38
    ## 30 7.989e-03  1.000e+00      9.67     3336.26 350891      358.38
    ## 27 1.429e-02  1.000e+00      4.82      901.86 350891      358.38
    ## 59 2.657e-02  1.000e+00     13.57     6895.64 350891      358.38
    ## 24 3.678e-02  1.000e+00      2.40      294.03 350891      358.38
    ## 52 4.587e-02  1.000e+00      1.16       45.78 350891      358.38
    ## 11 9.049e-02  1.000e+00      7.55     3732.87 350891      358.38
    ## 3  1.017e-01  1.000e+00      1.58      105.13 350891      358.38
    ## 69 1.135e-01  1.000e+00      1.14      117.76 350891      358.38
    ## 1  1.440e-01  1.000e+00      2.63      652.50 350891      358.38
    ## 25 1.481e-01  1.000e+00      1.80      157.11 350891      358.38
    ## 38 1.893e-01  1.000e+00     10.14     7065.67 350891      358.38
    ## 21 2.045e-01  1.000e+00      5.43     3057.60 350891      358.38
    ## 51 2.063e-01  1.000e+00      3.28     1529.91 350891      358.38
    ## 29 2.134e-01  1.000e+00      3.57     1557.34 350891      358.38
    ## 54 2.138e-01  1.000e+00      7.52     4744.94 350891      358.38

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

The control vs High pathway enrichment and Th2 low vs high pathway enrichment seem to be quite similar, many of the same pathways show up. Similar pathways come up when doing the control vs Low, but in different a order. This may not be surprising, as we only worked with about 500 genes.

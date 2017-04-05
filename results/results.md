# Results

In this file, we walk through our results, as well as the main analyses conducted.  We'll follow the order specified by our pipeline:  

![pipeline](https://github.com/STAT540-UBC/team_Undecided/blob/master/results/figures/teamUndecided_Pipeline.png "Pipeline")
We'll also give a link to the source code, and general description of the inputs and outputs of each section.  

## 0. Data Inspection
[Source code](https://github.com/STAT540-UBC/team_Undecided/blob/master/src/1_data_inspection%26qc/Cleaning_methylation_data.md)  

*Input*: the normalized RNA-seq counts and methylation data, as downloaded from GEO.  
*Output*: none, as we deemed that no correction was necessary.  
To begin, we conducted exploratory analysis of both the RNA-seq and methylation data to check if further cleaning or correction was necessary.
Namely, we performed PCA using limma of the RNA-seq data with respect to the various covariates, to see if they cluster, and obtained some figures demonstrating the p-values associated between the covariates and the PCs.  

![RNA-seq](https://github.com/STAT540-UBC/team_Undecided/blob/master/src/1_data_inspection%26qc/Cleaning_methylation_data_files/figure-markdown_github/unnamed-chunk-3-1.png)
![methylation](https://github.com/STAT540-UBC/team_Undecided/blob/master/src/1_data_inspection%26qc/Cleaning_methylation_data_files/figure-markdown_github/unnamed-chunk-4-1.png)  

Looking at the first figure, we can see that certain PCs do correlate with East Asian ethnicity and smoking status in the RNA-seq data, but because PCs 16, 36, and 81 explain little variance, we decided that additional batch correction was unnecessary.  In the other figure, we see that significant PCs do correlate with gender, age, and ethnicity.  To deal with this, these variables will be controlled for when we perform differential methylation analysis.   

## 1. K-Means Clustering for Patient Differentiation
[Source Code](https://github.com/STAT540-UBC/team_Undecided/blob/master/src/2_kmeans_clustering/Cluster.Rmd)  
*Input*: the normalized RNA-seq counts.  
*Output*: a list that assigns every patient to a cluster, depending on if they were designated Th2-high, Th2-low, or simply belonged in the control group.  

In this step, as we lack clinical data, we decided to classify asthma as either Th2-high or Th2-low based on the expression of three genes: CLCA1, serpinB2, and periostin.  These genes are often used to differentiate between asthma types in literature.  To do this, we ran k-means clustering on the normalized RNA-seq counts of asthma patients with k = 2.  The result is as shown: 
![3dscatterplot](https://github.com/STAT540-UBC/team_Undecided/blob/master/results/figures/3D.png "3D Scatterplot")
This designated 28 patients as Th2-high, and 29 patients as Th2-low.  

## 2. Differential Expression Analysis for RNA-seq Data
[Source Code]()  
*Input*: the clusters associated with each patient from the previous step.  
*Output*: a list of "interesting genes." (The most differentially expressed genes between the three groups)   

After, we performed differential expression analysis with the RNA-seq data using edgeR to extract a list of "interesting genes" (or genes with high differential expression between our three groups) based on p-value and FDR.  This step was vital, because performing differential network analysis is computationally heavy, so running our entire set of genes through it is infeasible.  Using edgeR, we obtained a list of genes, ordered by FDR.  Setting a threshold of FDR <= 0.01 filtered our list down to 571, which we calculated takes around 30 hours to run for our differential network analysis step.

## 3. Differential Methylation Analysis
[Source Code]()  
*Input*: the cleaned methylation data from post-data inspection.  
*Output*: three lists (corresponding to control vs high, control vs low, and high vs low) of weights associated with each gene, depending how differentially methylated they are (the more, the stronger the weight).  

## 4. Correlation Network Generation and Differential Network Analysis with Permutation Testing
*Input*: the list of "interesting genes" from the DEA step, and the lists of gene weights from the differential methylation analysis.  
*Output*: three gene-pair lists, with the associated weight of each edge.  

## 5. Network and Gene-Pair Visualization
*Input*: the three gene pair lists, with edge weights.  
*Output*: A meta-network used for visualizing significant gene pairs, and plots of how the expression of various significant gene pairs differs.  

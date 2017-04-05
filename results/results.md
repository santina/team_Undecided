# Results

In this file, we walk through our results, as well as the main analyses conducted.  We'll follow the order specified by our pipeline:
![pipeline](https://github.com/STAT540-UBC/team_Undecided/blob/master/results/figures/teamUndecided_Pipeline.png "Pipeline")
We'll also give a link to the source code, and general description of the inputs and outputs of each section.  

## 0. Data Inspection
*Input*: the normalized RNA-seq counts and methylation data, as downloaded from GEO.  
*Output*: none, as we deemed that no correction was necessary.  
To begin, we conducted exploratory analysis of both the RNA-seq and methylation data to check if further cleaning or correction was necessary.
Namely, we performed PCA of the RNA-seq data with respect to the various covariates, to see if they cluster.  
(P-value figures goes here)  
Looking at our figures, we can see that 

## 1. K-Means Clustering for Patient Differentiation
[Source code](https://github.com/STAT540-UBC/team_Undecided/blob/master/src/2_kmeans_clustering/Cluster.Rmd)  
*Input*: the normalized RNA-seq counts.  
*Output*: a list that assigns every patient to a cluster, depending on if they were designated Th2-high, Th2-low, or simply belonged in the control group.  

Here, as we lack clinical data, we decided to classify asthma as either Th2-high or Th2-low based on the expression of three genes: CLCA1, serpinB2, and periostin.  These genes are often used to differentiate between asthma types in literature.  To do this, we ran k-means clustering on the normalized RNA-seq counts of asthma patients with k = 2.  The result is as shown: 
![3dscatterplot](https://github.com/STAT540-UBC/team_Undecided/blob/master/results/figures/3D.png "3D Scatterplot")
This designated 28 patients as Th2-high, and 29 patients as Th2-low.  

## 2. Differential Expression Analysis for RNA-seq Data
*Input*: the clusters associated with each patient from the previous step.  
*Output*: three lists of "interesting genes." (Corresponding to control vs high, control vs low, and high vs low)  


## 3. Differential Methylation Analysis
*Input*: the cleaned methylation data from post-data inspection.  
*Output*: three lists (corresponding to control vs high, control vs low, and high vs low) of weights associated with each gene, depending how differentially methylated they are (the more, the stronger the weight).  

## 4. Correlation Network Generation and Differential Network Analysis with Permutation Testing
*Input*: the list of "interesting genes" from the DEA step, and the lists of gene weights from the differential methylation analysis.  
*Output*: three gene pair lists, with the associated weight of each edge.  

## 5. Network and Gene-Pair Visualization
*Input*: the three gene pair lists, with edge weights.  
*Output*: A meta-network used for visualizing significant gene pairs, and plots of how the expression of various significant gene pairs differs.  

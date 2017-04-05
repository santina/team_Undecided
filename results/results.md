# Results

In this document, we walk through our results, as well as the main analyses conducted.  We'll follow the order specified by our pipeline:
![pipeline](https://github.com/STAT540-UBC/team_Undecided/blob/master/results/figures/teamUndecided_Pipeline.png "Pipeline")

## 0. Data Inspection
To begin, we conducted exploratory analysis of both the RNA-seq and methylation data to check if further cleaning or correction was necessary.
Namely, we performed PCA of the RNA-seq data with respect to certain covariates, to see if they cluster.  
(P-value figure goes here)  
Looking at our figures, we can see that 

## 1. K-Means Clustering for Patient Differentiation
Input: the normalized RNA-seq counts.  
Output: a list that assigns every patient to a cluster, depending on if they were designated Th2-high, Th2-low, or simply belonged in the control group.  

## 2. Differential Expression Analysis for RNA-seq Data
Input: the clusters associated with each patient from the previous step.  
Output: three lists of "interesting genes." (Corresponding to control vs high, control vs low, and high vs low)  


## 3. Differential Methylation Analysis
Input: the cleaned methylation data from post-data inspection.  
Output: three lists (corresponding to control vs high, control vs low, and high vs low) of weights associated with each gene, depending how differentially methylated they are (the more, the stronger the weight).  

## 4. Correlation Network Generation and Differential Network Analysis with Permutation Testing
Input: the list of "interesting genes" from the DEA step, and the lists of gene weights from the differential methylation analysis.  
Output: three gene pair lists, with the associated weight of each edge.  

## 5. Network and Gene-Pair Visualization
Input: the three gene pair lists, with edge weights.  
Output: A meta-network used for visualizing significant gene pairs, and plots of how the expression of various significant gene pairs differs.  

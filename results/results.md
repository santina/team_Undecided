# Results

In this file, we walk through our results, as well as the main analyses conducted.  We'll follow the order specified by our pipeline:  

![pipeline](https://github.com/STAT540-UBC/team_Undecided/blob/master/results/figures/teamUndecided_Pipeline.png "Pipeline")
We'll also give a link to the source code (which is actually the .md file), and general description of the inputs and outputs of each section.  

As a general introduction to the data, we have two main datasets, both downloaded off GEO with accession number [GSE85568](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE85568).  They are both stored in the raw_data subfolder.  The first one consists of normalized RNA-seq counts (though we didn't manage to find exactly what steps they performed), which is 16,535 gene entries by 85 patients (57 asthmatics, 28 controls).  There were 71 asthmatics and 41 controls in the original, but only 85 patients have both RNA-seq and methylation data available.  Therefore, the methylation dataset was 327,271 probes by 85 patients, and had also undergone some correction by the original authors.  

## 0. Data Inspection
[Source code](https://github.com/STAT540-UBC/team_Undecided/tree/master/src/1_data_inspect_%26_4_diff%20met/Cleaning_methylation_data.md)  
*Input*: the normalized RNA-seq counts and methylation data, as downloaded from GEO.  
*Output*: none, as we deemed that no further correction was necessary.  
To begin, we conducted exploratory analysis of both the RNA-seq and methylation data to check if further cleaning or correction was necessary.
Namely, we performed PCA using limma of the RNA-seq data with respect to the various covariates, to see if they cluster, and obtained some figures demonstrating the p-values associated between the covariates and the PCs.  

![rna-seq](https://github.com/STAT540-UBC/team_Undecided/blob/master/results/figures/rnaseq.png "RNA-seq")
![methylation](https://github.com/STAT540-UBC/team_Undecided/blob/master/results/figures/methylation.png "Methylation")

Looking at the first figure, we can see that certain PCs do correlate with East Asian ethnicity and smoking status in the RNA-seq data, but because PCs 16, 36, and 81 explain little variance, we decided that additional batch correction was unnecessary.  In the other figure, we see that significant PCs do correlate with gender, age, and ethnicity.  To deal with this, these variables will be controlled for when we perform differential methylation analysis.   

As we have a very large number of probes for the methylation data (around 320,000), we filtered out some probes found to be non-differentially methylated in most cells.  This removed around 40,000 probes, shrinking our dataset.  

## 1. K-Means Clustering for Patient Differentiation
[Source Code](https://github.com/STAT540-UBC/team_Undecided/blob/master/src/2_kmeans_clustering/Cluster.md)  
*Input*: the normalized RNA-seq counts.  
*Output*: a list that assigns every patient to a cluster, depending on if they were designated Th2-high, Th2-low, or simply belonged in the control group.  

In this step, as we lack clinical data, we decided to classify asthma as either Th2-high or Th2-low based on the expression of three genes: CLCA1, serpinB2, and periostin.  These genes are often used to differentiate between asthma types in literature.  To do this, we ran k-means clustering on the normalized RNA-seq counts of asthma patients with k = 2.  The result is as shown: 
![3dscatterplot](https://github.com/STAT540-UBC/team_Undecided/blob/master/results/figures/3D.png "3D Scatterplot")
This designated 28 patients as having higher expression in the three genes (Th2-high), and 29 patients having lower expression (Th2-low).  This supports the findings in literature that suggests they can be used as biomarkers for asthma endotypes.  

## 2. Differential Expression Analysis for RNA-seq Data
[Source Code](https://github.com/STAT540-UBC/team_Undecided/blob/master/src/3_differential_expression/DE_ThreeGroup.md)  
*Input*: the normalized RNA-seq count data; the clusters associated with each patient from the previous step.  
*Output*: a list of "interesting genes." (The most differentially expressed genes between the three groups.)   

After, we performed differential expression analysis with the RNA-seq data using edgeR to extract a list of "interesting genes" (or genes with high differential expression between our three groups) based on p-value and FDR.  This step was vital, because performing differential network analysis is computationally heavy, so running our entire set of genes through it is infeasible.  Using edgeR, we obtained a list of genes, ordered by FDR.  Setting a threshold of FDR <= 0.01 filtered our list down to 571, which we calculated takes around 30 hours to run for our differential network analysis step.  

## 3. Differential Methylation Analysis
[Source Code](https://github.com/STAT540-UBC/team_Undecided/blob/master/src/1_data_inspect_%26_4_diff%20met/Cleaning_methylation_data.md#assessment-of-differentially-methylated-sites)  
*Input*: the cleaned methylation data from post-data inspection; the clusters associated with each patient from the previous step.  
*Output*: three lists (corresponding to control vs high, control vs low, and high vs low) of weights associated with each gene, depending on how differentially methylated they are between groups (the more, the stronger the weight).  

We then mapped each probe to their gene loci, to obtain the methylation count associated with each gene for each patient.  We then performed differential methylation analysis using edgeR to obtain "weights" by taking -log2(FDR) associated with each gene for each comparison.  We also saw here that more differentially methylated CpGs were found when comparing asthmatics to controls then between asthma subtypes, as the "weights" in the Th2-high to Th2-low comparison tend to be low (all weights < 3), compared to the Th2-high to control comparison, with two weights > 30, or even the Th2-low to control comparison, with several weights > 13.  This indicates a shared methylation profile among asthmatics, regardless of subtype.  

## 4. Correlation Network Generation and Differential Network Analysis with Permutation Testing
[Source Code](https://github.com/STAT540-UBC/team_Undecided/blob/master/src/5_weighted_corr_net_%26_diff_analysis/differential_coexpression_analysis_demonstration.md) (though this a toy example)  
*Input*: the list of "interesting genes" from the DEA step, and the lists of gene weights from the differential methylation analysis.  
*Output*: three gene-pair lists, with the associated weight of each edge, the unadjusted p-value, the fdr, and p-value after Bonferroni correction.  

First, we added the weights as previously determined via differential methylation and added them to the expression values of each gene.  Our reasoning for doing this is that we wanted a way to select for gene pairs with both high differential co-expression and co-methylation, as this is our way of incorporating both methylation and RNA-seq data, as was one of our goals mentioned back in our proposal.  With the final network, high value edges would correspond to this.  There are some disadvantages, the main being that since both the differential expression AND differential methylation values are fused into one, it's hard to tell how much each contributes to the final edge weight.  We continue by constructing three correlation matrices for the three groups (control, high, low).  Then we calculate the differential correlations between them by taking the absolute difference of each pairwise comparison between any two groups, representing our observations of the coexpression changes when you move from one group to another (this took around 30 hours to run).  To try and increase confidence in our significant observations being actually significant, we also performed permutation tests (1000 iterations) for each gene pair.  The result of this was that some gene pairs were noted to have an FDR or p-value of 0, as those pairs were always more significant than ones produced by random draw.  The p-values of each gene pair was stored; below, we have the p-value distributions associated with each comparison (control vs high, control vs low, and high vs low).

### Control vs. Th2-High
![controlvhigh](https://github.com/STAT540-UBC/team_Undecided/blob/master/results/figures/figure1_control_high.png "Control vs High")
### Control vs. Th2-Low
![controlvlow](https://github.com/STAT540-UBC/team_Undecided/blob/master/results/figures/figure2_control_low.png "Control vs Low")
### Th2-High vs. Th2-Low
![highvlow](https://github.com/STAT540-UBC/team_Undecided/blob/master/results/figures/figure3_high_low.png "High vs Low")  

We can see from these figures that the difference between the co-expression profiles of Th2-high and control is much larger than that of Th2-high and low, or Th2-low and control.  This is understandable, because asthma patients should share similar profiles regardless of subtype, and we know from published literature that the expression profiles of Th2-low asthma patients are closer to control than Th2-high are to controls.  

The outputted lists were then sent to the next stage for visualization.  

## 5. Network and Gene-Pair Visualization
[Source Code](https://github.com/STAT540-UBC/team_Undecided/blob/master/src/6_network_visualization/networkFilter.md) and [Source Code](https://github.com/STAT540-UBC/team_Undecided/blob/master/src/5_weighted_corr_net_%26_diff_analysis/differential_coexpression_analysis_demonstration.md#permutation-distributions-for-specific-gene-pairs)  
(Note: the first source is the processing step so that Cytoscape can use the data, the second source is looking closely at individual gene pairs)  
*Input*: the three gene pair lists, with edge weights.  
*Output*: A meta-network used for visualizing significant gene pairs, and plots of how the expression of various significant gene pairs differs.  

After processing the inputs and filtering out enough genes so that they can actually be visualized (using FDR = 0 and in the top 500 edge weights for each group), we put them into Cytoscape.  We show some of our figures below.  
### Control vs. Th2-High
![controlvhigh](https://github.com/STAT540-UBC/team_Undecided/blob/master/results/figures/ctlVshigh_cytoscape_subnetwork_CLU.png "Control vs High")
### Control vs. Th2-Low
![controlvlow](https://github.com/STAT540-UBC/team_Undecided/blob/master/results/figures/ctlVsLow_cytoscape_subnetwork_PACSIN1.png "Control vs Low")
### Th2-High vs. Th2-Low
![highvlow](https://github.com/STAT540-UBC/team_Undecided/blob/master/results/figures/lowVsHigh_cytoscape_subnetwork.png "High vs Low")  

Here, the size of the node corresponds to the degree (number of edges it has), while the thickness of the edge corresponds to its weight.  Red mean that the particular edge is in the top 20% of edges, when ordered by weight.  Like this, we can see that some genes are more highly connected than anothers, which could suggest that some could be important players in asthma-related pathways.  One such gene that appeared as well-connected with high-weight edges (meaning high differential co-expression and/orl co-methylation) in the control vs. Th2-high comparison network was CLU, or clusterin, which has been implicated in other literature as an asthma-related biomarker.  Using similar requirements (genes with a decent number of higher-weight edges), we picked a couple of genes we deemed interesting to probe further into.  

We then looked at permutation distributions for these genes.  For instance, below we have
### Control vs. Th2-High: CLU vs. USP54 distribution
![cluvusp54](https://github.com/STAT540-UBC/team_Undecided/blob/master/results/figures/figure4_distribution_control_high_clu_usp54.png "CLU vs USP54")
### Control vs. Th2-Low: PACSIN vs. SPATS distribution
![pacsinvspats](https://github.com/STAT540-UBC/team_Undecided/blob/master/results/figures/figure6_distribution_control_low_pacsin1_spats2.png "PACSIN vs SPATS")
### Th2-High vs Th2-Low: GPSM3 vs. CLU distribution
![gpsm3vclu](https://github.com/STAT540-UBC/team_Undecided/blob/master/results/figures/figure8_distribution_high_low_gpsm3_clu.png "GPSM3 vs CLU")
For each of these, note how much higher the "absolute correlation with methylation weighting" are compared to the randomly drawn pairs.  In fact, for each of these, the original is higher than 100% of the random draws, therefore giving the pair a p-value of 0, and suggests that their significance wasn't by chance.  

Now we can look at the co-expression plots for those same genes!
### Control vs. Th2-High: CLU vs. USP54 co-expression
![cluvusp54](https://github.com/STAT540-UBC/team_Undecided/blob/master/results/figures/figure5_expression_correlations_control_high_clu_usp54.png "CLU vs USP54")
### Control vs. Th2-Low: PACSIN vs. SPATS co-expression
![pacsinvspats](https://github.com/STAT540-UBC/team_Undecided/blob/master/results/figures/figure7_expression_correlations_control_low_pacsin1_spats2.png "PACSIN vs SPATS")
### Th2-High vs Th2-Low: GPSM3 vs. CLU co-expression
![gpsm3vclu](https://github.com/STAT540-UBC/team_Undecided/blob/master/results/figures/figure9_expression_correlations_high_low_gpsm3_clu.png "GPSM3 vs CLU")  
This gives a great picture of how the weighted differential expression of one gene changes with regards to the other gene.

## 6. Pathway Enrichment
[Source Code](https://github.com/STAT540-UBC/team_Undecided/blob/master/src/7_pathway_enrichment/PathwayEnrichment.md)  
*Input*: the three gene pair lists, with edge weights.  
*Output*: list of significant pathways, obtained by running our genes through KEGG.  


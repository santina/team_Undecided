## Links

Please see here for links to .md files containing details of our progress / preliminary analyses that we have done.

* [Aim 1 - PCA explorative analysis](https://github.com/STAT540-UBC/team_Undecided/blob/master/Emma/Quality_Control_RNAseq_Methylation.md)
* [Aim 2 - Patient clustering using k-means](https://github.com/STAT540-UBC/team_Undecided/blob/master/Allison_Scripts/Cluster.md)
* [Aim 3 - WGCNA modules analysis - Th2High](https://github.com/STAT540-UBC/team_Undecided/blob/master/Arjun_Scripts/Th2HighWGCNA.pdf)
* [Aim 3 - WGCNA modules analysis - Th2Low](https://github.com/STAT540-UBC/team_Undecided/blob/master/Arjun_Scripts/Th2LowWGCNA.pdf)
* [Aim 3 - WGCNA modules analysis - description and code](https://github.com/STAT540-UBC/team_Undecided/blob/master/Arjun_Scripts/WGCNA.Rmd)


## Overview

While major components of our analysis pipeline have remained the same from our initial propsal, we have made a small number of modifications, mostly surrounding the details of our methods. These ammendments have been informed by the progress of exploratory investigations as well as consultation with the course instructors. We have given a high level description of these changes and our current progress in this document. Links to individual .md files containing more details have been linked.

Our datasets have remained the same; we are analyzing two data sets in conjunction: a gene expression dataset and a methylation dataset. The data were collected from airway epithelial cells of asthma patients and healthy controls. 

Differential network analysis refers to the characterization of differences between two gene networks. These differences can be operationalized by, for example, assessing differences between correlation values, partial correlations, or connectivity to other genes in the networks. 

The main theme of our project is to investigate whether mechanisms of asthma can be elucidated by leveraging differential network analysis techniques. It has been previously established that two subtypes types of asthma can be delineated by the expression values of three genes: CLCA1, periostin, and serpinB2 (see our initial proposal for details). The two subtypes of asthma have critical differences in their clinical manifestations and responses to treatment. Differences between these two groups as well as with healthy controls in terms of gene regulation have important implications for the underlying diseae mechanism. Driven by the hypothesis that differences in gene correlation networks (built using expression and methylation data) of the different disease groups can lend insight into the disease etiology and potential thearpeutic targets, we aim to characterize these differences using differential network analysis methodologies. 

Our approach can be largely separated into 4 aims: 1. Data preprocessing, 2. Patient clustering, 3. Hypothesis testing, 4. Biological interpretations. We give a brief description of our aims in this document with links to more detailed .md files as well as plots that may be of relevance. This pipeline has been reorganized for more clarity. 


## Aim 1: Data preprocessing

In this step, our original goal was to ensure that the data sets were properly prepared for analysis. This included 
* aligning RNA reads in order to estimate expression values
* normalizing count matrices
* performing PCA to determine and potentially eliminate batch effects
* mapping methylation sites to genes


### Read alignment and count normalization

In our original proposal, we had planned to use STAR to align reads, subsequently estimating transcript abundances and performing any other normalization necessary. 

This step is no longer needed, as the authors provide a normalized count table that has already undergone alignment to the hg19 genome, read length normalization and removal of covariates through RUVSeq. 

Therefore, we will be using the data sets as they are, without any further modifications. This will give us a little bit more bandwith to focus on hypothesis testing (which is, of course, much more exciting!). 


### Exploratory PCA

[Aim 1 - PCA explorative analysis](https://github.com/STAT540-UBC/team_Undecided/blob/master/Emma/Quality_Control_RNAseq_Methylation.md)

Initial proposal: implement Hidden Covariates with Prior algorithm (Mostafavi, et al, written in matlab) to adjust for known (GC bias, age, sex) and unknown (comorbidities, batch effects) covariates.

When PCA is used to visualize the first two PCs of the data, we do not see any clusters based on known covariates. This means we don’t need to adjust for known covariates. We cannot use the HCP algorithm, because we do not have the unnormalized count table (we have the raw reads, however re-creating the count table from these reads would take extensive CPU time, which would detract from time allotted to other parts of the project). However, as can be seen in the results section, no visible clusters can be seen through PCA of our data. 


### Methylation sites preprocessing & mapping

The exact normalization methods the data deposited on GEO has undergone are not clear, so we are clarifying with the authors. 

The methylation data deposited on GEO contains probes that map to unique genes, so we do not need to remove promiscuous probes. If they have already used the SWAN algorithm to correct for probe bias, we will not do any further adjustment. PCA of all methylation data suggests that expression does not depend on known covariates. 

It is still possible for multiple probes to map to the same gene. In this case, we will obtain a summary value by taking, likely, the maximum value. Other options include taking the mean or the median. There is no clear consensus on the best approach in the current literature.  


## Aim 2: Patient clustering

[Aim 2 - Patient clustering using k-means](https://github.com/STAT540-UBC/team_Undecided/blob/master/Allison_Scripts/Cluster.md)

This step has remainted largely unchanged since our initial proposal. 

We are using k-means, since we have a predefined number of clusters (chose 2). CLCA1 shows very clear high/level expression clusters, the other genes less so. We have 6 Th2-high and 51 Th2-low (though some in the low-group are better choices than others). 

We're very excited to be able to find clusters corresponding to observations reported in the literature! This, of course, results in an imbalanced dataset. We will likely perform bootstrapping to resolve this issue. 


## Aim 3: Hypothesis testing

This aim is the key step in our analysis where we investigate and test our hypothesis formally. It consists of 3 distinct steps: gene filtering, network construction and differential network analysis.

To reiterate, our hypothesis is that molecular differences between healthy controls and the two subgroups of asthma patients can be characterized by leveraging differential network analyses on gene expression and methylation data. In combination with aim 4: biological interpretation, we attempt to elucidate the disease etiology of asthma and potentially learn about therapeutic targets. 

### Gene filtering

The number of genes must be filtered into something manageable. The human genome consists of ~20,000 genes. Differentially network analysis approaches typically require computing connectivity information at some global level for assessing individual genes and edges. For example, in the dna package, the difference between one single gene between two networks is assessed by taking the difference between their connectivity to the rest of the network. Dingo computes partial correlations, also by taking global connectivity into account. This entails that these algorithms can take a long time to compute. It is, therefore, important that we narrow down to the genes that we're interested in. 

There are a number of approaches to doing this. 

We can conduct differential expression analysis to find the genes that are differentially expressed between controls and asthmatics. This list of genes would likely contain genes that play important roles in the disease mechanism. However, our preliminary investigaton using limma suggest that very few genes are differentially expressed; at a significant value of p-value = 0.05, only 6 genes were found. This may not give us enough power to detect differences.

A more promising approach is to look at genes that are associated with relevant biological pathways. These pathways may be identified by quickly checking the Gene Ontology for the differentially expressed genes. More literature review will help us validate the role of these pathways in asthma. Once identified, we will focus on these genes in our subsequent analyses. 

We will attempt to maximize the number of genes we can look at given our resources, maybe ~100 genes or so. 


### Network construction

From aim 2: patient clustering, we have obtained two clusters of patients corresponding to two different subgroups of asthmatics. In combination with the control group, we have three groups in total. 

Separate correlation networks will be constructed for each group, using only the genes identified in the previous step. Each gene will be present as two nodes in the graph as it contains information from the expression as well as methylation data (methylation sites have been mapped to genes in aim 1); these correlation networks will contain hetergeous data types. The distance measure will simply be pearson's correlation in order to avoid biases introduced by the different scales in the two data types. As a result, the networks will be represented as correlation matrices with genes on both axes and each entry containing the correlation value of the corresponding two genes. 

As part of network construction, we also plan to cluster genes into various modules. While these modules are interesting in their own right, they can help us perform different network analysis. For example, investigating module preservation statistics.  

We have already done some preliminary exploration of performing WGCNA analysis. This is done using the entire expression dataset. See preliminary results: 

* [Aim 3 - WGCNA modules analysis - Th2High](https://github.com/STAT540-UBC/team_Undecided/blob/master/Arjun_Scripts/Th2HighWGCNA.pdf)
* [Aim 3 - WGCNA modules analysis - Th2Low](https://github.com/STAT540-UBC/team_Undecided/blob/master/Arjun_Scripts/Th2LowWGCNA.pdf)
* [Aim 3 - WGCNA modules analysis - description and code](https://github.com/STAT540-UBC/team_Undecided/blob/master/Arjun_Scripts/WGCNA.Rmd)


### Differential network analysis

Once we have the networks, we plan to carry out differential network analysis. 

So far, we have done some investigation using the dna package. This package assesses genes that are "differentially connected" to the rest of the network. In the preliminary analysis, we identified genes that are differentially expressed between control and patient groups and then used dna to find genes that have differential connectivity. It appears that in comparison to random genes, differentially expressed genes are more likely to also exhibit differential connectivity. The dna package is one way we will use to characterize the differences. 


In addition to dna, we will also develop our own differential netowrk analysis approach, specifically looking at differential correlations. In code, the networks will be represented as simply correlation matrices with genes on both axes and each entry containing the correlation value of the corresponding two genes. The networks constructed using the three groups will have identical dimensions on both axes. Therefore, a difference value can be easily computed. By bootstrapping, we can obtain the distribution of these differences and subsquently deriving p-values to assess these differences in correlation. 

### Multiple testing

Since this procedure can only be done on pairs of networks, we would only be able to assess pairwise comparisons. Therefore, to detect significant differences, we will adjust for multiple testing. Our choice will likely be FDR. 

These differential correlations will then inform our pathway analysis to arrive at appropriate biological interpretations.


## Aim 4: Biological interpretations

### Pathway analysis

The result of differential network analyses be a list of genes that are likely to play important roles in the asthma disease etiology. We will then perform pathway enrichment analyses to identify pathways corresponding to these genes. Informed by the existing literature, we can then attempt to interpret our findings in the context of biology. Pathway analyses can be done by accessing Gene Ontology, for example, and look at functions that are most represented by the genes we have identified. 


### Drug target identification (stretch goal)

It is our stretch (actually uber stretch) goal to look at potential therapeutic targets and drug availability. DrugBank.ca has an API which we can call. But we only get a free trial. Once we have identified genes relevant to asthma using this novel approach, and validated by pathway analysis, we can check against drug databases in an attempt to find potential drugs. We recognize that this part of the pipeline is still vague at this point. Details need to be work out in case we have the time and resource do so. 


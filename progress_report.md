## Overview

While major components of our analysis pipeline have remained the same from our initial propsal, we have made a small number of modifications, mostly surrounding the details of our methods. These ammendments have been informed by the progress of exploratory investigations as well as consultation with the course instructors. We have given a high level description of these changes and our current progress in this document. Links to individual .md files containing more details have been linked for each component.

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

Initial proposal: implement Hidden Covariates with Prior algorithm (Mostafavi, et al, written in matlab) to adjust for known (GC bias, age, sex) and unknown (comorbidities, batch effects) covariates.

When PCA is used to visualize the first two PCs of the data, we do not see any clusters based on known covariates. This means we don’t need to adjust for known covariates. We cannot use the HCP algorithm, because we do not have the unnormalized count table (we have the raw reads, however re-creating the count table from these reads would take extensive CPU time, which would detract from time allotted to other parts of the project). However, as can be seen in the results section, no visible clusters can be seen through PCA of our data. 


### Methylation sites preprocessing & mapping

The exact normalization methods the data deposited on GEO has undergone are not clear, so we are clarifying with the authors. 

The methylation data deposited on GEO contains probes that map to unique genes, so we do not need to remove promiscuous probes. If they have already used the SWAN algorithm to correct for probe bias, we will not do any further adjustment. PCA of all methylation data suggests that expression does not depend on known covariates. 

It is still possible for multiple probes to map to the same gene. In this case, we will obtain a summary value by taking, likely, the maximum value. Other options include taking the mean or the median. There is no clear consensus on the best approach in the current literature.  


## Aim 2: Patient clustering

This step has remainted largely unchanged since our initial proposal. 

We are using k-means, since we have a predefined number of clusters (chose 2). CLCA1 shows very clear high/level expression clusters, the other genes less so. We have 6 Th2-high and 51 Th2-low (though some in the low-group are better choices than others). 

This, of course, results in an imbalanced dataset. We will likely perform bootstrapping to resolve this issue. 


## Aim 3: Hypothesis testing

This aim is the key step in our analysis where we investigate and test our hypothesis formally. It consists of 2 separate steps: network construction and differential network analysis. Details below:

### Network construction

From aim 2: patient clustering, we have obtained two clusters of patients 



### Differential network analysis




## Aim 4: Biological interpretations


### Pathway analys


### Drug target identification (maybe)


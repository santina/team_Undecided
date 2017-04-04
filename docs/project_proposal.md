#### Using Differential Network Analysis (DiNA) to identify upstream regulators of CLCA1, periostin, and serpinB2 in Th2-low and high Asthma endotypes. 

Asthma, a disease characterized by chronic inflammation, affects over 235 million individuals worldwide (1). Broadly, asthma can be divided into two subgroups: patients that have high levels of T helper cells 2 (Th2) cytokines (Th2-high) and patients that have low levels of Th2 cytokines(Th2-low). Th2-high patients tend to show more severe symptoms, but respond better to therapies (2). Within airway epithelial cells (AECs), high expression of CLCA1, periostin, and serpinB2 genes in response to Th2 cytokines have been used as biomarkers to differentiate Th2-high and Th2-low populations. Exploring the biological differences between these two endotypes may illuminate upstream regulators that drive this immune response, which could become potential targets for drug therapies. 	

The integration of RNA- Seq and methylation data represents a potentially fruitful approach to obtaining a comprehensive understanding of the underlying biology of asthma. RNA-Seq provides a high-resolution snapshot of gene expression in a cell at any given time, while methylation data gives an indication of the epigenetic variation in asthma endotypes. In a previous study, methylation profiles from airway epithelial cells (AEC) were obtained from healthy and asthmatic patients (3). Co-expression network analysis was used to group differentially methylated CpGs (DMCs) into comethylation modules (3). Modules were highly correlated with phenotypic signatures, indicating DNA methylation plays a central role in promoting distinct molecular pathways of asthma pathogenesis, and therefore clinical heterogeneity (3). Here we propose to extend this observation by performing an integrated analysis of transcriptional and methylation profiles to identify key regulatory molecules that drive Th2 low and high endotypes.

Our goal is to use the RNA-Seq and methylation profiles of Th2-high, Th2-low, and control patients to generate co-expression networks with the R package, Weighted Co-expression Network Analysis (WGCNA) (4). We intend to define the Th2-high and Th2-low asthma subgroups by clustering on expression levels of CLCA2, periostin, and serapinB2. We then propose to apply DiNA, which is a recent class of network-based Bioinformatics algorithms that focuses on the differences in network topologies between two states (5). This approach builds on two lines of prior knowledge: differential expression analysis (DEA) and network expression analysis (NEA) (5). A combined approach allows one to examine changes between states in both single and groups of genes. We intend to use various R-implemented DiNA algorithms, including dna, mlDNA, and DINGO (6-8). The results of these methods can be compared for biological relevance, in order to reveal important regulators of asthma endotypes, as well as identify potential drug targets. 

These analyses will be conducted on methylation and gene expression data collected from samples of airway epithelial cells. AEC samples were obtained and isolated from 81 adult subjects by bronchoscopy (3). Of the 81 subjects, 54 were classified as asthmatic patients while 27 were non-asthmatics (3). Gene expression and methylation data was generated using Illumina HiSeq 2000 and Infinium Human Methylation 450K Bead Chip platforms, respectively (3). Ultimately, this study will provide a novel means of identified upstream regulators of CLCA1, Periostin, and SerpinB2, which can be potential drug targets. Furthermore, by probing the functions of upstream genes, we will learn more about the underlying etiology of asthma. 

##### References

1. Pawankar, R. (2014). Allergic diseases and asthma: a global public health concern and a call to action. World Allergy Organization Journal, 7(1), 12.
2. Wesolowska-Andersen, A., & Seibold, M. A. (2015). Airway molecular endotypes of asthma: dissecting the heterogeneity. Current opinion in allergy and clinical immunology, 15(2), 163.
3. Nicodemus-Johnson, J., Myers, R. A., Sakabe, N. J., Sobreira, D. R., Hogarth, D. K., Naureckas, E. T., ... & Nicolae, D. L. (2016). DNA methylation in lung cells is associated with asthma endotypes and genetic risk. JCI insight, 1(20).
4. Langfelder, P., & Horvath, S. (2008). WGCNA: an R package for weighted correlation network analysis. BMC bioinformatics, 9(1), 559.
5. Lichtblau, Y., Zimmermann, K., Haldemann, B., Lenze, D., Hummel, M., & Leser, U. (2016). Comparative assessment of differential network analysis methods. Briefings in bioinformatics, bbw061.
6. Gill, R., Datta, S., & Datta, S. (2010). A statistical framework for differential network analysis from microarray data. BMC bioinformatics, 11(1), 95.
7. Ma, C., Xin, M., Feldmann, K. A., & Wang, X. (2014). Machine Learning–Based Differential Network Analysis: A Study of Stress-Responsive Transcriptomes in Arabidopsis. The Plant Cell, 26(2), 520-537.
8. Ha, M. J., Baladandayuthapani, V., & Do, K. A. (2015). DINGO: differential network analysis in genomics. Bioinformatics, 31(21), 3413-3420.


##### Division of Labor 

Member | Background | Tasks | Contributions
 --- | --- | --- | ---
Emma Graham | Biochemistry and Biophysics | implement Hidden Covariates with Prior algorithm (Mostafavi, et al, written in matlab) to normalize RNAseq data and pre-process 450k methylation array using the SWAN algorithm in the minfi package | Alignment/preprocessing/normalization
Allison Tai | Computer science and Biochemistry | Cluster asthma patients and build separate networks integrating RNAseq and methylation data for each asthma subgroup and the controls using WGCNA | Cluster and build networks 
Eric Chu | Software engineering and Neuroscience | Identify preserved and non-preserved modules between each network and quantify the degree of module preservation through the calculation of module preservation statistics | Differential Network Analysis: identifying differentially expressed modules 
Arjun Baghela | Biochemistry | Perform pathway analysis on the modules identified in each network and check to see if existing drugs can target affect these genes | Differential Network Analysis: pathway analysis within preserved and non-preserved modules; Biological interpretation; Differential Network Analysis: pathway analysis

##### Detailed Aims/Methodologies 

1. Normalize RNAseq and methylation array data to account for batch effects and within-sample variation
	* RNAseq
		* align reads to  genome (STAR/TopHat)
		* implement Hidden Covariates with Prior algorithm (Mostafavi, et al, written in matlab) to adjust for known (GC bias, age, sex) and unknown (comorbidities, batch effects) covariates
	* Methylation probe array
		* Explore different normalization methods. Tentatively, 450k methylation array will be pre-processed using the SWAN algorithm in the minfi package. 

2. Cluster all asthma patients using expression profiles of CLCA1, periostin and serpinB2 to create two subgroups of patients: high Th2 and low Th2 
	* Subgroups will be defined by setting a boundary threshold (where above is high and below is low)

3. Build integrated expression/methylation co-expression networks for subgroups (WGCNA) 
	* Build separate gene expression and co-methylation networks for subgroup to try (WGCNA) 
	* Each network will have two different types of nodes: one type corresponding to CpG sites and other corresponding to the expression of a gene. 
	* We will also build a network within our control group to produce three networks in total
prune away uninformative nodes & edges

4. DiNA (Try many R packages)
	* Identify preserved and non-preserved modules within Th2-high, Th2-low and control patients
	* Quantify the degree of preservation through the calculation of module preservation statistics

5. Pathway Enrichment on Preserved Nodes/Modules. 

6. Look at Drug Databases. 
	* Not sure to what extent we will be able to do this. 
	* DrugBank.ca has an API which we can call. But we only get a free trial. 


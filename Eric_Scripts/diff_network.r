library(tidyverse)
library(dplyr)
# packages
library(dna)
library(limma)


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


# ----- main -----

setwd("/Users/ericchu/ws/team_Undecided/")

# raw data
rawMetadata <- read_csv("Raw_Data/GSE85566_metadata.txt")

rawCounts <- read.table("Raw_Data/GSE85567_RNASeq_normalizedcounts.txt", check.names = FALSE) %>% 
  rownames_to_column() %>% as_tibble() %>% 
  select(geneId = rowname, everything())

rawCountsT <- rawCounts %>% prepareTransposedCounts()
  


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
                              p.value = 0.2,
                              sort.by = "p")
  


# separate into 2 groups - control & asthma
asthmaSamplesMetadata <- rawMetadata %>% filter(Status == "Asthma")
controlSamplesMetadata <- rawMetadata %>% filter(Status == "Control")


# try running dna with differentially expressed genes

# find genes with differential connectivity! - individual genes!
asthmaExpression <- rawCountsT %>% filter(sample %in% asthmaSamplesMetadata$ID)
asthmaExpression <- asthmaExpression[rownames(diffGenes)] %>% as.data.frame()
controlExpression <- rawCountsT %>% filter(sample %in% controlSamplesMetadata$ID)
controlExpression <- controlExpression[rownames(diffGenes)] %>% as.data.frame()
geneLevelResult <- test.individual.genes(controlExpression, asthmaExpression, 
                                         scores = "cor", distance = "abs", 
                                         num.permutations = 1000)


# now try running it with random set of genes and see if we observe differences!

randomGenes <- rawCounts$geneId %>% sample(51)
asthmaExpressionRandom <- rawCountsT %>% filter(sample %in% asthmaSamplesMetadata$ID)
asthmaExpressionRandom <- asthmaExpressionRandom[randomGenes] %>% as.data.frame()
controlExpressionRandom <- rawCountsT %>% filter(sample %in% controlSamplesMetadata$ID)
controlExpressionRandom <- controlExpressionRandom[randomGenes] %>% as.data.frame()
geneLevelResultRandom <- test.individual.genes(controlExpressionRandom, asthmaExpressionRandom, 
                                         scores = "cor", distance = "abs", 
                                         num.permutations = 1000)



# ----- sandbox (for reference) -----

data("HeavyMice")
data("LeanMice")

args(test.individual.genes)

# comparing connectivity between individual genes
individualResult <- test.individual.genes(LeanMice, HeavyMice, 
                                scores = "cor", distance = "abs", 
                                rescale.scores = TRUE, num.permutations = 500)

# result %>% get.results() %>% rownames_to_column() %>% as_tibble() %>% filter(p.value == 0) %>% arrange(rowname) %>% View()

# comparing connectivity between gene sets

genes <- c("AA408451", "BC010552", "BC026585")
geneSetResult <- test.class.genes(LeanMice, HeavyMice, genelist = genes,
                                  scores = "cor", distance = "abs",
                                  rescale.scores = TRUE, num.permutations = 500)




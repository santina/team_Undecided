library(tidyverse)
library(dplyr)

library(dna)

# load the data
setwd("/Users/ericchu/ws/team_Undecided/")


# ----- functions -----







# ----- main -----






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




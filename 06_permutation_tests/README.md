# Permutation tests

We used phylogeny-wide permutation tests, similar to a bootstrap approach. Permutation tests determine the statistical significance of observed data by comparing it against an empirical null distribution generated from numerous random permutations of the original data (e.g., evolutionary rates of novel genes in a node, GOs in a node) with their labels (aquatic or terrestrial) reshuffled in each permutation. This allows us to test specific hypotheses about evolutionary patterns by breaking the observed correlation between variables while preserving their individual distributions.


In the first permutation test, we evaluated if the rate of emergence of novel genes (number of novel genes emerging per million year) in terrestrial nodes was significantly higher than aquatic nodes. Please find the R script in **permutation_novelHGcounts_and_timings_rates.R**, and the input file **aquatic_terrestrial_novelHGcounts_and_timings.csv**.

The second permutation test assessed if the biological functions found in terrestrial nodes are significantly different from those in other nodes. Please find the R script in **permutation_novel_GOterms_terrestrial_aquatic.R**, and the input file **novel_GOterms_terrestrial_aquatic_df.csv**.





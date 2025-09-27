# Read the TSV file
go_data <- read.table("merge_node11novel_bilateria_HGcount.tsv", header = TRUE, sep = "\t")

##structure of input tsv file example:
##TermID	Bilateria_Hsap_Number_of_Homology_Groups	node11expansion_Number_of_Homology_Groups
##GO:0008150	6309	24
##GO:0005575	6378	24
##GO:0005623	6323	24
##GO:0044464	6323	24
##GO:0003674	5432	24
##GO:0009987	5806	24
##GO:0005622	6158	22
##GO:0044424	6145	22
##GO:0043226	5693	20

# View the first few rows of the data to ensure it's loaded correctly
head(go_data)

total_bilateria <- 7195 # Total HGs in Hsap Bilateria; 
total_node11_novel <- 91 # Total HGs in Hsap node11_novel;

# store p-values from fisher's exact test for each GO term
p_values <- numeric(nrow(go_data))

# loop through each GO term
for(i in 1:nrow(go_data)) {
  a <- go_data$Bilateria_Hsap_Number_of_Homology_Groups[i]
  b <- go_data$node11novel_Number_of_Homology_Groups[i]
  c <- total_bilateria - a
  d <- total_node11_novel - b
  
#contingency table for the current GO term
  contingency_table <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE,
                              dimnames = list(c("With GO Term", "Without GO Term"),
                                              c("Bilateria", "Node11_novel")))
  
# perform fisher's exact test and store the p-value
  test_result <- fisher.test(contingency_table)
  p_values[i] <- test_result$p.value
}

# adjust p-values for multiple comparisons using Benjamini-Hochberg method
p_values_adjusted <- p.adjust(p_values, method = "BH")


results <- data.frame(TermID = go_data$TermID, Adjusted_P_Value = p_values_adjusted)


print(results)

# calculate enrichment ratios
enrichment_ratios <- numeric(nrow(go_data))
for(i in 1:nrow(go_data)) {
  a_ratio <- go_data$Bilateria_Hsap_Number_of_Homology_Groups[i] / total_bilateria
  b_ratio <- go_data$node11novel_Number_of_Homology_Groups[i] / total_node11_novel
  enrichment_ratios[i] <- b_ratio / a_ratio
}

# add enrichment ratios to the results
results$Enrichment_Ratio <- enrichment_ratios

write.csv(results, "results_p_ratio_all.csv", row.names = FALSE)

# Filter for significant enrichment in Node11_novel
# Filter both adjusted p-value threshold and enrichment ratio
enriched_in_node11_novel <- results[results$Adjusted_P_Value <= 0.05 & results$Enrichment_Ratio > 1, ]

print(enriched_in_node11_novel)

# Write the results to csv file
write.csv(enriched_in_node11_novel, "results_enrichedGO_node11novel_bilateria_back.csv", row.names = FALSE)


rm( list = ls( ) )

library(vegan)
library(ggplot2)

# Load data
setwd("path_to_your_input_file/")

df <- read.csv("novel_pfams_terrestrial_species_df.csv",
               stringsAsFactors = FALSE, check.names = FALSE)

species_names <- df[[1]]
groups <- as.factor(df[[ncol(df)]]) 

# Define feature columns
feature_cols <- 2:(ncol(df) - 1)
feature_cols
# Extract and coerce to numeric matrix
mat <- as.matrix(df[, feature_cols])
mat
# Coerce safely in case some columns came in as factors/characters
mat <- apply(mat, 2, function(x) as.numeric(as.character(x)))
rownames(mat) <- species_names
mat
# Clean and remove zero variance columns
vars <- apply(mat, 2, var, na.rm = TRUE)
if (any(vars == 0)) {
  message("Dropping ", sum(vars == 0), " zero-variance columns.")
  mat <- mat[, vars > 0, drop = FALSE]
}
mat

# PCoA with Jaccard distance 
# Calcualte pairwise Jaccard dissimilarities between all species
dist_jaccard <- vegdist(mat, method = "jaccard", binary = TRUE)
dist_jaccard

pcoa_jac <- cmdscale(dist_jaccard, eig = TRUE, k = 2)
pcoa_jac
eigvals <- pcoa_jac$eig
eigvals
prop_var <- round(100 * (eigvals / sum(eigvals))[1:2], 1)

# Build plotting df including groups
pcoa_df <- data.frame(
  PCo1 = pcoa_jac$points[, 1],
  PCo2 = pcoa_jac$points[, 2],
  groups = groups,
  label = rownames(pcoa_jac$points)
)
pcoa_df

# Combining groups for ellipses
groupA_members <- c("Group1_bdelloidea", "Group2_clitellata",
                    "Group4_nematoda", "Group5_tardigrada",
                    "Group6_onychophora")
groupA_members <- intersect(groupA_members, as.character(pcoa_df$groups))
pcoa_df$combined_groups <- ifelse(as.character(pcoa_df$groups) %in% groupA_members,
                                  "Group_A", "Group_B")
pcoa_df$combined_groups <- factor(pcoa_df$combined_groups, levels = c("Group_A", "Group_B"))
pcoa_df$combined_groups 

grp <- pcoa_df$combined_groups
feature_names <- colnames(mat)
feature_names

res_list <- lapply(feature_names, function(f) {
  tab <- table(mat[, f], grp)
  if (nrow(tab) < 2 || ncol(tab) < 2) return(NULL)  # skip features with no variation
  ft <- fisher.test(tab)
  propA <- sum(mat[grp == "Group_A", f]) / sum(grp == "Group_A")
  propB <- sum(mat[grp == "Group_B", f]) / sum(grp == "Group_B")
  data.frame(
    feature = f,
    p_value = ft$p.value,
    odds_ratio = if (!is.null(ft$estimate)) as.numeric(ft$estimate) else NA,
    prop_Group_A = propA,
    prop_Group_B = propB,
    abs_diff = abs(propA - propB),
    stringsAsFactors = FALSE
  )
})
res_list

fisher_df <- do.call(rbind, res_list)

# Multiple testing correction
fisher_df$padj <- p.adjust(fisher_df$p_value, method = "BH")

# Filter significant and substantial differences
# Label direction of enrichment
fisher_df$enriched_in <- ifelse(fisher_df$prop_Group_A > fisher_df$prop_Group_B,
                                "Group_A", "Group_B")
# Filter and test significance
sig <- subset(fisher_df, padj < 0.01 & abs_diff > 0.1)
sig <- sig[order(-sig$abs_diff), ]

# Show top 10 with key columns
head(sig[, c("feature", "enriched_in", "prop_Group_A", "prop_Group_B", "abs_diff", "padj", "odds_ratio")], 10)

# Save results
write.csv(fisher_df, "pfams_fisher_all_features.csv", row.names = FALSE)
write.csv(sig, "pfams_fisher_significant_enriched.csv", row.names = FALSE)

rm( list = ls( ) )

library(car)

setwd("novel_goterms_pfams/")
novel_goterms_df <- read.csv("novel_goterms_terrestrialspecies_df.csv")

species_names <- novel_goterms_df[, 1]

## there are 12117 columns
groups <- novel_goterms_df[, 12117]

##center the data
novelgoterms_df_centered <- scale(novel_goterms_df[, 2:12116], scale = FALSE, center = TRUE) 

##scale the data
novelgoterms_df_centered_scaled <- scale(novelgoterms_df_centered,  center = FALSE, scale = TRUE)

apply(novelgoterms_df_centered_scaled,2,sd)

pca_novelgoterms <- prcomp(novelgoterms_df_centered_scaled)

summary(pca_novelgoterms)

###write loadings to a csv file
loadings <- pca_novelgoterms$rotation
write.csv(loadings, "loadings.csv", row.names = TRUE)
###write scores to a csv file
pca_novelgoterms$x
scores <- pca_novelgoterms$x
# set species names as row names in the rotation matrix
if(length(species_names) == nrow(scores)) {
  rownames(scores) <- species_names
} else {
  warning("The number of species names does not match the number of rows in the rotation matrix.")
}
head(rownames(scores))
write.csv(scores, "scores.csv", row.names = TRUE)


plot(pca_novelgoterms, main = "")
barplot(pca_novelgoterms$rotation[,1], main = "")

pca_scores <- as.data.frame(pca_novelgoterms$x)


pca_scores$species = species_names
pca_scores$groups = groups

library(ggplot2)

my_colors <- c("Group1_bdelloidea"="#a6cee3", "Group2_clitellata"="#1f78b4", 
               "Group3_stylommatophora"="#b2df8a", "Group4_nematoda"="#33a02c", 
               "Group5_tardigrada"="#000000", "Group6_onychophora"="#e31a10", 
               "Group7_arachnida"="#fdbf6f", "Group8_myriapoda"="#b15928", 
               "Group9_armadillidium"="#f09cce", "Group10_hexapoda"="#9b59e3", 
               "Group11_tetrapoda"="#0ad1ca")

ggplot(pca_scores, aes(x = PC1, y = PC2, color = as.factor(groups))) +
  geom_point(alpha = 1, size = 2) +
  scale_color_manual(values = my_colors) +
  labs(title = "PCA of Novel GO Terms",
       x = "Principal Component 1 (PC1)",
       y = "Principal Component 2 (PC2)") +
  theme_minimal()


###statistical analysis for pc1
anova_pc1 <- aov(PC1 ~ groups, data = pca_scores)
summary(anova_pc1)
# perform Tukey's HSD test
tukey_results <- TukeyHSD(anova_pc1)
tukey_results

###statistical analysis for pc2
anova_pc2 <- aov(PC2 ~ groups, data = pca_scores)
summary(anova_pc2)
# Perform Tukey's HSD test
tukey_results <- TukeyHSD(anova_pc2)
tukey_results

# Combine specific groups into two broader groups for ellipses
pca_scores$combined_groups <- ifelse(pca_scores$groups %in% c("Group1_bdelloidea", "Group2_clitellata", 
                                                              "Group4_nematoda", "Group5_tardigrada", 
                                                              "Group6_onychophora"), "Group_A", "Group_B")

# Load the ggforce package
library(ggforce)
ellipse_colors <- c("Group_A" = "orange", "Group_B" = "darkgreen")

# Plot PCA with ellipses for combined groups
ggplot(pca_scores, aes(x = PC1, y = PC2, color = as.factor(groups))) +
  geom_point(alpha = 1, size = 3) +
  scale_color_manual(values = my_colors) +
  labs(title = "PCA of Novel GO Terms",
       x = "Principal Component 1 (PC1)",
       y = "Principal Component 2 (PC2)") +
  theme_minimal() +
  stat_ellipse(aes(group = combined_groups, color = combined_groups), type = "norm", linetype = 1, size = 1.5) +
  scale_color_manual(values = c(my_colors, ellipse_colors))

# perform MANOVA
manova_result <- manova(cbind(PC1, PC2) ~ combined_groups, data = pca_scores)
summary(manova_result)

# further look into the results
summary.aov(manova_result)

# perform pairwise comparison with Tukey HSD test
anova_pc1 <- aov(PC1 ~ combined_groups, data = pca_scores)
anova_pc2 <- aov(PC2 ~ combined_groups, data = pca_scores)

# perform Tukey's HSD test on both principal components
tukey_pc1 <- TukeyHSD(anova_pc1)
tukey_pc2 <- TukeyHSD(anova_pc2)

tukey_pc1
tukey_pc2



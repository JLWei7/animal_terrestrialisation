rm( list = ls( ) )

library(vegan)
library(ggplot2)

# ---- Load data ----
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
# # Coerce safely in case some columns came in as factors/characters
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



# Plot
my_colors <- c("Group1_bdelloidea"="#a6cee3", "Group2_clitellata"="#1f78b4", 
               "Group3_stylommatophora"="#b2df8a", "Group4_nematoda"="#33a02c", 
               "Group5_tardigrada"="#000000", "Group6_onychophora"="#e31a10", 
               "Group7_arachnida"="#fdbf6f", "Group8_myriapoda"="#b15928", 
               "Group9_armadillidium"="#f09cce", "Group10_hexapoda"="#9b59e3", 
               "Group11_tetrapoda"="#0ad1ca")

# Plot without labels
ggplot(pcoa_df, aes(x = PCo1, y = PCo2, color = as.factor(groups))) +
  geom_point(alpha = 0.9, size = 3) +
  scale_color_manual(values = my_colors) +
  labs(title = "PCA of Novel PFAMS",
       x = "Principal Component 1 (PC1)",
       y = "Principal Component 2 (PC2)") +
  theme_minimal()

# Plot with labels
ggplot(pcoa_df, aes(x = PCo1, y = PCo2, color = as.factor(groups))) +
  geom_point(size = 3, alpha = 0.9) +
  scale_color_manual(values = my_colors) +
  geom_text(aes(label = label), check_overlap = TRUE, size = 3, vjust = -0.7) +
  labs(
    title = "PCoA of PFAMS on Jaccard Distance",
    x = paste0("PCo1 (", prop_var[1], "%)"),
    y = paste0("PCo2 (", prop_var[2], "%)"),
    color = "Group"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "right")


# Combining groups for ellipses
groupA_members <- c("Group1_bdelloidea", "Group2_clitellata",
                    "Group4_nematoda", "Group5_tardigrada",
                    "Group6_onychophora")

groupA_members <- intersect(groupA_members, as.character(pcoa_df$groups))
pcoa_df$combined_groups <- ifelse(as.character(pcoa_df$groups) %in% groupA_members,
                                  "Group_A", "Group_B")
pcoa_df$combined_groups <- factor(pcoa_df$combined_groups, levels = c("Group_A", "Group_B"))
pcoa_df$combined_groups 

# Set Color
orig_groups <- levels(pcoa_df$groups)
base_colors <- setNames(scales::hue_pal()(length(orig_groups)), orig_groups)
ellipse_colors <- c("Group_A" = "orange", "Group_B" = "darkgreen")
color_mapping <- c(my_colors, ellipse_colors)


# Change the semi-terrstrial speceis shape to square
pcoa_df$shape_group <- ifelse(as.character(pcoa_df$groups) %in% groupA_members,
                              "Group_A", "Group_B")

shape_mapping <- c("Group_A" = 15,  # square
                   "Group_B" = 16)  # circle

# Plot with ellipses
p <- ggplot(pcoa_df, aes(x = PCo1, y = PCo2)) +
  geom_point(aes(color = groups, shape = shape_group), size = 3, alpha = 0.9) +
  stat_ellipse(aes(group = combined_groups, color = combined_groups),
               type = "norm", linetype = 1, size = 1.2, level = 0.95, show.legend = FALSE) +
  scale_color_manual(
    values = color_mapping,
    name = "Group"
  ) +
  scale_shape_manual(
    values = shape_mapping,
    guide = "none"
  ) +
  labs(
    title = "PCoA of Pfams in Novel Genes",
    x = paste0("PCo1 (", prop_var[1], "%)"),
    y = paste0("PCo2 (", prop_var[2], "%)")
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "right")

print(p)


# Variance
var_df <- data.frame(Axis = c("PCo1", "PCo2"), Proportion = prop_var)
var_df


library(vegan)

# Check dispersion homogeneity
bd <- betadisper(dist_jaccard, pcoa_df$combined_groups)
print(permutest(bd))

# PERMANOVA to check overall difference
set.seed(42)
perm_res <- adonis2(dist_jaccard ~ combined_groups, data = pcoa_df,
                    permutations = 10000, method = "jaccard")
print(perm_res)


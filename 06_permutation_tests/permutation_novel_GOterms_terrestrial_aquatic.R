rm( list = ls( ) )

library(car)
library(vegan) 

setwd("Path_to_the_csv_file")
novel_goterms_df <- read.csv("novel_GOterms_terrestrial_aquatic_df.csv")
str(novel_goterms_df)
colnames(novel_goterms_df)[1:10]
tail(colnames(novel_goterms_df))

# Prepare data
goterm_matrix <- as.matrix(novel_goterms_df[, 2:(ncol(novel_goterms_df)-1)])
# Get group labels 
group_labels <- as.factor(novel_goterms_df[[ncol(novel_goterms_df)]])
# Remove rows that have no annotation
non_empty <- rowSums(goterm_matrix) > 0   
cat(sum(!non_empty), "rows had zero GO terms and were removed\n")
goterm_matrix <- goterm_matrix[non_empty, ]
group_labels  <- group_labels [non_empty]

# Jaccard distance on group‐level presence/absence
# Compute observed group‐level binary profiles
group1_profile <- colMeans(goterm_matrix[group_labels=="aquatic", , drop=FALSE]) > 0
group2_profile <- colMeans(goterm_matrix[group_labels=="terrestrial", , drop=FALSE]) > 0

observed_jaccard <- vegdist(
  rbind(as.numeric(group1_profile), as.numeric(group2_profile)),
  method = "jaccard"
)[1]

# Permutation
n_perm <- 10000
perm_jaccard <- numeric(n_perm)
set.seed(42)
for(i in seq_len(n_perm)) {
  perm_labels <- sample(group_labels)
  p1 <- colMeans(goterm_matrix[perm_labels=="aquatic", , drop=FALSE]) > 0
  p2 <- colMeans(goterm_matrix[perm_labels=="terrestrial", , drop=FALSE]) > 0
  perm_jaccard[i] <- vegdist(rbind(as.numeric(p1), as.numeric(p2)),
                             method="jaccard")[1]
}

p_jaccard <- mean(perm_jaccard >= observed_jaccard)

cat(sprintf("Observed Jaccard = %.4f\n", observed_jaccard))
cat(sprintf("Permutation p-value (Jaccard) = %.10f\n", p_jaccard))

hist(perm_jaccard, breaks=100,
     main="Null distribution (Jaccard)", xlab="Jaccard Distance")
abline(v = observed_jaccard, col="red", lwd=2)
library(ggplot2)

# Put permutation values into a data frame
df_perm <- data.frame(Jaccard = perm_jaccard)

# Set binwidth
bw_fd <- 2 * IQR(df_perm$Jaccard) / length(df_perm$Jaccard)^(1/3)

# Plot
ggplot(df_perm, aes(x = Jaccard)) +
  geom_histogram(binwidth = bw_fd,
                 fill      = "#3182bd",   # blue fill
                 colour    = "white",
                 alpha     = 0.9) +
  geom_density(colour = "#08519c", size = 1, adjust = 1) +
  geom_vline(xintercept = observed_jaccard,
             colour = "red", size = 1.2) +
  geom_area(data = subset(df_perm, Jaccard >= observed_jaccard),
            aes(y = ..density..), stat = "bin",
            binwidth = bw_fd, fill = "red", alpha = 0.25) +
  annotate("text",
           x = observed_jaccard,
           y = Inf, vjust = 1.2, hjust = 1.05,
           colour = "red", size = 4.5,
           label = sprintf("Observed = %.3f\np = %.3g", 
                           observed_jaccard, p_jaccard)) +
  labs(title = "Permutation null distribution of Jaccard distance",
       x = "Jaccard distance",
       y = "Frequency") +
  theme_classic(base_size = 14) +
  theme(
    plot.title   = element_text(face = "bold"),
    panel.grid.major.y = element_line(colour = "grey90"),
    panel.grid.minor.y = element_blank()
  )



rm( list = ls( ) )

set.seed(2025)
# Read data
df <- read.csv("aquatic_terrestrial_novelHGcounts_and_timings.csv", comment.char = "#")

terr_ids <- df$node_id[df$habit_type == "terrestrial"]
terr_ids

# Set bin
df$exp_bin <- cut(df$exposure, breaks = 4)
df$exp_bin

# Set aquatic pool
bins_needed   <- unique(df$exp_bin[df$node_id %in% terr_ids])
aquatic_pool  <- df$node_id[df$habit_type == "aquatic" &
                              df$exp_bin %in% bins_needed]
aquatic_pool

# Define an exposure-scaled statistic
statistic <- function(ids){
  with(df[df$node_id %in% ids, ],
       sum(novel_HG_count) / sum(exposure))      # total gains per Myr
}

# Calculate observed profile of terrestrial nodes
S_obs <- statistic(terr_ids)
S_obs

# Permutation
B       <- 10000
S_boot  <- numeric(B)
for(b in seq_len(B)){
  # sample with replacement
  draw      <- sample(aquatic_pool, length(terr_ids), replace = TRUE)
  S_boot[b] <- statistic(draw)
}

# Empirical p-value and SES
p_boot <- (sum(S_boot >= S_obs) + 1) / (B + 1)
p_boot
SES_boot <- (S_obs - mean(S_boot)) / sd(S_boot)
SES_boot
cat(sprintf("Bootstrap p = %.5f; SES = %.2f\n", p_boot, SES_boot))


### Plot
library(ggplot2)

# Prepare the permutation results
df_perm <- data.frame(rate = S_boot)

# Set the number of bins
bins <- 100
max_rate <- ceiling(max(df_perm$rate) * 2) / 2
# Build the plot
ggplot(df_perm, aes(x = rate)) +
  ## histogram of all null draws
  geom_histogram(bins      = bins,
                 fill      = "#3182bd",
                 colour    = "white",
                 alpha     = 0.9) +
  ## vertical line at observed terrestrial rate
  geom_vline(xintercept = S_obs,
             colour      = "red",
             size        = 1.2) +
  ## annotation text in the top‐right corner
  annotate("text",
           x     = S_obs,
           y     = Inf, vjust = 1.1, hjust = 1.1,
           colour = "red", size = 4.5,
           label = sprintf("Observed = %.3f\np = %.3g",
                           S_obs, p_boot)) +
  scale_x_continuous(
    breaks = seq(0, max_rate, by = 0.5),
    limits = c(0, max_rate),
    expand = c(0, 0)
  ) +
  ## axis labels and title
  labs(title = "Distribution of Novel HG gain rate",
       x     = "Novel HG gains per Myr (random aquatic sets)",
       y     = "Count") +
  ## styling
  theme_classic(base_size = 14) +
  theme(
    plot.title         = element_text(face = "bold"),
    panel.grid.major.y = element_line(colour = "grey90"),
    panel.grid.minor.y = element_blank()
  )


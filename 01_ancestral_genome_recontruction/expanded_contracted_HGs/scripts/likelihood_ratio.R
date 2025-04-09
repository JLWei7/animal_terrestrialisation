# Install and load the lmtest package

library(lmtest)

# Set the final likelihoods for your models
lnL_split1_model1_twolambda <- -518263.54009812
lnL_split1_model2_threelambda <- -512685.1994883

lnL_split2_model1_twolambda <- -373630.16908639
lnL_split2_model2_threelambda <- -369121.97760656

lnL_split3_model1_twolambda <- -388264.08624134
lnL_split3_model2_threelambda <- -387866.17943892

# Calculate the likelihood ratio
lr_split1 <- -2 * (lnL_split1_model1_twolambda  - lnL_split1_model2_threelambda)
lr_split2 <- -2 * (lnL_split2_model1_twolambda  - lnL_split2_model2_threelambda)
lr_split3 <- -2 * (lnL_split3_model1_twolambda  - lnL_split3_model2_threelambda)

print(lr_split1)
print(lr_split2)
print(lr_split3)

# Calculate the degrees of freedom difference between the two models
df_difference <- 1  # Adjust this based on your models

# Perform a chi-squared test
p_value1 <- pchisq(lr_split1, df = df_difference, lower.tail = FALSE)
p_value2 <- pchisq(lr_split2, df = df_difference, lower.tail = FALSE)
p_value3 <- pchisq(lr_split3, df = df_difference, lower.tail = FALSE)

# Print the p-value with more decimal places
print(p_value1, digits = 22)
print(p_value2, digits = 22)
print(p_value3, digits = 22)

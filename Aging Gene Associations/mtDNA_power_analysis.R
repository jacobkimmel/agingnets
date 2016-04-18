# mtDNA change power analyses

library(pwr)

# Determine effect size

# Cohen's d: 
# mean(sample1) - mean(sample2) / sqrt( (stdev(sample1)**2 + stdev(sample2)**2) / 2 )

effect_size <- function(sample1, sample2){
  mean1 <- mean(sample1)
  mean2 <- mean(sample2)
  sd1 <- sd(sample1)
  sd2 <- sd(sample2)
  n1 <- length(sample1)
  n2 <- length(sample2)
  s <- sqrt( ((n1-1)*sd1**2 + (n2-1)*sd2**2) / (n1 + n2 - 2) )
  d <- (mean1-mean2) / s
  eff_return <- c(d, s)
  return(eff_return)
}

# Power level of pi = 0.8 is generally considered the baseline for adequacy 
# Implies a 4:1 ratio of Type II error (beta)-risk to Type I error (alpha)-risk
pwr_result <- pwr.t.test(d = d, power = 0.8, sig.level = 0.05)

# Heart mtDNA power analysis 

heart_WT_3 <- c(4077.892,
                2990.118,
                4841.29)

heart_KO_3 <- c(3998.493,
                2950.668,
                3903.672)


heart_WT_18 <- c(3901.340982,
              3539.554995,
              4829.158896,
              2756.447153)

heart_KO_18 <- c(2904.853352,
                 2034.729114,
                 2945.893062,
                 3336.36092,
                 2205.502637)

heart_WT_30 <- c(3346.330778,
                 4167.915005,
                 4832.042539,
                 4342.078213)

heart_KO_30 <- c(1916.645932,
                 2719.624265,
                 5362.861879,
                 3670.923019)

d_s_3 <- effect_size(heart_WT_3, heart_KO_3)
d_3 <- d_s_3[1]
d_s_18 <- effect_size(heart_WT_18, heart_KO_18)
d_18 <- d_s_18[1]
d_s_30 <- effect_size(heart_WT_30, heart_KO_30)
d_30 <- d_s_30[1]

heart3_pwr_analysis <- pwr.t.test(d = d_3, power = 0.8, sig.level = 0.05)
heart18_pwr_analysis <- pwr.t.test(d = d_18, power = 0.8, sig.level = 0.05)
heart30_pwr_analysis <- pwr.t.test(d = d_30, power = 0.8, sig.level = 0.05)

# Liver power analysis 

liver_WT_3 <- c(1296.9850, 1551.5900)
liver_KO_3 <- c(1026.3430, 1810.6270, 1523.9340)

liver_WT_18 <- c(773.9223,
                 538.4603,
                 1161.7410,
                 811.8618)
liver_KO_18 <- c(762.3389,
                 737.1492,
                 1090.5250,
                 283.6961,
                 850.7031)

liver_WT_30 <- c(686.0900,
                 1070.2380,
                 753.4509,
                 876.4587)
liver_KO_30 <- c(849.6493,
                 1033.9060,
                 738.9418,
                 859.7640)

d_s_3 <- effect_size(liver_WT_3, liver_KO_3)
d_3 <- d_s_3[1]
d_s_18 <- effect_size(liver_WT_18, liver_KO_18)
d_18 <- d_s_18[1]
d_s_30 <- effect_size(liver_WT_30, liver_KO_30)
d_30 <- d_s_30[1]

liver3_pwr_analysis <- pwr.t.test(d = d_3, power = 0.8, sig.level = 0.05)
liver18_pwr_analysis <- pwr.t.test(d = d_18, power = 0.8, sig.level = 0.05)
liver30_pwr_analysis <- pwr.t.test(d = d_30, power = 0.8, sig.level = 0.05)


# GAS muscle power analysis 

gas_WT_3 <- c(2176.308,
              1651.309,
              1684.544)
gas_KO_3 <- c(2792.995,
              2985.773,
              2985.773)

gas_WT_18 <- c(3695.635,
               6223.212,
               1482.934,
               2263.116)
gas_KO_18 <- c(1854.779,
               1878.490,
               1111.404,
               3242.842)

gas_WT_30 <- c(2985.773,
               3030.753,
               2205.331,
               2755.691)
gas_KO_30 <- c(3053.732,
               2596.875,
               2622.862,
               2226.798)

d_s_3 <- effect_size(gas_WT_3, gas_KO_3)
d_3 <- d_s_3[1]
d_s_18 <- effect_size(gas_WT_18, gas_KO_18)
d_18 <- d_s_18[1]
d_s_30 <- effect_size(gas_WT_30, gas_KO_30)
d_30 <- d_s_30[1]

gas3_pwr_analysis <- pwr.t.test(d = d_3, power = 0.8, sig.level = 0.05)
gas18_pwr_analysis <- pwr.t.test(d = d_18, power = 0.8, sig.level = 0.05)
gas30_pwr_analysis <- pwr.t.test(d = d_30, power = 0.8, sig.level = 0.05)

# Brain power analysis 

brain_WT_3 <- c(447.4635,
                315.6744)
brain_KO_3 <- c(287.5053,
                229.9257,
                325.4717)
brain_WT_18 <- c(3901.341,
                 3539.555,
                 4829.159,
                 2756.447)
brain_KO_18 <- c(2904.853,
                 2034.729,
                 2945.893,
                 3336.361,
                 2205.503)

d_s_3 <- effect_size(brain_WT_3, brain_KO_3)
d_3 <- d_s_3[1]
d_s_18 <- effect_size(brain_WT_18, brain_KO_18)
d_18 <- d_s_18[1]

brain3_pwr_analysis <- pwr.t.test(d = d_3, power = 0.8, sig.level = 0.05)
brain18_pwr_analysis <- pwr.t.test(d = d_18, power = 0.8, sig.level = 0.05)
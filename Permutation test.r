#### Permutation test  ####

# Load necessary libraries
library(ggplot2)
library(dplyr)
library(ropls)

# Conducting 200 permutations, each time fitting a new OPLS-DA model
# Initialize vectors for R2, Q2, and correlation coefficients
R2_perm <- sacurine.plsda@suppLs$permMN[,2]
Q2_perm <- sacurine.plsda@suppLs$permMN[,3]
Corr_Coef <- sacurine.plsda@suppLs$permMN[,7]

# Uncomment the following code block to perform permutation tests
# for(i in 1:P_M) {
#   Y_perm <- sample(GroupFc)
#   Corr_Coef[i] <- mean(Y_perm == GroupFc)
#   opls_model_perm <- opls(dataMatrix, Y_perm, predI = 1, orthoI = 1, log10L = TRUE)
#   R2_perm[i] <- opls_model_perm@summaryDF$`R2X(cum)`
#   Q2_perm[i] <- opls_model_perm@summaryDF$`Q2(cum)`
# }

# Calculate R2 and Q2 for the original model
R2_orig <- sacurine.plsda@suppLs$permMN[1,2]
Q2_orig <- sacurine.plsda@suppLs$permMN[1,3]

# Prepare data for plotting
dat_opls <- as.data.frame(rbind(cbind(Corr_Coef, R2_perm, rep('R', length(R2_perm))),
                                  cbind(Corr_Coef, Q2_perm, rep('Q', length(Q2_perm)))))
dat_opls <- dat_opls %>%
  mutate(Corr_Coef = as.numeric(Corr_Coef), RQ = as.numeric(ifelse(V3 == 'R', R2_perm, Q2_perm)),
         Group = as.factor(V3))

# Create a scatter plot with regression lines
p <- ggplot(dat_opls, aes(x = Corr_Coef, y = RQ, color = Group)) +
    geom_point() +
    scale_color_manual(values = c("R" = "#1b9e77", "Q" = "#d95f02")) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.25)) +
    scale_y_continuous(breaks = seq(-1.5, 1, by = 0.5)) +
    labs(title = "OPLS-DA", x = "Correlation Coefficient", y = "R2Y and Q2") +
    geom_hline(yintercept = 0, size = 0.5, color = "black") +
    geom_segment(aes(x = 0, y = fitR$coefficients[1], xend = 1, yend = R2_orig),
                 colour = 'black', linetype = 'dotdash') +
    geom_segment(aes(x = 0, y = fitQ$coefficients[1], xend = 1, yend = Q2_orig),
                 colour = 'black', linetype = 'dotdash')

# Save the scatter plot
filename_OPLS <- paste0('Results_', M_version, '/', label_type, '/OPLS_', label_pos, '.pdf')
ggsave(filename_OPLS, plot = p, width = 5, height = 4)

# Create and save a histogram for permutation distribution
p_histogram <- ggplot(dat_opls, aes(x = RQ, fill = Group)) +
    geom_histogram(aes(y = after_stat(density)), color = 'white', alpha = 1, position = 'identity', binwidth = 0.03) +
    scale_fill_manual(values = c("#1b9e77", "#d95f02")) +
    labs(title = "Permutation Distribution for R2 and Q2", x = "R2 or Q2") +
    scale_x_continuous(breaks = seq(-1, 1, by = 0.25)) +
    geom_vline(xintercept = R2_orig, size = 0.5, color = "#d95f02") +
    geom_vline(xintercept = Q2_orig, size = 0.5, color = "#1b9e77") +
    theme_classic() +
    theme(legend.position = c(0.1, 0.85))

# Save the histogram
filename_OPLS2 <- paste0('Results_', M_version, '/', label_type, '/OPLS_hist_', label_pos, '.pdf')
ggsave(filename_OPLS2, plot = p_histogram, width = 5, height = 4)


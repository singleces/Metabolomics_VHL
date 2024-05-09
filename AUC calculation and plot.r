#### AUC calculation ####

# Load necessary libraries
library(glmnet)
library(pROC)
library(ggplot2)
library(dplyr)

# Initialize parameters
n_CV = 5  # Number of cross-validation folds
AUC_I <- matrix(NA, ncol = length(names_S), nrow = n_CV)  
colnames(AUC_I) <- names_S
AUC_M <- rep(NA, length(names_S)) 
AUC.coef <- AUC_M  
N1 = ncol(data_VHL) 
N2 = ncol(data_CON)  
n_fold1 = floor(N1 / n_CV)  
n_fold2 = floor(N2 / n_CV)  
Y <- c(rep(1, N1), rep(0, N2))  
dat_plot_mean <- data.frame()  

# Cross-validation loop
for(pi in 1:length(names_S)){
    p_name <- names_S[pi]
    roc_fold <- list()  # List to store ROC performance objects for each fold
    for(cvi in 1:n_CV){
        # Define indices for the test set for each fold
        te1 <- ((cvi - 1) * n_fold1 + 1):(cvi * n_fold1)
        te2 <- ((cvi - 1) * n_fold2 + 1):(cvi * n_fold2)
        # Define training and test data
        X_vec <- data_C_V_S[-c(te1, N1 + te2), p_name]
        Y_vec <- Y[-c(te1, N1 + te2)]
        X_vec.pre <- data_C_V_S[c(te1, N1 + te2), p_name]
        Y_vec.pre <- Y[c(te1, N1 + te2)]
        
        # Fit logistic regression model using glmnet
        fit <- glmnet(cbind(1, X_vec), Y_vec, family='binomial', lambda = 0)
        # Predict response and create ROC curves
        Y_fit <- as.vector(as.numeric(predict(fit, newx = cbind(1, X_vec.pre), type = "response")))
        pred <- prediction(Y_fit, Y_vec.pre)
        perf <- performance(pred, "sens", "fpr")
        auc <- performance(pred, 'auc')
        roc_fold[[cvi]] <- perf
        AUC_I[cvi, pi] <- unlist(slot(auc, "y.values"))
    }
    AUC_M[pi] <- mean(AUC_I[, pi])
    X <- data_C_V_S[, p_name]
    fit <- glmnet(cbind(1, X), Y, family='binomial', lambda = 0)
    AUC.coef[pi] <- coef(fit)[-1][2]

    # Process combined analysis for selected predictors after running both plasma and tissue analyses
    filename_both <- paste0("Results_", M_version, "/Combine/AUCgreaterThan", 100 * AUC_low, "_", label_pos, "_both.csv")
    dat_both <- read.csv(filename_both)

    if (p_name %in% dat_both$MS2.name) {
        X_new <- seq(0, 1, 0.001)
        YY <- matrix(NA, nrow = 5, ncol = length(X_new))
        for (i in 1:5) {
            U <- unlist(roc_fold[[i]]@x.values)
            V <- unlist(roc_fold[[i]]@y.values)
            for (j in 1:length(X_new)) {
                k = 1
                while (U[k] <= X_new[j] && k <= length(U)) {
                    YY[i, j] <- V[k]
                    k = k + 1
                }
            }
        }
        meanY = colMeans(YY)
        dat_temp <- data.frame(rep(p_name, each = length(X_new) + 1), c(0, X_new), c(0, meanY))
        colnames(dat_temp) <- c('MS2.name', 'Specificity', 'Sensitivity')
        dat_plot_mean <- rbind(dat_plot_mean, dat_temp)
    }
}

# Save the AUC data
PP_AUC <- data.frame(names_S, AUC_M, AUC.coef)
colnames(PP_AUC) <- c('MS2.name', 'AUC', 'AUC.coef')
PP_AUC <- inner_join(PP_AUC, Wilcox_result, by = "MS2.name")
PP_AUC <- PP_AUC %>% filter(AUC > AUC_low) %>% arrange(desc(AUC))
filename = paste0('Results_', M_version, '/', label_type, '/AUC', 100 * AUC_low, '_', label_pos, '.csv')
write.csv(PP_AUC, file = filename)





#### AUC lollipop plot ####
# Loop to handle both NEG and POS label positions
for (posi in 1:2) {
  if (posi == 1) {
    label_pos <- 'NEG'
  } else {
    label_pos <- 'POS'
  }

  file_select <- paste0('Results_', M_version, '/Combine/AUCgreaterThan', 100 * AUC_low, '_', label_pos, '_both.csv')
  file_or <- paste0('Data/plasma/', label_pos, '-Mean.csv')

  data_select <- read.csv(file_select)
  data_original <- read.csv(file_or)
  data_original <- data_original %>% filter(MS2.name != '')

  filename <- paste0('Results_', M_version, '/Combine/AUC_', label_pos, '.pdf')

  list <- data_select$MS2.nam
  p <- ggplot(data_select) +
    geom_segment(aes(x = factor(MS2.name, levels = list), xend = factor(MS2.name, levels = list),
                     y = 0, yend = ifelse(AUC < AUC.1, AUC, AUC.1)), color = 'grey90') +
    geom_segment(aes(x = factor(MS2.name, levels = list), xend = factor(MS2.name, levels = list),
                     y = AUC, yend = AUC.1), color = 'grey90') +
    geom_point(aes(x = factor(MS2.name, levels = list), y = AUC, color = "Plasma"), size = 6) +
    geom_point(aes(x = factor(MS2.name, levels = list), y = AUC.1, color = "Tissue"), size = 6) +
    scale_color_manual(values = c("Plasma" = '#4584B3', "Tissue" = '#E22985'), 
                       name = "Source", 
                       breaks = c("Plasma", "Tissue"), 
                       labels = c("Plasma", "Tissue")) +
    geom_hline(yintercept = 0.65, linetype = "dashed", color = "grey") +
    theme_classic() + 
    theme(panel.grid.major.x = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA, size = 1),
          axis.ticks.x = element_blank(),
          axis.text.x = element_text(size = 15, colour = "black"),
          axis.text.y = element_text(size = 15, colour = "black"),
          legend.title = element_text(),
          legend.text = element_text(size = 12)) + 
    xlab("") +
    ylab("AUC") +
    coord_flip()  

  ggsave(filename, plot = p, width = 12, height = 6)
}
save.image("./Result_1105.Rdata")

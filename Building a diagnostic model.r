#### Building a diagnostic model ####

library(randomForest)
library(pROC)
library(ggplot2)
library(dplyr)

load("./Result_1105.Rdata")

data_pos <- read.csv('Data/plasma/POS-Mean.csv') 
data_neg <- read.csv('Data/plasma/NEG-Mean.csv') 

auc_both_pos <- read.csv(paste0('Results_', M_version, '/Combine/AUCgreaterThan65_POS_both.csv'))
auc_both_neg <- read.csv(paste0('Results_', M_version, '/Combine/AUCgreaterThan65_NEG_both.csv'))

data_neg_s <- data_neg %>% filter(MS2.name %in% auc_both_neg$MS2.name)
data_pos_s <- data_pos %>% filter(MS2.name %in% auc_both_pos$MS2.name)

data_pos_V <- as.matrix(data_pos_s %>% select(VHL1:VHL9))
rownames(data_pos_V) <- data_pos_s$MS2.name
data_pos_C <- as.matrix(data_pos_s %>% select(CON1:CON9))
rownames(data_pos_C) <- data_pos_s$MS2.name
data_neg_V <- as.matrix(data_neg_s %>% select(VHL1:VHL9))
rownames(data_neg_V) <- data_neg_s$MS2.name
data_neg_C <- as.matrix(data_neg_s %>% select(CON1:CON9))
rownames(data_neg_C) <- data_neg_s$MS2.name

data_pos_VC <- rbind(data_pos_V, data_pos_C) %>% t()
data_neg_VC <- rbind(data_neg_V, data_neg_C) %>% t()

data_all_VC <- cbind(as.data.frame(data_pos_VC, data_neg_VC), Group = c(rep(1, nrow(data_neg_V)), rep(0, nrow(data_neg_C))))

set.seed(184)
tr_ratio <- 0.8
NV_tr <- floor(nrow(data_neg_V) * tr_ratio)
NC_tr <- floor(nrow(data_neg_C) * tr_ratio)
tr_V <- sample(1:nrow(data_neg_V), NV_tr)
tr_C <- sample(1:nrow(data_neg_C), NC_tr)
train <- c(tr_V, nrow(data_neg_V) + tr_C)

fit.rf <- randomForest(Group~., data = data_all_VC[train, ], mtry = 13, importance = TRUE, replace = TRUE)

Y_pre <- predict(fit.rf, newdata = data_all_VC[-train, ])
roc_obj <- roc(data_all_VC[-train, ]$Group, Y_pre, ci = TRUE)
auc(roc_obj)

pdf(paste0('Results_', M_version, '/Diagnose_AUC.pdf'))
plot(roc_obj, print.auc = TRUE, print.thres = TRUE, identity.col = 'green')
dev.off()

imp = importance(fit.rf)
rownames(imp) <- colnames(data_all_VC)[1:23]  
xmin <- min(imp[, 1])
dotchart(imp[order(-imp[, 1]), 1], xlab = colnames(imp)[1], ylab = "")

filename <- paste0('Results_', M_version, '/Importance_lolipop.pdf')
ggsave(filename, width = 12, height = 6)

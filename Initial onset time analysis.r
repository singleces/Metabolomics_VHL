#### Analysis of Initial Onset Time by Different Parts ####
# Load necessary libraries
library(dplyr)
library(survival)
library(survminer)

names(C_index_part) <- c('CNS','RA','RCC','PCT','Pheo','Genital')
for(part in c(1,3:4)){
  if(part==1){
    part_label <- 'CNS'
    dat_VHL_A <- dat_VHL_All %>% mutate(status = !is.na(CNS),
                                        time = ifelse(status, as.numeric(CNS), 2023 - as.numeric(Birth.year)))
  }# CNS
  else if(part==2){
    part_label <- 'RA'
    dat_VHL_A <- dat_VHL_All %>% mutate(status = !is.na(RA),
                                        time = ifelse(status, as.numeric(RA), 2023 - as.numeric(Birth.year)))
  }# RA
  else if(part==3){
    part_label <- 'RCC'
    dat_VHL_A <- dat_VHL_All %>% mutate(status = !is.na(RCC),
                                        time = ifelse(status, as.numeric(RCC), 2023 - as.numeric(Birth.year)))
  }# RCC
  else if(part==4){
    part_label <- 'PCT'
    dat_VHL_A <- dat_VHL_All %>% mutate(status = !is.na(PCT),
                                        time = ifelse(status, as.numeric(PCT), 2023 - as.numeric(Birth.year)))
  }# PCT
  else if(part==5){
    part_label <- 'Pheo'
    dat_VHL_A <- dat_VHL_All %>% mutate(status = !is.na(Pheo),
                                        time = ifelse(status, as.numeric(Pheo), 2023 - as.numeric(Birth.year)))
  }# Pheo
  else if(part==6){
    part_label <- 'Genital'
    dat_VHL_A <- dat_VHL_All %>% mutate(status = !is.na(Shengzhi),
                                        time = ifelse(status, as.numeric(Shengzhi), 2023 - as.numeric(Birth.year)))
  }# Genital
  # combine data_selected to data_cox
  dat_VHL <- dat_VHL_A %>% select(Family.ID, Gender, BY, Family.history,
                                  PM_Exon1, PM_Exon2, PM_Exon3, typeIIA,
                                  typeIIB, Generation, time, status)

  if(part!=3){
    data_selected = data_selected_all
  }else{
    data_selected = data_selected_both
  }
  index_id <- match(dat_VHL$Family.ID, rownames(data_selected))
  data_selected <- data_selected[index_id,]
  data_s_01 <- data_selected %>% mutate_all(function(x){as.factor(x>median(x))})
  data_cox <- data.frame(dat_VHL, data_s_01)
  data_temp <- as.data.frame(lapply(data_cox[,-1],as.numeric))
  Family.ID <- data_cox$Family.ID
  data_cox <- data.frame(Family.ID, data_temp)

  if(test_Set){
    t_ratio <- 0.2
    set.seed(1234)
    NN <- nrow(data_cox)
    NN_test <- floor(NN*t_ratio)
    te <- sample(1:NN, NN_test)
    dat_cox <- data_cox[-te,]
    dat_cox_test <- data_cox[te,-1]
    C_index_part <- rep(NA,6)
  }
  M <- ncol(dat_cox)
  Cox_name <- colnames(data_cox)[c(2:10, 13:M)]
  colnames(dat_cox)[c(2:10, 13:M)] <- paste('X',1:(M-3),sep='')
  result <- data.frame(matrix(NA, ncol = 5, nrow = M - 3))
  for(ii in 1:(M-3)){
    X_name <- paste0('X',ii)
    formu <- as.formula(paste0('Surv(time, status) ~ X', ii))
    coxmodel <- coxph(formu,data = dat_cox)
    result[ii,c(1:5)] <- summary(coxmodel)$coefficients
  }
  rownames(result) <- Cox_name
  colnames(result) <- c("coef","exp(coef)","se(coef)","z","Pr(>|z|)")
  result2 <- data.frame(matrix(NA, ncol = 5, nrow = 2))
  data_cox2 <- mutate(data_cox, Deletion = 1 - (PM_Exon1 + PM_Exon2 + PM_Exon3),
                      typeI = 1 - (typeIIA + typeIIB))
  coxmodel2 <- coxph(Surv(time, status) ~ Deletion, data = data_cox2)
  result2[1,] <- summary(coxmodel2)$coefficients
  coxmodel2 <- coxph(Surv(time, status) ~ typeI, data = data_cox2)
  result2[2,] <- summary(coxmodel2)$coefficients
  rownames(result2) <- c('Deletion','typeI')
  colnames(result2) <- c("coef","exp(coef)","se(coef)","z","Pr(>|z|)")
  result_all <- rbind(result, result2)
  filename <- paste0('Results_',M_version,'/Surv/',part_label,'/Cox.csv')
  write.csv(result_all, file=filename)
  
  # multi variable cox
  name_inx <- Cox_name[result[,5] < 0.1]
  cox_data <- select(data_cox, time, status, all_of(name_inx))
  if(test_Set){
    cox_data <- select(data_cox[-te,], time, status, all_of(name_inx))
    multiCox_model <- coxph(Surv(time, status) ~ ., cox_data)
    multiCox_model <- MASS::stepAIC(multiCox_model)
    res.cox <- coxph(multiCox_model$formula, data =cox_data)
    C <- rcorrcens(Surv(time,status) ~ predict(res.cox, dat_cox_test), data = dat_cox_test)
    C_index_part[part] <- 1- C[1]
    cox <- coxph(multiCox_model$formula, data =cox_data) 
    filename_forest <- paste0('Results_', M_version, '/Surv/',part_label,'/Cox_forest.pdf')
    ggforest(cox, data = cox_data, main = "Hazard ratio")
    ggsave(filename_forest)
  }else{
    cox_data <- select(data_cox, time, status, all_of(name_inx))
    multiCox_model <- coxph(Surv(time, status) ~ ., cox_data)
    multiCox_model <- MASS::stepAIC(multiCox_model)
    cox <- coxph(multiCox_model$formula, data =cox_data) 
    filename_forest <- paste0('Results_', M_version, '/Surv/',part_label,'/Cox_forest.pdf')
    ggforest(cox, data = cox_data, main = "Hazard ratio")
    ggsave(filename_forest)

  }
  result_multicox <- as.data.frame(summary(multiCox_model)$coefficients)
  if(ncol(result_multicox) != 0){
    filename_multicox <- paste0('Results_', M_version, '/Surv/',part_label,'/multiCox.csv')
    write.csv(result_multicox, file=filename_multicox)
  }
  for(ii in 1:(M-3)){
    if(Cox_name[ii] %in% rownames(result_multicox)){
      if(ii >= 10){
        lname = Cox_name[ii]
        llab = c('Order1', 'Others')
      }else{
        if(ii == 2){
          lname = 'Birth year'
          llab = c('<=1980', '>1980')
        }else{
          if(ii == 1){
            lname = 'Gender'
            llab = c('Male', 'Female')
          }else{
            lname = 'Type'
            llab = c(Cox_name[ii], 'Others')
          }
        }
      }
      X_name <- paste0('X',ii)
      formu <- as.formula(paste0('Surv(time, status) ~ X', ii))
      fit <- surv_fit(formu,data = dat_cox)
      filename_surv <- paste0('Results_', M_version, '/Surv/',part_label,'/KM_',ii,'.pdf')
      ggsurvplot(fit,  fun = "event",
                 legend.tiltle = lname,
                 legend.labs = llab, title = paste0(part_label,'_',lname),
                 ylab='Penetrance(%)',xlab = 'Age (Years)',
                 pval = TRUE, ggtheme = theme_bw(),
                 palette = c('red', 'blue'))
      ggsave(filename_surv, width = 4, height=3)
    }
  }
}

if(test_Set){
  C_index_all <- c(C_index_A, C_index_part)
  names(C_index_all)[1] <- 'ALL'
  write.csv(C_index_all, file = paste0('Results_', M_version, '/Surv/C_index.csv'))
  C_index_all
}



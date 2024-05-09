#### Comparative box plots ####

# Load necessary libraries
library(ggplot2)
library(dplyr)
library(ggpubr)

{
  # Combining Analysis: Merging AUC data from plasma and tissue samples
  filename_plasma = paste0('Results_', M_version, '/plasma/AUC', 100 * AUC_low, '_', label_pos, '.csv')
  filename_tissue = paste0('Results_', M_version, '/tissue/AUC', 100 * AUC_low, '_', label_pos, '.csv')
  
  auc_plasma <- read.csv(filename_plasma)
  auc_tissue <- read.csv(filename_tissue)
  
  idx1 <- match(auc_tissue$MS2.name, auc_plasma$MS2.name)
  auc_both <- data.frame(auc_plasma[idx1, -1]) 
  auc_both <- auc_both %>% filter(!is.na(MS2.name))  
  
  idx2 <- match(auc_both$MS2.name, auc_tissue$MS2.name)
  tmp <- data.frame(auc_tissue[idx2, -(1:2)])  
  auc_both <- data.frame(auc_both, tmp)  

  filename_both <- paste0("Results_", M_version, "/Combine/AUCgreaterThan", 100 * AUC_low, "_", label_pos, "_both.csv")
  write.csv(auc_both, file = filename_both)


  for (ppi in 1:nrow(auc_both)) {
    p_name <- auc_both$MS2.name[ppi]
    
    safe_p_name <- gsub("[[:punct:]]", "_", p_name)  

    Title_box <- paste0("Box Plots of ", p_name)
    filename_box = paste0('Results_', M_version, '/Combine/', label_pos, '_Boxboth_', ppi, '_', safe_p_name, '.pdf')

    dat_boxplot = data_plasma[, c(p_name, 'GroupFc')]
    colnames(dat_boxplot) = c('Intensity', 'Group')
    dat_boxplot$Intensity <- as.numeric(dat_boxplot$Intensity)
    dat_boxplot$Group <- as.factor(dat_boxplot$Group)

    dat_boxplot2 = data_tissue[, c(p_name, 'GroupFc')]
    colnames(dat_boxplot2) = c('Intensity', 'Group')
    dat_boxplot2$Intensity <- as.numeric(dat_boxplot2$Intensity)
    dat_boxplot2$Group <- as.factor(dat_boxplot2$Group)
    
    col1 = c("#4584B3", "#E22985")  # Colors for plasma
    col2 = c('blue', "#e31a1c")  # Colors for tissue

    # Create the plasma box plot
    p1 = ggplot(dat_boxplot, aes(x = Group, y = Intensity, color = Group)) + 
      geom_violin(alpha = 0.4, size = 0.5, width = 0.6) +
      geom_boxplot(position = "dodge", alpha = 1, outlier.size = 0.3, size = 0.5, width = 0.2) +
      theme_classic() +
      labs(title = Title_box) +
      geom_jitter(size = 1, width = 0.15, alpha = 0.5) +
      scale_color_manual(values = col1) +
      theme(legend.position = "none") +
      theme(axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8))

    # Create the tissue box plot
    p2 = ggplot(dat_boxplot2, aes(x = Group, y = Intensity, color = Group)) +
      geom_violin(alpha = 0.4, size = 0.5, width = 0.6) +
      geom_boxplot(position = "dodge", alpha = 1, outlier.size = 0.3, size = 0.5, width = 0.2) +
      theme_classic() +
      labs(title = "Tissue") +
      geom_jitter(size = 1, width = 0.2, alpha = 0.5) +
      scale_color_manual(values = col2) +
      theme(legend.position = "none")

    # Combine the plots side by side
    p <- p1 | p2
    p <- p + plot_annotation(title = Title_box, theme = theme(plot.title = element_text(hjust = 0.5)))

    # Save the combined plot to a PDF
    ggsave(filename_box, width = 6, height = 5)
  }
}


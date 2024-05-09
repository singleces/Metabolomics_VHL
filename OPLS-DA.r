# Load necessary libraries
library(ggplot2)
library(dplyr)
library(ropls)

# Set a random seed for reproducibility
set.seed(1)

# Data preprocessing
P_M <- 200  # Set the number of permutation tests
dataMatrix <- data_C_V %>% 
  select(-c('Group'))  
GroupFc <- data_C_V$Group  

# Perform OPLS-DA analysis
sacurine.plsda <- opls(
  dataMatrix, 
  y = GroupFc, 
  predI = 1, 
  orthoI = 1, 
  log10L = TRUE, 
  permI = P_M
)

# Compute sample scores
sample.score <- as.data.frame(sacurine.plsda@scoreMN) %>%
  mutate(
    GroupFc = data_C_V$Group, 
    o1 = sacurine.plsda@orthoScoreMN[,1]
  )

# Plot scores 
p <- ggplot(sample.score, aes(p1, o1, color = GroupFc)) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", size = 0.5) +
  geom_point() +
  labs(x = 'P1', y = 'to1') +
  stat_ellipse(level = 0.95, linetype = 'solid', size = 1, show.legend = FALSE) +
  scale_color_manual(values = c('#6BA0C5', '#DB706F')) +
  theme_bw() +
  theme(
    legend.position = c(0.15, 0.85),
    legend.title = element_blank(),
    legend.text = element_text(color = 'black', size = 12, face = 'plain'),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.text = element_text(color = 'black', size = 15, face = 'plain'),
    axis.title = element_text(color = 'black', size = 15, face = 'plain'),
    axis.ticks = element_line(color = 'black')
  )
print(p)

# Save the plot 
filename_OPLSS <- paste0('Results_', M_version, '/', label_type, '/OPLS_Score_', label_pos, '.pdf')
ggsave(filename_OPLSS, width = 5, height = 4)

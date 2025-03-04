# Load necessary libraries
library(ggplot2)
library(reshape2)
library(dplyr)

# Read the data
data <- read.csv("/data/WHRI-Bioinformatics/Nicholls/BP_Gene_Prioritisation/BP_Gene_Prior_v2/2_machine_learning/Multiclass/Output/fold_balanced_accuracies.csv")

# Melt the data for ggplot2
data_melted <- melt(data, id.vars = 'X', variable.name = 'Fold', value.name = 'Balanced_Accuracy')

# Rename the 'X' column to 'Model' for clarity
data_melted <- data_melted %>%
  rename(Model = X)

# Create the violin plot with updated style and features
plot <- data_melted %>%
  ggplot(aes(x = Model, y = Balanced_Accuracy, color = Model)) + 
  theme_bw() +
  geom_violin(trim = FALSE) +
  geom_point(position = position_dodge(width = 0.9), size = 1.5) +  # Use position_dodge for alignment
  scale_color_discrete(name = "Model") +
  theme(
    axis.text.x = element_text(color = "grey20", size = 16, hjust = .5, vjust = .5, face = "plain"),
    axis.text.y = element_text(color = "grey20", size = 14, angle = 0, hjust = 1, vjust = 0, face = "plain"),  
    axis.title.x = element_text(color = "grey20", size = 16, angle = 0, hjust = .5, vjust = 0, face = "plain"),
    axis.title.y = element_text(color = "grey20", size = 16, angle = 90, hjust = .5, vjust = .5, face = "plain"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14)
  ) +
  ggtitle("Model Performance Across Folds") + 
  xlab("Model") +
  ylab("Balanced Accuracy")

# Save the plot
ggsave("/data/WHRI-Bioinformatics/Nicholls/BP_Gene_Prioritisation/BP_Gene_Prior_v2/2_machine_learning/Multiclass/Output/model_performance_violin_plot.png",
       plot = plot, width = 18, height = 8, dpi = 300, limitsize = FALSE)

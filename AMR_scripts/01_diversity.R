# Set current directory to source file location
current_working_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(current_working_dir)

# Load necessary libraries
library(ggplot2)
library(reshape2)
library(tidyverse)
library(ggsignif)
library(vegan)
library(dunn.test)
library(ggpubr)

# import data
kma_out <- read.csv("kma_out.csv", check.names = FALSE)

# import meta
kma_meta <- read.csv("Metafile.csv", check.names = FALSE)

# add info on region_grouping (1 – Rural; 0 – Urban)
kma_out <- kma_out %>%
  left_join(
    kma_meta %>% dplyr::select(Group, Group.Bi) %>% unique(),
    by = c("tribe" = "Group")
  )
kma_out <- kma_out %>% 
  mutate(Group.Bi = factor(Group.Bi, levels = c(0, 1), labels = c("Urban", "Rural"))) %>% 
  relocate(Group.Bi, .after = tribe)

# Distribution (Tribe) ----

# Perform Shapiro-Wilk test on log-transformed values for each tribe
for (tribe in unique(kma_out$tribe)) {
  tribe_data <- kma_out$RPKM[kma_out$tribe == tribe]
  cat("Shapiro-Wilk test for RPKM in", tribe, ":\n")
  print(shapiro.test(tribe_data))
  cat("\n")
}

# kruskal willis
kruskal_result <- kruskal.test(RPKM ~ tribe, data = kma_out)
print(kruskal_result)

# Extracting the chi-squared statistic and p-value
kruskal_chisq <- kruskal_result$statistic
kruskal_pvalue <- kruskal_result$p.value

# dunn test
dunn_test_result <- dunn.test(kma_out$RPKM, kma_out$tribe, method = "bonferroni")

# re-arrange tribe according to least urban to most urban
kma_out$tribe <- factor(kma_out$tribe, levels = c("Jahai", "Temiar", "Temuan", "Malay"))

# visualization
rpkm_boxPlot_tribes <- ggplot(kma_out, aes(x = tribe, y = RPKM, fill = tribe)) +
  geom_boxplot() +
  scale_y_log10() +
  theme_classic() +
  scale_fill_manual(values = c("Jahai" = "darkgreen", "Temiar" = "skyblue",
                               "Temuan" = "orange", "Malay" = "pink")) +
  theme(axis.text.x = element_text()) +
  labs(title = "Boxplot of RPKM (log scale) by Tribe", y = "RPKM (log scale)", x = "Sample ID") +
  ggsignif::geom_signif(
    comparisons = list(c("Temiar", "Temuan"), 
                       c("Jahai", "Temiar"), 
                       c("Temiar", "Malay")),
    map_signif_level = TRUE,
    textsize = 3,
    y_position = c(7, 6, 6)
  )

plot1 <- rpkm_boxPlot_tribes + 
  annotate("text", x = 1, y = 1e09, 
           label = paste("Kruskal-Wallis chi-squared = ", round(kruskal_chisq, 2),
                         "\np-value = ", format(kruskal_pvalue, scientific = TRUE)), 
           hjust = 0, vjust = 1.5, size = 3.2, color = "red")
plot1

# Save as PNG
# ggsave("resfinder_plot1_distributionBoxPlot.png", plot = plot1, width = 10, height = 6, units = "in")

# Save as PDF
# ggsave("resfinder_plot1_distributionBoxPlot.pdf", plot = plot1, width = 10, height = 6, units = "in")

# Shannon Index (Tribe)----

## data preparation
# subset the data
df <- kma_out[c("sampleID", "refSequence", "RPKM")]

# Ensure the RPKM column is numeric
df$RPKM <- as.numeric(df$RPKM)
data_pivot <- pivot_wider(df, names_from = sampleID, values_from = RPKM, values_fill = list(RPKM = 0))

# Convert tibble to data frame
data_pivot <- as.data.frame(data_pivot)

# Set row names to the values in the first column (`refSequence`)
rownames(data_pivot) <- data_pivot[, 1]

# Remove the `refSequence` column since it's now the row names
data_pivot <- data_pivot[, -1]
data_matrix <- as.matrix(data_pivot)

# Transpose the matrix to calculate Shannon index per sample
data_matrix_t <- t(data_matrix)

# Calculate Shannon index per sample
shannon_index <- vegan::diversity(data_matrix_t)

# Add sample names back to the results
names(shannon_index) <- colnames(data_matrix)

# Convert to dataframe
shannon_index_df <- data.frame(
  sampleID = names(shannon_index),
  Shannon_Index = shannon_index
)

# Function to determine tribe based on sampleID prefix
get_tribe <- function(sampleID) {
  prefix <- substring(sampleID, 1, regexpr("[0-9]", sampleID) - 1)
  tribe <- switch(prefix,
                  "J"   = "Jahai", "TRk" = "Temiar",
                  "TM"  = "Temuan", "MLY" = "Malay")
  return(tribe)
}

# add tribe column
shannon_index_df$tribe <- sapply(shannon_index_df$sampleID, get_tribe)

## Normality test
for (tribe in unique(shannon_index_df$tribe)) {
  tribe_data <- shannon_index_df$Shannon_Index[shannon_index_df$tribe == tribe]
  cat("Shapiro-Wilk test for Shannon_Index in", tribe, ":\n")
  print(shapiro.test(tribe_data))
  cat("\n")
}

# Perform one-way ANOVA
anova_result <- aov(Shannon_Index ~ tribe, data = shannon_index_df)
summary(anova_result)

# adhoc test
TukeyHSD(anova_result)

# re-arrange tribe according to least urban to most urban
shannon_index_df$tribe <- factor(shannon_index_df$tribe,
                                 levels = c("Jahai", "Temiar", "Temuan", "Malay"))

# Create the boxplot
boxplot <- ggplot(shannon_index_df, aes(x = tribe, y = Shannon_Index, fill = tribe)) +
  geom_boxplot(outlier.shape = NA) +
  labs(title = "Shannon Index Distribution by Group",
       x = "Group",
       y = "Shannon Index",
       fill = "Group") +
  scale_fill_manual(values = c("Jahai" = "darkgreen", "Temiar" = "skyblue",
                               "Temuan" = "orange", "Malay" = "pink")) +
  geom_jitter(width = 0.2, color = "blue") +
  theme_classic() +
  theme(axis.text.x = element_text()) +
  ggsignif::geom_signif(
    comparisons = list(c("Jahai", "Temuan")),
    annotations = "p = 0.042",
    map_signif_level = TRUE,
    textsize = 3.5,
    y_position = c(3, 3.2)
  ) +
  annotate("text", x = 1.5, y = 3.5, 
           label = "ANOVA: p = 0.0402", size = 4, color = "red")
boxplot

# Save as PNG
# ggsave("resfinder_plot2_shannonBoxPlot_revised.png", plot = boxplot, width = 10, height = 6, units = "in")

# Save as PDF
# ggsave("resfinder_plot2_shannonBoxPlot_revised.pdf", plot = boxplot, width = 10, height = 6, units = "in")

# Distribution (Region Group) ----

# Perform Shapiro-Wilk test on log-transformed values for each group
for (Group.Bi in unique(kma_out$Group.Bi)) {
  Group.Bi_data <- kma_out$RPKM[kma_out$Group.Bi == Group.Bi]
  cat("Shapiro-Wilk test for RPKM in", Group.Bi, ":\n")
  print(shapiro.test(Group.Bi_data))
  cat("\n")
}

# perform wilcoxon test since not normal dist
wilcox_result <- wilcox.test(RPKM ~ Group.Bi, data = kma_out)
wilcox_result

# keep value for visualization
wilcox_W <- wilcox_result$statistic
wilcox_p <- wilcox_result$p.value

# visualization
rpkm_boxPlot_Group <- ggplot(kma_out, aes(x = Group.Bi, y = RPKM, fill = Group.Bi)) +
  geom_boxplot() +
  scale_y_log10() +
  theme_classic() +
  scale_fill_manual(values = c("Rural" = "darkgreen", "Urban"= "pink")) +
  theme(axis.text.x = element_text()) +
  labs(title = "Boxplot of RPKM (log scale) by Region Group", 
       y = "RPKM (log scale)", 
       x = "Region Group",
       fill = "Region Group")
rpkm_boxPlot_Group

rpkm_boxPlot_Group_v2 <- rpkm_boxPlot_Group + 
  annotate("text", x = 1.2, y = 1e06, 
           label = paste0("Wilcoxon rank-sum test\n",
                          "W = ", wilcox_W,
                          "\n",
                          "p = ", formatC(wilcox_p, format = "e", digits = 2)),
           hjust = 0, vjust = 1.5, size = 3.2, color = "red")
rpkm_boxPlot_Group_v2

# Save as PNG
# ggsave("resfinder_plot3_distributionBoxPlot_regionGroup.png", plot = rpkm_boxPlot_Group_v2, width = 10, height = 6, units = "in")

# Save as PDF
# ggsave("resfinder_plot3_distributionBoxPlot_regionGroup.pdf", plot = rpkm_boxPlot_Group_v2, width = 10, height = 6, units = "in")

# Shannon Index (Region Group)----

## data preparation
# subset the data
df <- kma_out[c("sampleID", "refSequence", "RPKM")]

# Ensure the RPKM column is numeric
df$RPKM <- as.numeric(df$RPKM)
data_pivot <- pivot_wider(df, names_from = sampleID, values_from = RPKM, values_fill = list(RPKM = 0))

# Convert tibble to data frame
data_pivot <- as.data.frame(data_pivot)

# Set row names to the values in the first column (`refSequence`)
rownames(data_pivot) <- data_pivot[, 1]

# Remove the `refSequence` column since it's now the row names
data_pivot <- data_pivot[, -1]
data_matrix <- as.matrix(data_pivot)

# Transpose the matrix to calculate Shannon index per sample
data_matrix_t <- t(data_matrix)

# Calculate Shannon index per sample
shannon_index <- vegan::diversity(data_matrix_t)

# Add sample names back to the results
names(shannon_index) <- colnames(data_matrix)

# Convert to dataframe
shannon_index_df <- data.frame(
  sampleID = names(shannon_index),
  Shannon_Index = shannon_index
)

# Function to determine tribe based on sampleID prefix
get_tribe <- function(sampleID) {
  prefix <- substring(sampleID, 1, regexpr("[0-9]", sampleID) - 1)
  tribe <- switch(prefix,
                  "J"   = "Jahai", "TRk" = "Temiar",
                  "TM"  = "Temuan", "MLY" = "Malay")
  return(tribe)
}

# add tribe column
shannon_index_df$tribe <- sapply(shannon_index_df$sampleID, get_tribe)

# add info on region_grouping (1 – Rural; 0 – Urban)
shannon_index_df <- shannon_index_df %>%
  left_join(
    kma_meta %>% dplyr::select(Group, Group.Bi) %>% unique(),
    by = c("tribe" = "Group")
  )
shannon_index_df <- shannon_index_df %>% 
  mutate(Group.Bi = factor(Group.Bi, levels = c(0, 1), labels = c("Urban", "Rural")))

# shapiro test
for (group in unique(shannon_index_df$Group.Bi)) {
  Group.Bi_data <- shannon_index_df$Shannon_Index[shannon_index_df$Group.Bi == group]
  cat("Shapiro-Wilk test for Shannon Index in", Group.Bi, ":\n")
  print(shapiro.test(Group.Bi_data))
  cat("\n")
}

# perform t test since Shannon index passed normal dist
t_test <- t.test(Shannon_Index ~ Group.Bi, data = shannon_index_df)
t_test

# save value for viz
p_val <- signif(t_test$p.value, 3)

# visualization
set.seed(888) # keep jitter same
shannon_boxPlot_Group <- ggplot(shannon_index_df, aes(x = Group.Bi, y = Shannon_Index, fill = Group.Bi)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, color = "blue") +
  theme_classic() +
  scale_fill_manual(values = c("Rural" = "darkgreen", "Urban"= "pink")) +
  theme(axis.text.x = element_text()) +
  labs(title = "Boxplot of Shannon Index by Region Group", 
       y = "Shannon Index", 
       x = "Region Group",
       fill = "Region Group")
shannon_boxPlot_Group

set.seed(888) # keep jitter same
shannon_boxPlot_Group_v2 <- shannon_boxPlot_Group + 
  annotate("text", x = 1.5, y = 3, 
           label = paste("t-test, p =", p_val), 
           size = 4, color = "red")
shannon_boxPlot_Group_v2

# Save as PNG
# ggsave("resfinder_plot4_shannonBoxPlot.png", plot = shannon_boxPlot_Group_v2, width = 10, height = 6, units = "in")

# Save as PDF
# ggsave("resfinder_plot4_shannonBoxPlot.pdf", plot = shannon_boxPlot_Group_v2, width = 10, height = 6, units = "in")

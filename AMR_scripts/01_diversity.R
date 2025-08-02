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
library(scales)

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

# Perform Shapiro-Wilk test on RPKM of each tribe
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
data.frame(
  comparison = dunn_test_result$comparisons,
  p_value = dunn_test_result$P.adjusted
)

# re-arrange tribe according to least urban to most urban
kma_out$tribe <- factor(kma_out$tribe, levels = c("Jahai", "Temiar", "Temuan", "Malay"))

# visualization
rpkm_plot_groups <- ggplot(kma_out, aes(x = tribe, y = RPKM, fill = tribe)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1, fill = "white") +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_classic(base_family = "sans") +
  theme(legend.position = "none",
        plot.title = element_text(
          face = "bold",
          size = 15),
        text = element_text(size = 16)) +
  scale_fill_manual(values = c("Jahai" = "darkgreen", "Temiar" = "skyblue",
                               "Temuan" = "orange", "Malay" = "pink")) +
  labs(title = "A",
       x = NULL,
       y = "ARG Load (Log scaled RPKM)" 
       ) +
  ggsignif::geom_signif(
    comparisons = list(c("Temiar", "Temuan"),
                       c("Jahai", "Temiar"),
                       c("Temiar", "Malay")),
    map_signif_level = TRUE,
    annotations = c("p = 3.0484e-05",
                    "p = 1.8472e-04",
                    "p = 1.6855e-03"),
    textsize = 4,
    y_position = c(7, 6.5, 6.5)
  )
rpkm_plot_groups

# Save as PNG
# ggsave("resfinder_rpkmPlot_groups.png", plot = rpkm_plot_groups, width = 10, height = 6, units = "in")

# Save as PDF
# ggsave("resfinder_rpkmPlot_groups.pdf", plot = rpkm_plot_groups, width = 10, height = 6, units = "in")

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
shannon_plot_groups <- ggplot(shannon_index_df, aes(x = tribe, y = Shannon_Index, fill = tribe)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.05, fill = "white", outlier.shape = NA) +
  theme_classic(base_family = "sans") +
  theme(legend.position = "none",
        plot.title = element_text(
          face = "bold",
          size = 15),
        text = element_text(size = 16)) +
  labs(title = "B",
       x = NULL,
       y = "ARG Diversity (Shannon Index)",
       ) +
  scale_fill_manual(values = c("Jahai" = "darkgreen", "Temiar" = "skyblue",
                               "Temuan" = "orange", "Malay" = "pink")) +
  ggsignif::geom_signif(
    comparisons = list(c("Jahai", "Temuan")),
    annotations = "p = 0.042",
    map_signif_level = TRUE,
    textsize = 4,
    y_position = c(4.3, 3.2)
  )
shannon_plot_groups

# Save as PNG
# ggsave("resfinder_shannonPlot_groups.png", plot = shannon_plot_groups, width = 10, height = 6, units = "in")

# Save as PDF
# ggsave("resfinder_shannonPlot_groups.pdf", plot = shannon_plot_groups, width = 10, height = 6, units = "in")

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

# adjust sequence for plot
kma_out$Group.Bi <- factor(kma_out$Group.Bi, levels = c("Rural", "Urban"))

# visualization
rpkm_plot_Region <- ggplot(kma_out, aes(x = Group.Bi, y = RPKM, fill = Group.Bi)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.05, fill = "white", outlier.shape = NA) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_classic(base_family = "sans") +
  theme(legend.position = "none",
        plot.title = element_text(
          face = "bold",
          size = 15),
        text = element_text(size = 16)) +
  labs(title = "C",
       x = NULL,
       y = "ARG Load (Log scaled RPKM)" 
       # fill = "Region Groups"
       ) +
  ggsignif::geom_signif(
    comparisons = list(c("Urban", "Rural")),
    annotations = "p = 0.02707",
    map_signif_level = TRUE,
    textsize = 4,
    y_position = 6.5) +
  scale_fill_manual(values = c("Rural" = "darkgreen", "Urban"= "pink"))
rpkm_plot_Region

# Save as PNG
# ggsave("resfinder_rpkmPlot_Region.png", plot = rpkm_plot_Region, width = 10, height = 6, units = "in")

# Save as PDF
# ggsave("resfinder_rpkmPlot_Region.pdf", plot = rpkm_plot_Region, width = 10, height = 6, units = "in")

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
  cat("Shapiro-Wilk test for Shannon Index in", group, ":\n")
  print(shapiro.test(Group.Bi_data))
  cat("\n")
}

# perform t test since Shannon index passed normal dist
t_test <- t.test(Shannon_Index ~ Group.Bi, data = shannon_index_df)
t_test

# save value for viz
p_val <- signif(t_test$p.value, 3)

# adjust sequence for plot
shannon_index_df$Group.Bi <- factor(shannon_index_df$Group.Bi, levels = c("Rural", "Urban"))

# visualization
shannon_plot_Region <- ggplot(shannon_index_df, aes(x = Group.Bi, y = Shannon_Index, fill = Group.Bi)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.05, fill = "white", outlier.size = 0.5) +
  theme_classic(base_family = "sans") +
  theme(legend.position = "none",
        plot.title = element_text(
          face = "bold",
          size = 15),
        text = element_text(size = 16)) +
  labs(title = "D",
       x = NULL,
       y = "ARG Diversity (Shannon Index)" 
       # fill = "Region Groups"
  ) +
  ggsignif::geom_signif(
    comparisons = list(c("Urban", "Rural")),
    annotations = p_val,
    map_signif_level = TRUE,
    textsize = 4,
    y_position = 4.25) +
  scale_fill_manual(values = c("Rural" = "darkgreen", "Urban"= "pink"))
shannon_plot_Region

# Save as PNG
# ggsave("resfinder_shannonPlot_Region.png", plot = shannon_plot_Region, width = 10, height = 6, units = "in")

# Save as PDF
# ggsave("resfinder_shannonPlot_Region.pdf", plot = shannon_plot_Region, width = 10, height = 6, units = "in")

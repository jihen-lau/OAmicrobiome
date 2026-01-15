# Set current directory to source file location
current_working_dir <- dirname(rstudioapi::getActiveDocumentContext()$path) 
setwd(current_working_dir)

# Load necessary libraries
library(tidyverse) 
library(UpSetR)

# Data Preparation ----

# import data
kma_out <- read.csv("kma_out.csv", check.names = FALSE)

# UpSet Plot ----

# Split data by tribe and extract Gene name
gene_list <- kma_out %>%
  group_by(tribe) %>%
  summarise(refSequence = list(unique(refSequence))) %>% 
  deframe()

# make it a dataframe for plotting
gene_df <- fromList(gene_list)

# plot
# Save as pdf
# pdf("resfinder_upsetPlot.pdf", width = 14, height = 10)

# Uncomment the next line to export the plot
# png("resfinder_upsetPlot.png", width = 10, height = 5, units = "in", res = 300)

upset(gene_df, nsets = 4, point.size = 4, text.scale = 1.8, 
      sets = c("Jahai", "Temiar", "Temuan", "Malay"), 
      keep.order = TRUE, order.by = "freq",
      mainbar.y.label = "Number of ARGs in\nSet Intersections", 
      sets.x.label = "Number of ARGs")

# dev.off()

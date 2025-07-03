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
pdf("resfinder_plot5_upsetPlot.pdf", width = 14, height = 10)
upset(gene_df, nsets = 4, point.size = 4, text.scale = 1.8, 
      sets = c("Jahai", "Temiar", "Temuan", "Malay"), 
      keep.order = TRUE, order.by = "freq",
      mainbar.y.label = "Number of AMR Genes in\nSet Intersections", 
      sets.x.label = "Number of AMR Genes")
grid::grid.text("AMR Gene Distribution Across Groups", x = 0.6, y = 0.95, gp = grid::gpar(fontsize = 16))
dev.off()

# Set current directory to source file location
current_working_dir <- dirname(rstudioapi::getActiveDocumentContext()$path) 
setwd(current_working_dir)

# Load necessary libraries
library(tidyverse)
library(ggplot2)

# Data Preparation ----
# # Create a list of gene mapping data files in the current directory
# fileList <- list.files(pattern = "*.mapstat")
# 
# # Initialize an empty list to store all processed gene mapping data
# mapstat_conc <- lapply(fileList, function(file) {
#   # Extract sampleID from the file name
#   sampleID <- sub("\\.mapstat$", "", file)
#   
#   # Read the gene mapping data from file into a data frame
#   data <- read.delim(file, header = TRUE, check.names = FALSE, skip = 6)
#   
#   # Add sampleID column
#   data <- cbind(sampleID, data)
# })
# 
# # Combine all data frames into a single data frame
# mapstat_conc <- do.call(rbind, mapstat_conc)
# 
# # rename column
# colnames(mapstat_conc)[colnames(mapstat_conc) == "# refSequence"] <- "refSequence"
# 
# # add gene length
# resfinderFG_meta <- read.csv("resfinderFG_db.csv")
# mapstat_conc <- left_join(mapstat_conc, resfinderFG_meta, by = "refSequence")
# 
# # Split the first column into four new columns
# kma_resfinderFG <- mapstat_conc %>%
#   separate(refSequence, into = c("feature", "accession", "source", "antibiotic"), sep = "\\|")
# 
# # select columns
# kma_resfinderFG <- kma_resfinderFG %>% select(sampleID, feature, accession, source, antibiotic, readCount, Length)
# 
# # remove "_FG" suffix at sampleID
# kma_resfinderFG$sampleID <- gsub("_FG$", "", kma_resfinderFG$sampleID)
# 
# # Calculate total reads per sample
# sample_totalReads <- kma_resfinderFG %>%
#   group_by(sampleID) %>%
#   summarise(total_reads = sum(readCount, na.rm = TRUE)) %>%
#   ungroup()
# 
# # Join the total reads with the original dataframe
# kma_resfinderFG <- kma_resfinderFG %>%
#   left_join(sample_totalReads, by = "sampleID")
# 
# # Calculate RPKM
# kma_resfinderFG <- kma_resfinderFG %>%
#   mutate(RPKM = (readCount * 1e9) / (Length * total_reads))
# 
# # Function to determine tribe based on sampleID prefix
# get_tribe <- function(sampleID) {
#   prefix <- substring(sampleID, 1, regexpr("[0-9]", sampleID) - 1)
#   tribe <- switch(prefix,
#                   "J"   = "Jahai",
#                   "TRk" = "Temiar",
#                   "TM"  = "Temuan",
#                   "MLY" = "Malay")
#   return(tribe)
# }
# 
# # Add tribe column to the merged data frame
# kma_resfinderFG$tribe <- sapply(kma_resfinderFG$sampleID, get_tribe)
# kma_resfinderFG <- kma_resfinderFG %>% relocate("tribe", .after="sampleID")
# 
# head(kma_resfinderFG) %>% knitr::kable() 

# export to csv
# write.csv(kma_resfinderFG, file = "kma_resfinderFG.csv", row.names = FALSE)

# import data
kma_resfinderFG <- read.csv("kma_resfinderFG.csv", check.names = FALSE)

# Stacked Bar Chart ----

# calculate percentage within each antibiotic family
amr_summary <- kma_resfinderFG %>%
  group_by(antibiotic, tribe) %>%  
  summarise(
    total_RPKM = sum(RPKM, na.rm = TRUE)
  ) %>%
  group_by(antibiotic) %>%
  mutate(
    percentage = (total_RPKM / sum(total_RPKM)) * 100
  ) %>%
  ungroup()

amr_summary$tribe <- factor(amr_summary$tribe, levels = c("Malay", "Temiar", 
                                                          "Temuan", "Jahai"))

# Create a named vector for the labels
antibiotic_labels <- c(
  "AMC" = "AMC (Amoxicillin/Clavulanic acid)", "AMP" = "AMP (Ampicillin)",
  "AMX" = "AMX (Amoxicillin)", "ATM" = "ATM (Aztreonam)",
  "CAR" = "CAR (Carbenicillin)", "CAZ" = "CAZ (Ceftazidime)",
  "CHL" = "CHL (Chloramphenicol)", "CIP" = "CIP (Ciprofloxacin)",
  "CTX" = "CTX (Cefotaxime)", "CYC" = "CYC (Cycloserine)",
  "FEP" = "FEP (Cefepime)", "GEN" = "GEN (Gentamicin)",
  "KAN" = "KAN (Kanamycin)", "MIN" = "MIN (Minocycline)",
  "OXY" = "OXY (Oxytetracycline)", "PEN" = "PEN (Penicillin)",
  "PIP" = "PIP (Piperacillin)", "SIS" = "SIS (Sisomicin)",
  "SMZ" = "SMZ (Sulfamethoxazole)", "SPT" = "SPT (Spectinomycin)",
  "SXT" = "SXT (Trimethoprim/Sulfamethoxazole)", "TET" = "TET (Tetracycline)",
  "TGC" = "TGC (Tigecycline)", "TMP" = "TMP (Trimethoprim)"
)

# re-arrange tribe according to most urban to least urban
amr_summary$tribe <- factor(amr_summary$tribe, levels = c("Malay", "Temuan", "Temiar", "Jahai"))

resfinderFG_stackedBar <- ggplot(amr_summary, aes(x = antibiotic, 
                                                  y = percentage, 
                                                  fill = tribe)) + 
  geom_bar(stat = "identity") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.margin = margin(l = 60, r = 30, 
                             t = 20, b = 20,
                             unit = "pt")) +
  scale_fill_manual(values = c("Jahai" = "darkgreen", "Temiar" = "skyblue",
                               "Temuan" = "orange", "Malay" = "pink")) +
  scale_x_discrete(labels = antibiotic_labels) +
  labs(x = "Antibiotic Family", y = "Proportion in Total Pool",
       fill = "Tribe", title = "Distribution of functional genes in each tribes (in %)")

resfinderFG_stackedBar

# Save as PNG
# ggsave("resfinderFG_plot1_stackedBar.png", plot = resfinderFG_stackedBar, width = 10, height = 6, units = "in")

# Save as PDF
# ggsave("resfinderFG_plot1_stackedBar.pdf", plot = resfinderFG_stackedBar, width = 10, height = 8, units = "in")

# Heatmap ----

# Aggregate RPKM of individual antibiotic family within each sample
anFam_data <- kma_resfinderFG %>%
  group_by(antibiotic, sampleID) %>%  
  summarise(
    total_RPKM = sum(RPKM, na.rm = TRUE)
  ) %>%
  ungroup()

# Create a complete data frame with all gene families
complete_data <- anFam_data %>%
  complete(sampleID, antibiotic, fill = list(total_RPKM = 0))

# Spread data into a wide format
heatmap_data <- complete_data %>%
  pivot_wider(names_from = antibiotic, values_from = total_RPKM, values_fill = 0)

# Convert to matrix format and set row names
heatmap_matrix <- as.matrix(heatmap_data %>% select(-sampleID))
rownames(heatmap_matrix) <- heatmap_data$sampleID

# Log-transform data
heatmap_matrix_log <- log10(heatmap_matrix + 1)  # Adding 1 to avoid log10(0)
rownames(heatmap_matrix_log) <- heatmap_data$sampleID

# Transpose matrix
heatmap_matrix_log_T <- t(heatmap_matrix_log)

# Define sample order
ordered_sample_ids <- c(
  # Jahai samples
  "J02", "J04", "J09", "J10", "J14", "J16", "J18", "J24", "J26", "J31", "J38", "J41",
  # Temiar samples
  "TRk011F", "TRk025M", "TRk041F", "TRk064F", "TRk069F", "TRk094M", "TRk123F", "TRk125M", "TRk136M",
  # Temuan samples
  "TM016M", "TM017F", "TM037F", "TM039M", "TM056F", "TM114M", "TM123M", "TM125F", "TM167F", "TM168M", "TM169M", "TM175F",
  # Malay samples
  "MLY001", "MLY003", "MLY004", "MLY005", "MLY006", "MLY007", "MLY008", "MLY009", "MLY010"
)

# Reorder matrix columns to match sample order
heatmap_matrix_log_T <- heatmap_matrix_log_T[, ordered_sample_ids]

# Replace rownames using the lookup
rownames(heatmap_matrix_log_T) <- antibiotic_labels[rownames(heatmap_matrix_log_T)]

# Define tribe groups
tribes <- c(
  rep("Jahai", 12),
  rep("Temiar", 9), 
  rep("Temuan", 12), 
  rep("Malay", 9) 
)

# Column annotation for tribes
column_anno <- HeatmapAnnotation(
  foo = anno_empty(border = FALSE), # add space for degree of urbanisation
  Tribes = tribes,
  col = list(Tribes = c("Jahai" = "darkgreen",
                        "Temiar" = "skyblue",
                        "Temuan" = "orange",
                        "Malay" = "pink")),
  show_annotation_name = FALSE,
  show_legend = FALSE
)

# Create custom legends
tribes_legend <- Legend(
  labels = c("Jahai", "Temiar", "Temuan", "Malay"),
  legend_gp = gpar(fill = c("darkgreen", "skyblue", "orange", "pink")),
  title = "Group"
)

# Create color mapping for heatmap
breaks <- c(0, 1, 2, 3, 4, 5, 6)
color_gradient <- colorRampPalette(c("#ffcbd1", "lightpink", "#ee6b6e", "#ff2c2c", "#CC0000", "#990000", "#660000"))(length(breaks) - 1)
colors <- c("white", color_gradient)
col_fun <- colorRamp2(breaks, colors)

# Uncomment the next line to export the plot
# png("resfinderFG_heatmap.png", width = 10, height = 8, units = "in", res = 1200)

# Uncomment the next line to export the plot as pdf
# pdf("resfinderFG_heatmap.pdf", width = 10, height = 8)

# Heatmap
draw(Heatmap(heatmap_matrix_log_T,
             name = "ARG Load (Log10 RPKM)",
             row_title = "Antibiotic Classes",
             column_title = "Group",
             # Formatting parameters
             row_names_gp = gpar(fontsize = 10), 
             row_title_gp = gpar(fontsize = 15),
             column_title_gp = gpar(fontsize = 15), 
             # Display options
             show_row_names = TRUE,
             show_column_names = FALSE,
             cluster_rows = FALSE,
             cluster_columns = FALSE,
             column_names_side = "top",
             row_names_side = "left",
             show_heatmap_legend = TRUE,
             # Color scale
             col = col_fun,
             # Annotations
             top_annotation = column_anno,
             # Grouping
             column_split = factor(tribes, levels = c("Jahai", "Temiar", "Temuan", "Malay")),
             # Cell borders
             cell_fun = function(j, i, x, y, width, height, fill) {
               grid.rect(x = x, y = y, width = width, height = height, 
                         gp = gpar(col = "black", lwd = 0.5, fill = NA))
             }),
     # Add custom legends
     annotation_legend_list = list(tribes_legend),
     merge_legends = TRUE
)

# Uncomment the next line if exporting
# dev.off()

# Set current directory to source file location
current_working_dir <- dirname(rstudioapi::getActiveDocumentContext()$path) 
setwd(current_working_dir)

# Load necessary libraries
library(readxl)
library(tidyverse) 
library(ComplexHeatmap)
library(circlize)
library(reshape2)  
library(RColorBrewer)  
library(viridis)

# Data Preparation and Tranformation ----

# import data
kma_out <- read.csv("kma_out.csv", check.names = FALSE)

# Prepare gene that involve in heatmap
refSequence_meta <- read_excel("unique_refSeq.xlsx", sheet = "Sheet2")

gene_of_interest <- refSequence_meta %>% 
  filter(includeHeatmap == "yes") %>% 
  pull(refSequence)

kma_filtered <- kma_out %>% 
  filter(refSequence %in% gene_of_interest)

# Join filtered_data with refSequence_meta to add only the geneFam column
kma_filtered <- kma_filtered %>%
  left_join(refSequence_meta %>% select(refSequence, geneFam), by = c("refSequence" = "refSequence"))

kma_filtered <- kma_filtered %>%
  relocate(geneFam, .after = "refSequence")

head(kma_filtered) %>% 
  knitr::kable()

# Aggregating RPKM values
geneFam_data <- kma_filtered %>%
  group_by(sampleID, geneFam) %>%
  summarize(geneFam_RPKM = sum(RPKM, na.rm = TRUE)) %>%
  ungroup()

# Define all possible gene families (including missing ones)
all_geneFams <- c(
  "aac(3)", "aac(6')", "aac(6')-aph(2'')", "aph(2'')", "aph(3')", "aph(3'')", "aph(6)",
  "blaACT", "blaCTX-M", "blaOXA", "blaSHV", "blaTEM",
  "mcr", "erm",
  "qnrB", "qnrS", 
  "sul1", "sul2", "sul3",
  "tet"
) 

# Create a complete data frame with all gene families
complete_data <- geneFam_data %>%
  complete(sampleID, geneFam = all_geneFams, fill = list(geneFam_RPKM = 0))

# Spread data into a wide format
heatmap_data <- complete_data %>%
  pivot_wider(names_from = geneFam, values_from = geneFam_RPKM, values_fill = 0)

# Matrix Preparation ----

# Convert to matrix format and set row names
heatmap_matrix <- as.matrix(heatmap_data %>% select(-sampleID))
rownames(heatmap_matrix) <- heatmap_data$sampleID

# Log-transform data
heatmap_matrix_log <- log10(heatmap_matrix + 1)  # Adding 1 to avoid log10(0)
rownames(heatmap_matrix_log) <- heatmap_data$sampleID

# Transpose matrix and reorder rows
heatmap_matrix_log_T <- t(heatmap_matrix_log)
heatmap_matrix_log_T <- heatmap_matrix_log_T[all_geneFams, ]

# Define sample order by tribe
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

# Annotation and Color ----

# Define tribe groups
tribes <- c(
  rep("Jahai", 12),
  rep("Temiar", 9), 
  rep("Temuan", 12), 
  rep("Malay", 9) 
)

# Define drug classes for each gene family
geneFam_category <- c(
  rep("Aminoglycosides", 7),
  rep("Beta-lactamase", 5),
  "Colistin", "Erythromycin",
  rep("Quinolone", 2),
  rep("Sulfonamide", 3),
  "Tetracycline"
)

# Create color mapping for heatmap
breaks <- c(0, 0.3169162, 1, 2, 3, 4, 5, 6)
color_gradient <- colorRampPalette(c("#ffcbd1", "lightpink", "#ee6b6e", "#ff2c2c", "#CC0000", "#990000", "#660000"))(length(breaks) - 1)
colors <- c("white", color_gradient)
col_fun <- colorRamp2(breaks, colors)

# Annotations and Legends ----

# Row annotation for drug classes
row_anno <- rowAnnotation(
  `Drug Class` = anno_simple(geneFam_category, 
                             col = c("Aminoglycosides" = "blue",
                                     "Beta-lactamase" = "khaki3",
                                     "Colistin" = "yellow",
                                     "Erythromycin" = "red",
                                     "Quinolone" = "aquamarine",
                                     "Sulfonamide" = "purple",
                                     "Tetracycline" = "magenta1")),
  show_annotation_name = FALSE
)

# # Create group-wise urbanisation score
# urbanisation_values <- c(
#   rep(0.0, 12), 
#   rep(0.33, 9),   
#   rep(0.67, 12),    
#   rep(1.0, 9)      
# )
# 
# # Column annotation for tribes
# column_anno <- HeatmapAnnotation(
#   `Degree of urbanisation` = urbanisation_values,
#   Tribes = tribes,
#   col = list(`Degree of urbanisation` = colorRamp2(c(0, 0.33, 0.67, 1), 
#                                        c("#00A65A", "#F4E04D", "#FF9F40", "#D80027")
#                                        ),
#              Tribes = c("Jahai" = "darkgreen", 
#                         "Temiar" = "skyblue", 
#                         "Temuan" = "orange", 
#                         "Malay" = "pink")
#              ),
#   show_annotation_name = FALSE,
#   show_legend = FALSE
# )

# One label per tribe
tribe_labels <- c("Jahai", "Temiar", "Temuan", "Malay")

# One fill per tribe
tribe_colors <- c("darkgreen", "skyblue", "orange", "pink")

# Column annotation for tribes
column_anno <- HeatmapAnnotation(
  foo = anno_empty(border = FALSE), # add space for degree of urbanisation
  Tribes = anno_block(
    gp = gpar(fill = tribe_colors, col = NA),
    labels = tribe_labels,
    labels_gp = gpar(col = "white", fontsize = 12, fontface = "bold"),
  ),
  col = list(Tribes = c("Jahai" = "darkgreen",
                        "Temiar" = "skyblue",
                        "Temuan" = "orange",
                        "Malay" = "pink")),
  show_annotation_name = FALSE,
  show_legend = FALSE
)

# Create custom legends
drugClass_legend <- Legend(
  labels = c("Aminoglycosides", "Beta-lactamase", "Colistin", "Erythromycin", 
             "Quinolone", "Sulfonamide", "Tetracycline"),
  legend_gp = gpar(fill = c("blue", "khaki3", "yellow", "red", 
                            "aquamarine", "purple", "magenta1")),
  title = "Drug Class"
)

tribes_legend <- Legend(
  labels = c("Jahai", "Temiar", "Temuan", "Malay"),
  legend_gp = gpar(fill = c("darkgreen", "skyblue", "orange", "pink")),
  title = "Group"
)


# Heatmap and Export ----

# Uncomment the next line to export the plot
png("resfinder_ARGheatmap.png", width = 10, height = 8, units = "in", res = 300)

# Uncomment the next line to export the plot as pdf
pdf("resfinder_ARGheatmap.pdf", width = 10, height = 8)

# Draw the heatmap with legends
draw(Heatmap(heatmap_matrix_log_T,
             name = "ARG Load (Log10 RPKM)",
             row_title = "ARG Classes",
             column_title = "Group",
             # Formatting parameters
             row_names_gp = gpar(fontsize = 12), 
             column_names_gp = gpar(fontsize = 12),
             row_title_gp = gpar(fontsize = 15),
             column_title_gp = gpar(fontsize = 15), 
             # Display options
             show_row_names = TRUE,
             show_column_names = FALSE,
             cluster_rows = FALSE,
             cluster_columns = FALSE,
             column_names_side = "top",
             row_names_side = "left",
             show_heatmap_legend = FALSE,
             # Color scale
             col = col_fun,
             # Annotations
             top_annotation = column_anno,
             left_annotation = row_anno,
             # Grouping
             column_split = factor(tribes, levels = c("Jahai", "Temiar", "Temuan", "Malay")),
             row_split = geneFam_category,
             # Cell borders
             cell_fun = function(j, i, x, y, width, height, fill) {
               grid.rect(x = x, y = y, width = width, height = height, 
                         gp = gpar(col = "black", lwd = 0.5, fill = NA))
             }),
     # Add custom legends
     annotation_legend_list = list(drugClass_legend),
     merge_legends = TRUE, 
     heatmap_legend_side = "bottom"
)

# Uncomment the next line if exporting
dev.off()


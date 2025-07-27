# Set current directory to source file location
current_working_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(current_working_dir)

# Load necessary libraries
library(tidyverse) # data processing
library(vegan) # distance
library(ape) # pcoa
library(ggplot2)

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

# re-shape data
kma_wide <- kma_out %>%
  select(sampleID, refSequence, RPKM) %>% 
  pivot_wider(names_from = refSequence, values_from = RPKM, values_fill = 0) %>% 
  as.data.frame()

# use ID as rownames
rownames(kma_wide) <- kma_wide$sampleID

# drop ID column and save as matrix
kma_mat <- kma_wide %>% 
  select(everything(), -sampleID) %>% 
  as.matrix()

# log transorm RPKM
kma_logRPKM <- log10(kma_mat + 1)

# calculate distance (https://rdrr.io/cran/vegan/man/vegdist.html)
kma_dist_mat <- vegdist(kma_logRPKM)

# pcoa (https://rdrr.io/cran/ape/man/pcoa.html)
res <- pcoa(kma_dist_mat)
res$values
biplot(res)

# ggplot
# obtain coordination
coords <- res$vectors

# Extract first two axes
pcoa_df <- as.data.frame(coords[, 1:2])
pcoa_df$sampleID <- rownames(pcoa_df)

# Merge with metadata
meta <- kma_out %>% select(sampleID, tribe, Group.Bi) %>% distinct()
pcoa_merged <- pcoa_df %>% 
  left_join(meta, by = "sampleID") %>% 
  select(sampleID, tribe, Group.Bi, Axis.1, Axis.2)

# plot
PCoA <- ggplot(pcoa_merged, aes(x = Axis.1, y = Axis.2, color = tribe, shape = Group.Bi)) +
  geom_point(size = 3) +
  stat_ellipse(type = "t", level = 0.95, 
               linetype = 1, linewidth = 0.5) +
  theme_minimal() +
  scale_color_manual(values = c("Jahai" = "darkgreen", "Temiar" = "skyblue",
                               "Temuan" = "orange", "Malay" = "pink")) +
  labs(x = "PCoA1", y = "PCoA2",
       color = "Groups", shape = "Region Groups") +
  coord_equal()
PCoA

# Save as PNG
# ggsave("resfinder_rpkmPlot_PCoA.png", plot = PCoA, width = 10, height = 6, units = "in")

# Save as PDF
# ggsave("resfinder_rpkmPlot_PCoA.pdf", plot = PCoA, width = 10, height = 6, units = "in")


##################################

# Claude

# Set current directory to source file location
current_working_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(current_working_dir)

# Load necessary libraries
library(tidyverse) # data processing
library(vegan) # distance
library(ape) # pcoa
library(ggplot2)

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

# re-shape data
kma_wide <- kma_out %>%
  select(sampleID, refSequence, RPKM) %>% 
  pivot_wider(names_from = refSequence, values_from = RPKM, values_fill = 0) %>% 
  as.data.frame()

# use ID as rownames
rownames(kma_wide) <- kma_wide$sampleID

# drop ID column and save as matrix
kma_mat <- kma_wide %>% 
  select(everything(), -sampleID) %>% 
  as.matrix()

# log transform RPKM
kma_logRPKM <- log10(kma_mat + 1)

# calculate distance
kma_dist_mat <- vegdist(kma_logRPKM)

# pcoa
res <- pcoa(kma_dist_mat)

# obtain coordination
coords <- res$vectors
pcoa_df <- as.data.frame(coords[, 1:2])
pcoa_df$sampleID <- rownames(pcoa_df)

# Merge with metadata
meta <- kma_out %>% select(sampleID, tribe, Group.Bi) %>% distinct()
pcoa_merged <- pcoa_df %>% 
  left_join(meta, by = "sampleID") %>% 
  select(sampleID, tribe, Group.Bi, Axis.1, Axis.2)


# Calculate the direction of maximum variation for each group
group_directions <- pcoa_merged %>% 
  group_by(tribe) %>% 
  summarise(
    mean_x = mean(Axis.1),
    mean_y = mean(Axis.2),
    var_x = var(Axis.1),
    var_y = var(Axis.2),
    .groups = 'drop'
  ) %>% 
  mutate(
    # Direction vector (simplified)
    dir_x = mean_x + sign(mean_x) * sqrt(var_x) * 0.5,
    dir_y = mean_y + sign(mean_y) * sqrt(var_y) * 0.5
  )

p3 <- ggplot(pcoa_merged, aes(x = Axis.1, y = Axis.2, color = tribe, shape = Group.Bi)) +
  geom_point(size = 3) +
  stat_ellipse(type = "t", level = 0.95, linetype = 1, linewidth = 1) +
  # Add directional arrows
  geom_segment(data = group_directions,
               aes(x = mean_x, y = mean_y, xend = dir_x, yend = dir_y, color = tribe),
               arrow = arrow(length = unit(0.3, "cm")),
               linewidth = 1, inherit.aes = FALSE) +
  theme_minimal() +
  scale_color_manual(values = c("J" = "darkgreen", "Tr" = "skyblue",
                                "Tn" = "orange", "My" = "pink")) +
  labs(x = "PCoA1", y = "PCoA2",
       color = "Groups", shape = "Region Groups",
       title = "PCoA with Group Variation Vectors") +
  coord_equal()

print(p3)

# Display eigenvalues for context
cat("PCoA Eigenvalues (proportion of variance explained):\n")
print(res$values$Relative_eig[1:5])

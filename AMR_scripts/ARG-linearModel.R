#Linear regression for differential abundance analysis of ARG in 
##Urban vs Rural groups

# Set current directory to source file location
current_working_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(current_working_dir)

# Load necessary libraries
library(ggplot2)
library(tidyverse)
library(compositions)

# import data
kma_out <- read.csv("kma_out.csv", check.names = FALSE)

# import meta
kma_meta <- read.csv("Metafile.csv", check.names = FALSE)

# add info on region grouping (1 – Rural; 0 – Urban)
kma_out <- kma_out %>%
  left_join(
    kma_meta %>% dplyr::select(Group, Group.Bi) %>% unique(),
    by = c("tribe" = "Group")
  )
kma_out <- kma_out %>% 
  mutate(Group.Bi = factor(Group.Bi, levels = c(0, 1), labels = c("Urban", "Rural"))) %>% 
  relocate(Group.Bi, .after = tribe)

# Create a complete data frame
complete_data <- kma_out %>%
  pivot_wider(id_cols = c(sampleID, Group.Bi),
              names_from = refSequence, 
              values_from = RPKM, 
              values_fill = 0)

# extract only counts
arg_matrix <- complete_data %>%
  select(-sampleID, -Group.Bi) %>%
  as.data.frame()

# add pseudocount
arg_matrix_pseudo <- arg_matrix + 1

# clr transformation 
clr_matrix <- clr(acomp(arg_matrix_pseudo)) # acomp to inform compositional data

# re-create dataframe
clr_df <- as.data.frame(clr_matrix)
clr_df$sampleID <- complete_data$sampleID
clr_df$Group.Bi <- complete_data$Group.Bi
clr_df <- clr_df %>% 
  select(sampleID, Group.Bi, everything())

# Function to perform linear model for each ARG
perform_lm_daa <- function(clr_data) {
  # Get ARG names (exclude metadata columns)
  arg_names <- setdiff(names(clr_data), c("Group.Bi", "sampleID"))
  
  results <- data.frame(
    ARG = character(),
    estimate = numeric(),
    std_error = numeric(),
    t_value = numeric(),
    p_value = numeric(),
    adj_p_value = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (arg in arg_names) {
    # Fit linear model
    formula_str <- paste0("`", arg, "` ~ Group.Bi")
    model <- lm(as.formula(formula_str), data = clr_data)
    
    # Extract results
    coef_summary <- summary(model)$coefficients
    
    if (nrow(coef_summary) > 1) {  # Check if Group.Bi coefficient exists
      results <- rbind(results, data.frame(
        ARG = arg,
        estimate = coef_summary[2, "Estimate"],
        std_error = coef_summary[2, "Std. Error"],
        t_value = coef_summary[2, "t value"],
        p_value = coef_summary[2, "Pr(>|t|)"]
      ))
    }
  }
  
  # Multiple testing correction
  results$adj_p_value <- p.adjust(results$p_value, method = "BH")
  
  return(results)
}

# Run the analysis
daa_results <- perform_lm_daa(clr_df)

# export
# write.csv(daa_results, "ARGs_DAA_results.csv", row.names = FALSE)

# save RData
# save.image(file = "ARGs_DAA.RData")

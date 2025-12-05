# R Script for Plotting RNA Statistical Potentials (Score vs. Distance)
# Requires the 'tidyverse' package (which includes ggplot2 and dplyr)

# ----------------------------------------------------------------------
# 1. CONFIGURATION & USER INPUT
# ----------------------------------------------------------------------

# Load the required plotting library
# Ensure tidyverse is installed (install.packages("tidyverse") if needed)
library(tidyverse)

# Prompt the user for the directory path
cat("Please enter the full directory path where the 10 potential TXT files (e.g., AG.txt) are located:\n")
POTENTIAL_FILES_PATH <- readline(prompt = "> ")

# Remove potential trailing slashes for consistency
POTENTIAL_FILES_PATH <- normalizePath(POTENTIAL_FILES_PATH, mustWork = FALSE)

# Define constants
BIN_WIDTH <- 1.0 # Angstrom, matches Python script
MAX_DIST <- 20.0 # Angstrom

# ----------------------------------------------------------------------
# 2. LOAD DATA
# ----------------------------------------------------------------------

# Get a list of all potential files (e.g., "AG.txt", "CU.txt")
file_list <- list.files(path = POTENTIAL_FILES_PATH, pattern = "\\.txt$", full.names = TRUE)

if (length(file_list) == 0) {
  stop(paste("ERROR: No .txt potential files found in the specified path:", POTENTIAL_FILES_PATH))
}

cat(paste("Found", length(file_list), "potential files. Combining data...\n"))

# Read and combine all data files
all_data <- map_dfr(file_list, function(file_path) {
  
  # Extract the residue pair name from the filename (e.g., "AG" from ".../AG.txt")
  filename <- basename(file_path)
  pair_name <- gsub("\\.txt$", "", filename)
  
  # Read the single column of scores
  scores <- read_delim(
    file_path, 
    col_names = "Score", 
    delim = "\n", 
    col_types = cols(Score = col_double()),
    show_col_types = FALSE
  )
  
  # Add metadata: Pair name and Distance
  scores %>%
    mutate(
      Pair = factor(pair_name), # Keep pairs as factors for easy plotting
      # Calculate the distance center for the bin (0.5, 1.5, 2.5, ...)
      Distance = (row_number() - 0.5) * BIN_WIDTH
    )
})

# ----------------------------------------------------------------------
# 3. GENERATE PLOT
# ----------------------------------------------------------------------

potential_plot <- all_data %>%
  ggplot(aes(x = Distance, y = Score, color = Pair)) +
  
  # Plotting the discrete points
  geom_point(size = 2, alpha = 0.8) +
  
  # Connecting the points with lines
  geom_line(linewidth = 1) +
  
  # Use facet_wrap to separate the profiles for each pair
  facet_wrap(~ Pair, scales = "free_y", ncol = 4) +
  
  # Horizontal line at 0 (the reference energy level)
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", alpha = 0.6) +
  
  # Customize Axis and Labels
  scale_x_continuous(
    name = "Distance (Å)", 
    breaks = seq(0, MAX_DIST, 5), 
    limits = c(0, MAX_DIST)
  ) +
  labs(
    title = "RNA Statistical Interaction Potentials (C3' - C3')",
    y = "Pseudo-Energy Score (-ln(f_obs / f_ref))",
    caption = paste("Data loaded from:", POTENTIAL_FILES_PATH)
  ) +
  
  # Apply a clean theme for aesthetics
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none", # Hide redundant legend
    plot.title = element_text(face = "bold", hjust = 0.5),
    strip.background = element_rect(fill = "#e0e0e0", color = NA),
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

# Print the plot
print(potential_plot)
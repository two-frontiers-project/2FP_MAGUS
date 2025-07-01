#!/usr/bin/env Rscript

#### example script for how you can parse output of find-euks to get a list of potential euk bin filenames

suppressMessages(library(tidyverse))

# Parse arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
  stop("Usage: Rscript filter_euks.R <input_csv> <output_file> <min_bin_size,contamination_cutoff>")
}

input_csv <- args[1]
output_file <- args[2]
thresholds <- str_split(args[3], ",")[[1]]

if (length(thresholds) != 2) {
  stop("Thresholds must be two comma-separated values: <min_bin_size,contamination_cutoff>")
}

min_bin_size <- as.numeric(thresholds[1])
contam_cutoff <- as.numeric(thresholds[2])

# Read input file
df <- read_csv(input_csv, show_col_types = FALSE)

# Ensure numeric columns
df <- df %>%
  mutate(
    bin_size = as.numeric(bin_size),
    checkm_Contamination = as.numeric(checkm_Contamination)
  )

# Filter
filtered_bins <- df %>%
  filter(
    bin_size > min_bin_size,
    checkm_Contamination > contam_cutoff
  ) %>%
  distinct(bin_name)

# Write output
write_lines(filtered_bins$bin_name, output_file)

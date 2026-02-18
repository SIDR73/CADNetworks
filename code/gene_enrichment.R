# CAD Interactor Enrichment Analysis 
# Tests SMC, OMIM CAD, and Aragam GWAS enrichment using BioPlex Interactions vs. random sampling permutations

library(tidyverse)

dir.create("./figures/enrichment_combined", recursive = TRUE, showWarnings = FALSE)
set.seed(42)

# load data
network <- read_csv("./data/master_network_full_annotated.csv", show_col_types = FALSE)
bioplex <- read_tsv("./data/BioPlex_293T_Network_10K_Dec_2019.tsv", show_col_types = FALSE)

cat(sprintf("Experimental interactions: %d\n", nrow(network)))
cat(sprintf("BioPlex interactions: %d\n\n", nrow(bioplex)))

# load reference gene lists
# SMC genes
smc_genes <- read_tsv("./data/smc_expressed_other_enriched_HPA_tissue_category_rna_Any_Tissue.tsv",
                      show_col_types = FALSE)
smc_gene_list <- unique(smc_genes$Gene)

# OMIM CAD genes
omim_data <- read_tsv("./data/cad_gene_filtering/OMIM-Gene-Map-Search_coronary_artery_disease.tsv",
                      skip = 4, show_col_types = FALSE)
omim_cad_genes <- unique(omim_data$`Approved Symbol`)
omim_cad_genes <- omim_cad_genes[!is.na(omim_cad_genes)]
cat(sprintf("OMIM CAD genes: %d\n", length(omim_cad_genes)))

# Aragam 2022 GWAS genes
aragam_data <- read_csv("./data/Aragam_2022_S31_CAD_GWAS.csv", show_col_types = FALSE)
aragam_cad_genes <- unique(aragam_data$most_likely_causal_gene)
aragam_cad_genes <- aragam_cad_genes[!is.na(aragam_cad_genes)]
cat(sprintf("  Aragam GWAS genes: %d\n\n", length(aragam_cad_genes)))

# store lists
reference_gene_sets <- list(
  SMC    = smc_gene_list,
  OMIM   = omim_cad_genes,
  Aragam = aragam_cad_genes
)

# Bioplex interaction map
bioplex_map <- rbind(
  bioplex %>% select(gene = SymbolA, interactor = SymbolB),
  bioplex %>% select(gene = SymbolB, interactor = SymbolA)
) %>%
  filter(
    gene != interactor,
    !is.na(gene),
    !is.na(interactor)
  ) %>%
  distinct()

all_genes_in_bioplex <- unique(c(bioplex_map$gene, bioplex_map$interactor))

cat(sprintf("BioPlex interactions - bidirectional: %d\n", nrow(bioplex_map)))
cat(sprintf("Unique genes in BioPlex: %d\n\n", length(all_genes_in_bioplex)))

# 1st Degree Interactions
network_with_genes <- network %>%
  mutate(Gene = str_remove(source_harmonized, "-.*$"))

cad_interactions <- network_with_genes %>%
  select(Gene, Isoform = source_harmonized, Interactor = target_harmonized,
         Specificity = specificity) %>%
  distinct()

first_degree_all <- unique(cad_interactions$Interactor)
cat(sprintf("1st degree interactors: %d\n", length(first_degree_all)))

# 2nd Degree Interactors
second_degree_raw <- bioplex_map %>%
  filter(gene %in% first_degree_all) %>%
  select(FirstDegree = gene, SecondDegree = interactor) %>%
  distinct()

second_degree_all <- unique(second_degree_raw$SecondDegree)

cat(sprintf("2nd degree candidates: %d\n", nrow(second_degree_raw)))
cat(sprintf("2nd degree unique genes: %d\n\n", length(second_degree_all)))

# Find total network size
total_network_size <- length(unique(c(first_degree_all, second_degree_all)))
cat(sprintf("Total network size (1st + 2nd degree): %d\n\n", total_network_size))

# Filter interactor sets for each reference set
# Create configurations for all combinations
all_configs <- list()

for (ref_name in names(reference_gene_sets)) {
  ref_genes <- reference_gene_sets[[ref_name]]
  
  first_filtered  <- intersect(first_degree_all, ref_genes)
  second_filtered <- intersect(second_degree_all, ref_genes)
  combined        <- unique(c(first_filtered, second_filtered))
  
  all_configs[[paste0("1st_degree_", ref_name)]]       <- first_filtered
  all_configs[[paste0("2nd_degree_", ref_name)]]       <- second_filtered
  all_configs[[paste0("1st_2nd_combined_", ref_name)]] <- combined
  
  cat(sprintf("%s:\n", ref_name))
  cat(sprintf("  1st degree: %d genes\n", length(first_filtered)))
  cat(sprintf("  2nd degree: %d genes\n", length(second_filtered)))
  cat(sprintf("  Combined: %d genes\n\n", length(combined)))
}

# Permutation testing
run_permutation_test <- function(observed_genes,
                                 universe_genes,
                                 target_gene_list,
                                 total_network_size,
                                 n_permutations = 10000) {
  
  n_genes      <- length(observed_genes)
  obs_count    <- sum(observed_genes %in% target_gene_list)
  obs_fraction <- obs_count / total_network_size  # Fraction of total network
  
  null_counts <- replicate(n_permutations, {
    random_genes <- sample(universe_genes, n_genes, replace = FALSE)
    sum(random_genes %in% target_gene_list)
  })
  
  expected <- mean(null_counts)
  p_value  <- mean(null_counts >= obs_count)
  
  list(
    n_genes         = n_genes,
    obs_count       = obs_count,
    obs_fraction    = obs_fraction,
    expected_count  = expected,
    expected_frac   = expected / total_network_size,  # Fraction of TOTAL network
    sd              = sd(null_counts),
    fold_enrichment = obs_count / expected,
    p_value         = p_value,
    null_counts     = null_counts,
    total_network_size = total_network_size
  )
}

# Permutation tests for each configuration
n_permutations <- 10000
all_results    <- list()
perm_data      <- list()

for (config_name in names(all_configs)) {
  
  gene_set <- all_configs[[config_name]]
  
  # Determine which reference set to use
  ref_name <- case_when(
    grepl("_SMC$", config_name)    ~ "SMC",
    grepl("_OMIM$", config_name)   ~ "OMIM",
    grepl("_Aragam$", config_name) ~ "Aragam",
    TRUE ~ NA_character_
  )
  
  target_genes <- reference_gene_sets[[ref_name]]
  
  result <- run_permutation_test(
    observed_genes     = gene_set,
    universe_genes     = all_genes_in_bioplex,
    target_gene_list   = target_genes,
    total_network_size = total_network_size,
    n_permutations     = n_permutations
  )
  
  cat(sprintf("%s genes observed: %d / %d (%.1f%% of total network)\n",
              ref_name, result$obs_count, total_network_size, result$obs_fraction * 100))
  cat(sprintf("Expected: %.2f / %d (%.1f%% of total network)\n",
              result$expected_count, total_network_size, result$expected_frac * 100))
  cat(sprintf("Fold enrichment: %.2fx\n", result$fold_enrichment))
  cat(sprintf("Z-score: %.2f\n",
              (result$obs_count - result$expected_count) / result$sd))
  cat(sprintf("Permutation p-value: %.4f %s\n\n",
              result$p_value,
              ifelse(result$p_value < 0.05, "***", "")))
  
  all_results[[config_name]] <- data.frame(
    Configuration       = config_name,
    Reference_Set       = ref_name,
    Total_Network_Size  = total_network_size,
    N_Genes_Tested      = result$n_genes,
    Obs_Count           = result$obs_count,
    Obs_Fraction_Pct    = result$obs_fraction * 100,
    Expected_Count      = result$expected_count,
    Expected_Frac_Pct   = result$expected_frac * 100,
    SD                  = result$sd,
    Z_Score             = (result$obs_count - result$expected_count) / result$sd,
    Fold_Enrichment     = result$fold_enrichment,
    P_Value             = result$p_value,
    Significant         = result$p_value < 0.05
  )
  
  perm_data[[config_name]] <- result
}

# Save results

results_df <- bind_rows(all_results)

write.csv(results_df,
          "./figures/enrichment_combined/enrichment_permutation_results_all.csv",
          row.names = FALSE)

cat("Results saved to: enrichment_permutation_results_all.csv\n\n")

# Plot with each configuration with broken axes

for (config_name in names(perm_data)) {
  
  result   <- perm_data[[config_name]]
  ref_name <- results_df$Reference_Set[results_df$Configuration == config_name]
  
  filename <- sprintf("./figures/enrichment_combined/permutation_%s.pdf", config_name)
  pdf(filename, width = 10, height = 8)
  
  # Convert counts to fractions
  null_fractions <- result$null_counts / result$total_network_size
  obs_fraction   <- result$obs_count / result$total_network_size
  exp_fraction   <- result$expected_count / result$total_network_size
  
  h <- hist(null_fractions, breaks = 50, plot = FALSE)
  
  y_max <- max(h$counts) * 1.25
  
  # broken axis plot
  null_max <- max(null_fractions)
  null_width <- null_max - min(h$mids)
  gap_width <- null_width * 0.1
  obs_width <- obs_fraction * 0.05
  
  plot(NA, xlim = c(0, null_width + gap_width + obs_width * 2), 
       ylim = c(0, y_max),
       xlab = sprintf("Fraction of total network (%s genes)", ref_name),
       ylab = "Frequency",
       main = sprintf("%s\nObserved: %.1f%% | Expected: %.1f%% ± %.2f%%",
                      config_name, 
                      obs_fraction * 100,
                      exp_fraction * 100,
                      (result$sd / result$total_network_size) * 100,
                      result$fold_enrichment),
       xaxt = "n",
       cex.main = 0.9)
  
  # plot null distribution
  null_x_transformed <- h$mids - min(h$mids)
  lines(null_x_transformed, h$counts, lwd = 2, col = "gray60")
  
  # plot expected line
  exp_x_transformed <- exp_fraction - min(h$mids)
  abline(v = exp_x_transformed, col = "#1f77b4", lwd = 2, lty = 2)
  
  # observed line
  obs_x_transformed <- null_width + gap_width + obs_width
  abline(v = obs_x_transformed, col = "black", lwd = 3)
  
  # break symbols
  break_x <- null_width + gap_width / 2
  segments(break_x - gap_width * 0.15, -y_max * 0.02, 
           break_x - gap_width * 0.05, y_max * 0.02, lwd = 2)
  segments(break_x + gap_width * 0.05, -y_max * 0.02, 
           break_x + gap_width * 0.15, y_max * 0.02, lwd = 2)
  
  # Draw axis lines with breaks
  # Left section
  segments(0, 0, null_width, 0, lwd = 1)
  # Right section  
  segments(null_width + gap_width, 0, 
           null_width + gap_width + obs_width * 2, 0, lwd = 1)
  
  # Custom x-axis ticks
  # Custom x-axis ticks with adaptive precision
  left_ticks <- pretty(c(min(h$mids), null_max), n = 5)
  left_ticks <- left_ticks[left_ticks >= min(h$mids) & left_ticks <= null_max]
  
  # Use 4 decimal places for fractions
  axis(1, at = left_ticks - min(h$mids), 
       labels = sprintf("%.4f", left_ticks))
  
  axis(1, at = obs_x_transformed, 
       labels = sprintf("%.4f", obs_fraction))
  
  
  # Add text annotations
  text(exp_x_transformed, y_max * 0.85,
       sprintf("Expected = %.1f%%",
               exp_fraction * 100),
       pos = 4, col = "#1f77b4", font = 2, cex = 1.0)
  
  text(obs_x_transformed, y_max * 0.85,
       sprintf("Observed = %.1f%%\nFold = %.2fx\nP = %.4f",
               obs_fraction * 100, 
               result$fold_enrichment,
               result$p_value),
       pos = 2, col = "black", font = 2, cex = 1.0)
  
  dev.off()
  cat(sprintf("  Saved: %s\n", basename(filename)))
}

# Summary Bar Plots
for (ref_name in names(reference_gene_sets)) {
  
  subset_data <- results_df %>%
    filter(Reference_Set == ref_name) %>%
    arrange(Configuration)
  
  if (nrow(subset_data) == 0) next
  
  pdf(sprintf("./figures/enrichment_combined/summary_%s.pdf", ref_name), 
      width = 9, height = 6)
  
  # Create shorter configuration labels
  short_labels <- gsub("1st_2nd_combined_", "Combined_", subset_data$Configuration)
  short_labels <- gsub("1st_degree_", "1st_", short_labels)
  short_labels <- gsub("2nd_degree_", "2nd_", short_labels)
  short_labels <- gsub(paste0("_", ref_name), "", short_labels)
  
  bar_vals <- matrix(
    c(subset_data$Obs_Fraction_Pct, subset_data$Expected_Frac_Pct),
    nrow = 2, byrow = TRUE,
    dimnames = list(c("Observed", "Expected"), short_labels)
  )
  
  bp <- barplot(bar_vals,
                beside = TRUE,
                col = c("#d62728", "#1f77b4"),
                ylim = c(0, max(bar_vals) * 1.3),
                ylab = sprintf("%% of total network (%s genes)", ref_name),
                main = sprintf("%s Enrichment vs Random Sampling", ref_name),
                las = 1,
                cex.names = 0.85,
                legend.text = rownames(bar_vals),
                args.legend = list(x = "topright", bty = "n", cex = 0.9))
  
  # significance annotations
  for (i in seq_along(subset_data$Configuration)) {
    if (subset_data$Significant[i]) {
      text(mean(bp[, i]), max(bar_vals[, i]) * 1.08, "***", cex = 1.4)
    }
  }
  
  dev.off()
}
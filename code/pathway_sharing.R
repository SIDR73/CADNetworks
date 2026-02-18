# Pathway Sharing Analysis

library(tidyverse)
library(msigdbr)
library(pheatmap)
library(RColorBrewer)

dir.create("./figures/pathway_sharing", recursive = TRUE, showWarnings = FALSE)
set.seed(42)

network <- read_csv("./data/master_network_full_annotated.csv", show_col_types = FALSE)
biogrid <- read_tsv("./data/biogrid.txt", col_types = cols(.default = "c"))

cat(sprintf("Experimental interactions: %d\n", nrow(network)))
cat(sprintf("BioGRID interactions: %d\n\n", nrow(biogrid)))

# Biogrid Map
biogrid_physical <- biogrid %>%
  select(
    gene_a = `Official Symbol Interactor A`,
    gene_b = `Official Symbol Interactor B`,
    system_type = `Experimental System Type`
  ) %>%
  filter(
    system_type == "physical",
    gene_a != gene_b,
    !is.na(gene_a),
    !is.na(gene_b)
  )

biogrid_map <- rbind(
  biogrid_physical %>% select(gene = gene_a, interactor = gene_b),
  biogrid_physical %>% select(gene = gene_b, interactor = gene_a)
) %>%
  distinct()

cat(sprintf("BioGRID interactions: %d\n\n", nrow(biogrid_map)))

# Interactor Sets
network_with_genes <- network %>%
  mutate(Gene = str_remove(source_harmonized, "-.*$"))

# 1st Degree
first_degree <- network_with_genes %>%
  filter(specificity %in% c("risk-specific", "nonrisk-specific")) %>%
  select(Gene, Isoform = source_harmonized, Interactor = target_harmonized, Specificity = specificity) %>%
  distinct()

cat(sprintf("1st degree (experimental): %d pairs\n", nrow(first_degree)))
cat(sprintf("Risk: %d\n", sum(first_degree$Specificity == "risk-specific")))
cat(sprintf("Non-risk: %d\n", sum(first_degree$Specificity == "nonrisk-specific")))

# 2nd Degree
risk_1st_interactors <- first_degree %>%
  filter(Specificity == "risk-specific") %>%
  pull(Interactor) %>%
  unique()

risk_2nd <- biogrid_map %>%
  filter(gene %in% risk_1st_interactors) %>%
  select(FirstDegree = gene, SecondDegree = interactor) %>%
  distinct()

risk_2nd_with_info <- risk_2nd %>%
  left_join(
    first_degree %>% 
      filter(Specificity == "risk-specific") %>%
      select(Gene, Isoform, FirstDegree = Interactor),
    by = "FirstDegree",
    relationship = "many-to-many"
  ) %>%
  select(Gene, Isoform, Interactor = SecondDegree) %>%
  mutate(Specificity = "risk-specific") %>%
  distinct()

nonrisk_1st_interactors <- first_degree %>%
  filter(Specificity == "nonrisk-specific") %>%
  pull(Interactor) %>%
  unique()

nonrisk_2nd <- biogrid_map %>%
  filter(gene %in% nonrisk_1st_interactors) %>%
  select(FirstDegree = gene, SecondDegree = interactor) %>%
  distinct()

nonrisk_2nd_with_info <- nonrisk_2nd %>%
  left_join(
    first_degree %>% 
      filter(Specificity == "nonrisk-specific") %>%
      select(Gene, Isoform, FirstDegree = Interactor),
    by = "FirstDegree",
    relationship = "many-to-many"
  ) %>%
  select(Gene, Isoform, Interactor = SecondDegree) %>%
  mutate(Specificity = "nonrisk-specific") %>%
  distinct()

second_degree <- rbind(risk_2nd_with_info, nonrisk_2nd_with_info)

cat(sprintf("2nd degree (BioGRID): %d pairs\n", nrow(second_degree)))
cat(sprintf("Risk: %d\n", sum(second_degree$Specificity == "risk-specific")))
cat(sprintf("Non-risk: %d\n\n", sum(second_degree$Specificity == "nonrisk-specific")))

# Load pathway datasets from MSIGDB
msigdb_h <- msigdbr(species = "Homo sapiens", collection = "H")
msigdb_c2 <- msigdbr(species = "Homo sapiens", collection = "C2")
msigdb_c5 <- msigdbr(species = "Homo sapiens", collection = "C5")

pathway_databases <- list(
  Reactome = msigdb_c2 %>% filter(gs_subcollection == "CP:REACTOME"),
  WikiPathways = msigdb_c2 %>% filter(gs_subcollection == "CP:WIKIPATHWAYS"),
  GO_BP = msigdb_c5 %>% filter(gs_subcollection == "GO:BP"),
  GO_MF = msigdb_c5 %>% filter(gs_subcollection == "GO:MF"),
  GO_CC = msigdb_c5 %>% filter(gs_subcollection == "GO:CC")
)

pathway_databases <- pathway_databases[sapply(pathway_databases, nrow) > 0]

# Find the top N pathways for a given gene
get_top_pathways <- function(gene, db_data, n_top = 1) {
  pathways <- db_data %>%
    filter(gene_symbol == gene) %>%
    pull(gs_name) %>%
    unique()
  
  if (length(pathways) == 0) return(character(0))
  head(pathways, n_top)
}

# Check pathway sharing
check_sharing <- function(interactor_data, db_data, n_top = 1) {
  results <- list()
  
  for (i in 1:nrow(interactor_data)) {
    gene <- interactor_data$Gene[i]
    interactor <- interactor_data$Interactor[i]
    
    gene_pathways <- get_top_pathways(gene, db_data, n_top)
    interactor_pathways <- get_top_pathways(interactor, db_data, n_top)
    
    shares <- FALSE
    if (length(gene_pathways) > 0 && length(interactor_pathways) > 0) {
      shared <- intersect(gene_pathways, interactor_pathways)
      shares <- length(shared) > 0
    }
    
    results[[i]] <- data.frame(
      Gene = gene,
      Isoform = interactor_data$Isoform[i],
      Interactor = interactor,
      Specificity = interactor_data$Specificity[i],
      Shares = shares
    )
  }
  
  bind_rows(results)
}

# Global Analysis: Risk vs. Non-risk sharing 
n_top_list <- c(1, 3, 5, 10)
degree_configs <- list(
  "1st_only" = first_degree,
  "1st_2nd_combined" = rbind(first_degree, second_degree)
)

all_results <- list()

for (degree_name in names(degree_configs)) {
  interactor_data <- degree_configs[[degree_name]]
  
  for (n_top in n_top_list) {
    cat(sprintf("Testing TOP %d pathway(s)...\n", n_top))
    
    for (db_name in names(pathway_databases)) {
      db_data <- pathway_databases[[db_name]]
      
      sharing_results <- check_sharing(interactor_data, db_data, n_top)
      
      contingency <- sharing_results %>%
        group_by(Specificity, Shares) %>%
        summarise(n = n(), .groups = "drop") %>%
        pivot_wider(names_from = Shares, values_from = n, values_fill = 0)
      
      risk_share <- contingency %>% 
        filter(Specificity == "risk-specific") %>% 
        pull(if("TRUE" %in% names(.)) `TRUE` else 2) %>% 
        {if(length(.) == 0) 0 else .}
      
      risk_no_share <- contingency %>% 
        filter(Specificity == "risk-specific") %>% 
        pull(if("FALSE" %in% names(.)) `FALSE` else 2) %>% 
        {if(length(.) == 0) 0 else .}
      
      nonrisk_share <- contingency %>% 
        filter(Specificity == "nonrisk-specific") %>% 
        pull(if("TRUE" %in% names(.)) `TRUE` else 2) %>% 
        {if(length(.) == 0) 0 else .}
      
      nonrisk_no_share <- contingency %>% 
        filter(Specificity == "nonrisk-specific") %>% 
        pull(if("FALSE" %in% names(.)) `FALSE` else 2) %>% 
        {if(length(.) == 0) 0 else .}
      
      fisher_matrix <- matrix(c(risk_share, risk_no_share, 
                                nonrisk_share, nonrisk_no_share),
                              nrow = 2, byrow = TRUE)
      
      fisher_test <- fisher.test(fisher_matrix)
      
      risk_total <- risk_share + risk_no_share
      nonrisk_total <- nonrisk_share + nonrisk_no_share
      
      risk_pct <- ifelse(risk_total > 0, risk_share / risk_total * 100, 0)
      nonrisk_pct <- ifelse(nonrisk_total > 0, nonrisk_share / nonrisk_total * 100, 0)
      
      all_results[[length(all_results) + 1]] <- data.frame(
        Configuration = degree_name,
        N_Top = n_top,
        Database = db_name,
        Risk_Share = risk_share,
        Risk_NoShare = risk_no_share,
        Risk_Total = risk_total,
        Risk_Pct = risk_pct,
        NonRisk_Share = nonrisk_share,
        NonRisk_NoShare = nonrisk_no_share,
        NonRisk_Total = nonrisk_total,
        NonRisk_Pct = nonrisk_pct,
        P_Value = fisher_test$p.value,
        Odds_Ratio = fisher_test$estimate,
        Significant = fisher_test$p.value < 0.05
      )
      
      if (fisher_test$p.value < 0.05 || n_top == 1) {
        cat(sprintf("%s (TOP %d):\n", db_name, n_top))
        cat(sprintf("Risk: %d/%d (%.1f%%) share pathways\n", 
                    risk_share, risk_total, risk_pct))
        cat(sprintf("Non-risk: %d/%d (%.1f%%) share pathways\n", 
                    nonrisk_share, nonrisk_total, nonrisk_pct))
        cat(sprintf("Fisher's p-value: %.4f %s\n", 
                    fisher_test$p.value,
                    ifelse(fisher_test$p.value < 0.05, "***", "")))
        cat(sprintf("Odds ratio: %.3f\n\n", fisher_test$estimate))
      }
    }
  }
}

results_df <- bind_rows(all_results)
write.csv(results_df, "./figures/pathway_sharing/pathway_sharing_fisher_exact.csv", row.names = FALSE)

# Permutation Analysis
# essentially shuffle the "risk" vs "non-risk" labels randomly

n_permutations <- 10000

sig_results <- results_df %>%
  filter(Significant == TRUE, Odds_Ratio > 1) %>%  # Only where risk > non-risk
  arrange(P_Value)

cat(sprintf("Found %d significant results to analyze\n\n", nrow(sig_results)))

permutation_results <- list()

for (i in 1:nrow(sig_results)) {
  result <- sig_results[i, ]
  
  cat(sprintf("Permutation test %d/%d: %s - %s (TOP %d)\n", 
              i, nrow(sig_results), 
              result$Configuration, result$Database, result$N_Top))
  
  # Get the data for this configuration
  interactor_data <- degree_configs[[result$Configuration]]
  db_data <- pathway_databases[[result$Database]]
  
  # Get sharing results
  sharing_results <- check_sharing(interactor_data, db_data, result$N_Top)
  
  # Observed difference
  obs_risk_pct <- result$Risk_Pct
  obs_nonrisk_pct <- result$NonRisk_Pct
  obs_diff <- obs_risk_pct - obs_nonrisk_pct
  
  # Permutation test
  null_diffs <- numeric(n_permutations)
  
  for (perm in 1:n_permutations) {
    # Shuffle specificity labels
    shuffled <- sharing_results
    shuffled$Specificity <- sample(shuffled$Specificity)
    
    # Recalculate sharing percentages
    perm_contingency <- shuffled %>%
      group_by(Specificity, Shares) %>%
      summarise(n = n(), .groups = "drop") %>%
      pivot_wider(names_from = Shares, values_from = n, values_fill = 0)
    
    perm_risk_share <- perm_contingency %>% 
      filter(Specificity == "risk-specific") %>% 
      pull(if("TRUE" %in% names(.)) `TRUE` else 2) %>% 
      {if(length(.) == 0) 0 else .}
    
    perm_risk_total <- perm_contingency %>%
      filter(Specificity == "risk-specific") %>%
      summarise(total = sum(c_across(where(is.numeric)))) %>%
      pull(total)
    
    perm_nonrisk_share <- perm_contingency %>% 
      filter(Specificity == "nonrisk-specific") %>% 
      pull(if("TRUE" %in% names(.)) `TRUE` else 2) %>% 
      {if(length(.) == 0) 0 else .}
    
    perm_nonrisk_total <- perm_contingency %>%
      filter(Specificity == "nonrisk-specific") %>%
      summarise(total = sum(c_across(where(is.numeric)))) %>%
      pull(total)
    
    perm_risk_pct <- ifelse(perm_risk_total > 0, perm_risk_share / perm_risk_total * 100, 0)
    perm_nonrisk_pct <- ifelse(perm_nonrisk_total > 0, perm_nonrisk_share / perm_nonrisk_total * 100, 0)
    
    null_diffs[perm] <- perm_risk_pct - perm_nonrisk_pct
  }
  
  p_perm <- mean(null_diffs >= obs_diff)
  
  cat(sprintf("Observed difference: %.2f%% (Risk) - %.2f%% (Non-risk) = %.2f%%\n",
              obs_risk_pct, obs_nonrisk_pct, obs_diff))
  cat(sprintf("Permutation p-value: %.4f\n\n", p_perm))
  
  permutation_results[[i]] <- list(
    configuration = result$Configuration,
    database = result$Database,
    n_top = result$N_Top,
    obs_risk_pct = obs_risk_pct,
    obs_nonrisk_pct = obs_nonrisk_pct,
    obs_diff = obs_diff,
    null_diffs = null_diffs,
    p_fisher = result$P_Value,
    p_perm = p_perm,
    odds_ratio = result$Odds_Ratio
  )
}

# Create plots for each significant result 

for (i in 1:length(permutation_results)) {
  perm_res <- permutation_results[[i]]
  
  filename <- sprintf("./figures/pathway_sharing/permutation_%s_%s_top%d.pdf",
                      perm_res$configuration,
                      perm_res$database,
                      perm_res$n_top)
  
  pdf(filename, width = 10, height = 8)
  
  # create smooth density distribution from histogram
  h <- hist(perm_res$null_diffs, breaks = 50, plot = FALSE)
  
  config_name <- ifelse(perm_res$configuration == "1st_only",
                        "1st Degree Only",
                        "1st + 2nd Degree Combined")
  
  plot(h$mids, h$counts,
       type = "l",
       lwd = 2,
       col = "gray60",
       main = sprintf("%s: %s (TOP %d)\nRisk: %.2f%% | Non-Risk: %.2f%% | Odds Ratio: %.2f",
                      config_name,
                      perm_res$database,
                      perm_res$n_top,
                      perm_res$obs_risk_pct,
                      perm_res$obs_nonrisk_pct,
                      perm_res$odds_ratio),
       xlab = "Difference in Pathway Sharing % (Risk - Non-Risk)",
       ylab = "Frequency",
       xlim = range(c(perm_res$null_diffs, perm_res$obs_diff)),
       ylim = c(0, max(h$counts) * 1.2),
       cex.main = 0.95)
  
  # Certical line for observed
  abline(v = perm_res$obs_diff, col = "#d62728", lwd = 3)
  
  label_x <- perm_res$obs_diff
  label_y <- max(h$counts) * 1.1
  
  # text position
  text_pos <- ifelse(perm_res$obs_diff > mean(range(h$mids)), 2, 4)
  
  text(label_x,
       label_y,
       sprintf("Observed = %.2f%%\nFisher's p = %.4f",
               perm_res$obs_diff,
               perm_res$p_fisher),
       pos = text_pos,
       col = "#d62728",
       font = 2,
       cex = 1.2,
       adj = ifelse(text_pos == 2, 1, 0))
  
  dev.off()
  
  cat(sprintf("Saved: %s\n", basename(filename)))
}
# Read your network
network <- read_csv("./data/master_network_full_annotated.csv")

# Risk-specific genes (RED nodes)
risk_genes <- network %>%
  filter(specificity == "risk-specific") %>%
  pull(target_harmonized) %>%
  unique()

# Non-risk-specific genes (BLUE nodes)
nonrisk_genes <- network %>%
  filter(specificity == "nonrisk-specific") %>%
  pull(target_harmonized) %>%
  unique()

# Save to file
write.table(risk_genes, "risk_gene_list.txt", 
            quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(nonrisk_genes, "nonrisk_gene_list.txt", 
            quote = FALSE, row.names = FALSE, col.names = FALSE)

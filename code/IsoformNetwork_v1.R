# CAD Isoform Network Visualization v1

library(tidyverse)
library(igraph)
library(ggraph)

# Read data
network <- read_csv("./data/master_network_full_annotated.csv")

# ============================================================================
# 1. SUMMARY TABLE (2 minutes)
# ============================================================================

summary_table <- network %>%
  group_by(source_harmonized) %>%
  summarize(
    Total = n(),
    Risk = sum(specificity == "risk-specific"),
    NonRisk = sum(specificity == "nonrisk-specific"),
    Shared = sum(specificity == "shared"),
    Jaccard = round(Shared / (Risk + NonRisk + Shared), 2)
  ) %>%
  arrange(desc(Total))

write_csv(summary_table, "QUICK_summary_table.csv")
print(summary_table)

# ============================================================================
# 2. SINGLE GENE NETWORK - ITGA1 (3 minutes)
# ============================================================================

plot_simple_network <- function(data, gene) {
  
  # Filter
  df <- data %>% filter(source_harmonized == gene)
  
  # Create graph
  g <- graph_from_data_frame(
    df %>% select(source_harmonized, target_harmonized),
    directed = FALSE
  )
  
  # Colors
  node_colors <- c(gene = "#9370DB")  # Bait = purple
  for(i in 1:nrow(df)) {
    target <- df$target_harmonized[i]
    spec <- df$specificity[i]
    node_colors[target] <- case_when(
      spec == "risk-specific" ~ "#E74C3C",      # Red
      spec == "nonrisk-specific" ~ "#3498DB",   # Blue
      spec == "shared" ~ "#95A5A6"              # Gray
    )
  }
  
  # Layout
  set.seed(123)
  l <- layout_with_fr(g)
  
  # Node sizes
  node_sizes <- ifelse(V(g)$name == gene, 15, 8)
  
  # Plot
  pdf(paste0("QUICK_", gene, "_network.pdf"), width = 10, height = 8)
  plot(g, 
       vertex.color = node_colors[V(g)$name],
       vertex.size = node_sizes,
       vertex.label.cex = 0.7,
       vertex.label.color = "black",
       vertex.frame.color = "white",
       edge.color = "gray70",
       edge.width = 1,
       layout = l,
       main = paste(gene, "Interaction Network"))
  legend("topright", 
         legend = c("Bait", "Risk-specific", "Non-risk-specific", "Shared"),
         col = c("#9370DB", "#E74C3C", "#3498DB", "#95A5A6"),
         pch = 19, pt.cex = 2, cex = 0.8)
  dev.off()
  
  cat(paste("✓", gene, "network saved\n"))
}

# Make networks for key genes
plot_simple_network(network, "ITGA1")
plot_simple_network(network, "SRP54")
plot_simple_network(network, "ADAMTS7")

# ============================================================================
# 3. QUICK BAR CHART
# ============================================================================

pdf("QUICK_interaction_barplot.pdf", width = 10, height = 6)
barplot(
  t(as.matrix(summary_table[, c("Risk", "NonRisk", "Shared")])),
  names.arg = summary_table$source_harmonized,
  col = c("#E74C3C", "#3498DB", "#95A5A6"),
  legend = c("Risk-specific", "Non-risk-specific", "Shared"),
  main = "Interaction Profile by Gene",
  xlab = "Gene",
  ylab = "Number of Interactions",
  las = 2,
  cex.names = 0.8
)
dev.off()

cat("\n✓ All quick visualizations complete!\n")
cat("Files created:\n")
cat("  - QUICK_summary_table.csv\n")
cat("  - QUICK_ITGA1_network.pdf\n")
cat("  - QUICK_SRP54_network.pdf\n")
cat("  - QUICK_ADAMTS7_network.pdf\n")
cat("  - QUICK_interaction_barplot.pdf\n")
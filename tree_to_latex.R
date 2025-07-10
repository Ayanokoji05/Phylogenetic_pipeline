#!/usr/bin/env Rscript

# Load required libraries
if (!requireNamespace("ggtree", quietly = TRUE)) {
    install.packages("BiocManager")
    BiocManager::install("ggtree")
}
if (!requireNamespace("treeio", quietly = TRUE)) {
    BiocManager::install("treeio")
}
if (!requireNamespace("ggplot2", quietly = TRUE)) {
    install.packages("ggplot2")
}

library(ggtree)
library(treeio)
library(ggplot2)

# === CONFIG ===
newick_file <- "annotated_tree.nwk"
output_pdf <- "tree_plot.pdf"
dpi <- 300

# === READ TREE ===
tree <- read.tree(newick_file)

# === PARSE TIP LABELS ===
split_labels <- strsplit(tree$tip.label, "\\|")
label_df <- data.frame(
  original = tree$tip.label,
  transcript_id = sapply(split_labels, `[`, 1),
  gene_symbol = sapply(split_labels, `[`, 2),
  species = sapply(split_labels, `[`, 3),
  stringsAsFactors = FALSE
)

# Build final label
label_df$new_label <- paste(label_df$transcript_id,
                            label_df$gene_symbol,
                            label_df$species,
                            sep = " | ")

# === DETERMINE LABEL WIDTH IN cm ===
# Estimate PDF width in cm (not tree scale!)
max_label_chars <- max(nchar(label_df$new_label))
label_width_cm <- max_label_chars * 0.18  # ~0.18cm per character at 10-12pt font
tree_panel_width_cm <- 6  # space for tree structure

plot_width_cm <- tree_panel_width_cm + label_width_cm  # total width
n_tips <- length(tree$tip.label)
line_height_cm <- 0.25
plot_height_cm <- max(5, n_tips * line_height_cm)

# === DRAW TREE WITH ALIGNED LABELS ===
p <- ggtree(tree) %<+% label_df +
  geom_tiplab(aes(label = new_label),
              align = TRUE,
              linetype = "dotted",
              linesize = 0.3,
              size = 3,
              offset = 0,      # required with align=TRUE
              clip = "off") +
  theme(plot.margin = margin(2, 2, 2, 2, unit = "cm"))

# === EXPORT ===
ggsave(output_pdf,
       plot = p,
       width = plot_width_cm,
       height = plot_height_cm,
       units = "cm",
       dpi = dpi,
       device = cairo_pdf,
       limitsize = FALSE)

cat("✅ Tree with labels saved to:", output_pdf, "\n")
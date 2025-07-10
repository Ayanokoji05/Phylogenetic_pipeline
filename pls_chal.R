#!/usr/bin/env Rscript

#' Create a PDF with a large phylogenetic tree that may span multiple pages
#'
#' @param tree_file Path to input Newick tree file (e.g., "annotated_tree.nwk")
#' @param output_file Path to output PDF file (default: "large_tree.pdf")
#' @param caption Figure caption (default: "Phylogenetic tree")
#' @param label Figure label for LaTeX (default: "fig:tree")
#' @param width Page width in inches (default: 8.5)
#' @param height Page height in inches (default: 11)
#' @param dpi Resolution for the output (default: 300)
#' 
#' @return Generates a PDF file and prints LaTeX code to include it
create_large_tree_pdf <- function(tree_file = "annotated_tree.nwk",
                                  output_file = "large_tree.pdf",
                                  caption = "Phylogenetic tree",
                                  label = "fig:tree",
                                  width = 8.5,
                                  height = 11,
                                  dpi = 300) {
  
  # Load required packages
  if (!requireNamespace("ape", quietly = TRUE)) install.packages("ape")
  if (!requireNamespace("ggtree", quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
    BiocManager::install("ggtree")
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
  if (!requireNamespace("gridExtra", quietly = TRUE)) install.packages("gridExtra")
  
  library(ape)
  library(ggtree)
  library(ggplot2)
  library(gridExtra)
  
  # Read the tree file
  tree <- read.tree(tree_file)
  
  # Create the tree plot
  p <- ggtree(tree) + 
    geom_tiplab(size = 2) +  # Adjust size as needed
    theme_tree2() +
    ggplot2::labs(caption = caption)
  
  # Calculate the aspect ratio of the tree
  tree_ratio <- length(tree$tip.label) / max(node.depth(tree))
  
  # Standard paper ratio
  paper_ratio <- width / height
  
  # Determine if we need to adjust the plot dimensions
  if (tree_ratio > paper_ratio) {
    # Tree is wider than tall relative to paper - adjust width
    plot_width <- width
    plot_height <- width / tree_ratio
  } else {
    # Tree is taller than wide relative to paper - adjust height
    plot_height <- height
    plot_width <- height * tree_ratio
  }
  
  # Save to PDF
  ggsave(output_file, p, 
         width = plot_width, 
         height = plot_height, 
         units = "in",
         dpi = dpi,
         limitsize = FALSE)
  
  # Generate LaTeX code for inclusion
  latex_code <- sprintf('
\\begin{figure}[p]
    \\centering
    \\includegraphics[width=\\textwidth,height=\\textheight,keepaspectratio]{%s}
    \\caption{%s}
    \\label{%s}
\\end{figure}
', output_file, caption, label)
  
  cat("\nLaTeX code to include in your document:\n")
  cat(latex_code)
  
  # Return the plot object in case further modification is needed
  invisible(p)
}

# Example usage:
# create_large_tree_pdf(tree_file = "annotated_tree.nwk")
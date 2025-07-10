#!/usr/bin/env Rscript

# Load libraries with error handling
suppressPackageStartupMessages({
  if (!require(ape, quietly = TRUE)) {
    stop("Package 'ape' is required but not installed. Please install it with: install.packages('ape')")
  }
  if (!require(ggtree, quietly = TRUE)) {
    stop("Package 'ggtree' is required but not installed. Please install it with: BiocManager::install('ggtree')")
  }
  if (!require(ggplot2, quietly = TRUE)) {
    stop("Package 'ggplot2' is required but not installed. Please install it with: install.packages('ggplot2')")
  }
})

# Check command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  cat("Usage: Rscript script.R <input_tree_file> <output_file>\n")
  cat("Example: Rscript script.R final_tree.nwk output.pdf\n")
  cat("Note: Recommended output formats are PDF or SVG for large trees\n")
  quit(status = 1)
}

infile  <- args[1]  # e.g., "final_tree.nwk"
outfile <- args[2]  # e.g., "output.pdf"

# Validate input file exists
if (!file.exists(infile)) {
  stop(paste("Input file does not exist:", infile))
}

# Load tree with error handling
tryCatch({
  tree <- read.tree(infile)
}, error = function(e) {
  stop(paste("Error reading tree file:", e$message))
})

# Validate tree object
if (is.null(tree) || length(tree$tip.label) == 0) {
  stop("Invalid tree or empty tree loaded")
}

# Validate and root the tree if needed
if (!is.rooted(tree)) {
  cat("Tree is not rooted. Rooting at midpoint...\n")
  tree <- midpoint.root(tree)
}

# Remove any problematic characters from tip labels
tree$tip.label <- gsub("'", "", tree$tip.label)  # Remove single quotes
tree$tip.label <- gsub("^\\s+|\\s+$", "", tree$tip.label)  # Remove leading/trailing spaces

# Calculate dimensions
max_label_length <- max(nchar(tree$tip.label))
tip_count <- length(tree$tip.label)
tree_depth <- max(node.depth.edgelength(tree))

# Create the tree plot
p <- ggtree(tree, layout = "rectangular", branch.length = "branch.length") +
  geom_tiplab(align = TRUE, 
              linetype = "dotted",
              size = 3, 
              offset = tree_depth * 0.05,
              family = "mono",
              hjust = 0) +
  theme_tree2() +
  theme(
    plot.margin = margin(30, 50, 30, 30),
    axis.text.x = element_text(size = 10),
    axis.title.x = element_text(size = 12)
  ) +
  xlim(-tree_depth * 0.1, tree_depth * 1.8) +
  xlab("Branch Length")

# Smart dimension calculation for large trees
calculate_dimensions <- function(tip_count, max_label_length, tree_depth) {
  # Base dimensions
  base_height <- 8
  base_width <- 15
  
  # Dynamic scaling factors (adjusted for large trees)
  height_scale <- ifelse(tip_count > 500, 0.12, 0.2)
  width_scale <- ifelse(tip_count > 500, 0.06, 0.12)
  
  # Calculate dimensions with caps
  img_height <- min(120, base_height + (tip_count * height_scale))
  img_width <- min(48, base_width + (tree_depth * 0.8) + (max_label_length * width_scale))
  
  return(list(width = img_width, height = img_height))
}

dims <- calculate_dimensions(tip_count, max_label_length, tree_depth)
img_width <- dims$width
img_height <- dims$height

cat(paste("Tree has", tip_count, "tips\n"))
cat(paste("Maximum label length:", max_label_length, "characters\n"))
cat(paste("Tree depth:", round(tree_depth, 4), "\n"))
cat(paste("Output image dimensions:", round(img_width, 2), "x", round(img_height, 2), "inches\n"))

# Clean up any existing graphics devices
while (length(dev.list()) > 0) {
  dev.off()
}

# Determine output format
output_format <- tools::file_ext(outfile)

# Save the plot with format-specific handling
tryCatch({
  if (tolower(output_format) %in% c("pdf", "svg")) {
    # Vector formats (recommended for large trees)
    ggsave(outfile, plot = p,
           width = img_width, height = img_height,
           units = "in", limitsize = FALSE)
    cat(paste("Successfully saved vector graphic:", outfile, "\n"))
    
  } else if (tolower(output_format) == "png") {
    # Raster format with safety checks
    max_pixels <- 32768
    required_pixels <- max(img_width * 300, img_height * 300)
    
    if (required_pixels > max_pixels) {
      adjusted_dpi <- floor(max_pixels / max(img_width, img_height))
      warning(paste("Requested PNG dimensions too large. Reducing DPI to", adjusted_dpi))
      ggsave(outfile, plot = p,
             width = img_width, height = img_height,
             units = "in", dpi = adjusted_dpi, limitsize = FALSE)
    } else {
      ggsave(outfile, plot = p,
             width = img_width, height = img_height,
             units = "in", dpi = 300, limitsize = FALSE)
    }
    cat(paste("PNG saved successfully:", outfile, "\n"))
    
  } else {
    # Default to PDF if format not recognized
    pdf_file <- sub(paste0("\\.", output_format, "$"), ".pdf", outfile)
    warning(paste("Unsupported format", output_format, "- saving as PDF instead"))
    ggsave(pdf_file, plot = p,
           width = img_width, height = img_height,
           units = "in", limitsize = FALSE)
    cat(paste("Saved as PDF:", pdf_file, "\n"))
  }
}, error = function(e) {
  stop(paste("Failed to save output:", e$message))
})

cat("Script completed successfully!\n")
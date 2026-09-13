#!/usr/bin/env Rscript

# Plot structural-variant length distributions from two VCF files.
# Usage: plot_distributions.R <file1.vcf> <file2.vcf> [-o <output-directory>]

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)
usage <- "Usage: plot_distributions.R <file1.vcf> <file2.vcf> [-o <output-directory>]"

if (length(args) == 1 && args[[1]] %in% c("-h", "--help")) {
  cat(usage, "\n")
  quit(status = 0)
}

output_dir <- "."
if (length(args) == 4 && args[[3]] %in% c("-o", "--output")) {
  output_dir <- args[[4]]
  args <- args[1:2]
}
if (length(args) != 2) stop(usage, call. = FALSE)

vcf_paths <- args
missing_files <- vcf_paths[!file.exists(vcf_paths)]
if (length(missing_files) > 0) {
  stop("VCF file not found: ", paste(missing_files, collapse = ", "), call. = FALSE)
}

read_sv_lengths <- function(path) {
  vcf <- tryCatch(
    read.delim(
      path, comment.char = "#", header = FALSE, sep = "\t", quote = "",
      stringsAsFactors = FALSE, fill = TRUE
    ),
    error = function(e) stop("Could not read ", path, ": ", conditionMessage(e), call. = FALSE)
  )

  if (nrow(vcf) == 0 || ncol(vcf) < 8) {
    stop("VCF has no records with the eight required columns: ", path, call. = FALSE)
  }

  info <- vcf[[8]]
  matched <- regexec("(?:^|;)SVLEN=([^;]+)", info, perl = TRUE)
  raw_values <- regmatches(info, matched)
  raw_values <- unlist(lapply(raw_values, function(value) {
    if (length(value) < 2) return(character())
    strsplit(value[[2]], ",", fixed = TRUE)[[1]]
  }), use.names = FALSE)
  lengths <- suppressWarnings(abs(as.numeric(raw_values)))
  lengths <- lengths[is.finite(lengths)]

  if (length(lengths) == 0) {
    stop("VCF contains no numeric SVLEN values: ", path, call. = FALSE)
  }

  lengths
}

make_distribution_plot <- function(lengths, color, show_y_label, x_limits = NULL) {
  data <- data.frame(length = lengths)

  ggplot(data, aes(x = length)) +
    geom_histogram(
      aes(y = after_stat(density)), bins = 60,
      fill = color, color = "white", alpha = 0.45, linewidth = 0.2
    ) +
    geom_density(color = color, linewidth = 1, adjust = 0.8) +
    labs(
      x = "Absolute SV length (bp)",
      y = if (show_y_label) "Density" else NULL
    ) +
    scale_x_continuous(
      breaks = scales::breaks_pretty(n = 5),
      labels = scales::label_comma()
    ) +
    scale_y_continuous(labels = NULL) +
    coord_cartesian(xlim = x_limits) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      axis.title.x = element_text(margin = margin(t = 10))
    )
}

length_sets <- lapply(vcf_paths, read_sv_lengths)
# Okabe-Ito colorblind-safe palette: blue and vermilion.
plot_colors <- c("#4B0092", "#16AD16")
plots <- Map(
  make_distribution_plot,
  length_sets,
  plot_colors,
  c(TRUE, FALSE),
  list(c(0, 50000), NULL)
)
combined_plot <- plots[[1]] + plots[[2]] + plot_layout(ncol = 2)

if (!dir.exists(output_dir) && !dir.create(output_dir, recursive = TRUE)) {
  stop("Could not create output directory: ", output_dir, call. = FALSE)
}
output_path <- file.path(output_dir, "distributions.png")
ggsave(output_path, combined_plot, width = 12, height = 5.5, dpi = 300)

message("Saved distribution plot to ", normalizePath(output_path))

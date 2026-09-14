#!/usr/bin/env Rscript
# Run from the main project folder:
# pixi run Rscript fig3-dotplot/compare_dotplots.R fig3-dotplot/reference-vs-reference.paf fig3-dotplot/simulated-vs-reference.paf -o fig3-dotplot/control-vs-simulated-dotplot.png

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)
usage <- paste(
  "Usage: compare_dotplots.R <control.paf> <simulated.paf>",
  "[-o <output.png>] [--width <inches>] [--height <inches>] [--dpi <number>]"
)

if (length(args) == 1L && args[[1]] %in% c("-h", "--help")) {
  cat(usage, "\n")
  quit(status = 0L)
}

options <- list(
  output = "fig3-dotplot/control-vs-simulated-dotplot.png",
  width = 18,
  height = 9,
  dpi = 200
)
inputs <- character()
i <- 1L
while (i <= length(args)) {
  option_names <- c("-o" = "output", "--output" = "output",
                    "--width" = "width", "--height" = "height", "--dpi" = "dpi")
  if (args[[i]] %in% names(option_names)) {
    if (i == length(args)) stop("Missing value for ", args[[i]], call. = FALSE)
    name <- unname(option_names[args[[i]]])
    value <- args[[i + 1L]]
    if (name != "output") {
      value <- suppressWarnings(as.numeric(value))
      if (!is.finite(value) || value <= 0) stop("Invalid value for ", args[[i]], call. = FALSE)
    }
    options[[name]] <- value
    i <- i + 2L
  } else {
    if (startsWith(args[[i]], "-")) stop("Unknown option: ", args[[i]], call. = FALSE)
    inputs <- c(inputs, args[[i]])
    i <- i + 1L
  }
}

if (length(inputs) != 2L) stop(usage, call. = FALSE)
if (any(!file.exists(inputs))) {
  stop("PAF file not found: ", paste(inputs[!file.exists(inputs)], collapse = ", "), call. = FALSE)
}
if (!grepl("\\.png$", options$output, ignore.case = TRUE)) options$output <- paste0(options$output, ".png")

# Locate the shared plotting functions whether this script is called from the
# project folder or through an absolute path.
script_arg <- commandArgs()[startsWith(commandArgs(), "--file=")][1]
script_path <- normalizePath(sub("^--file=", "", script_arg))
source(file.path(dirname(dirname(script_path)), "scripts", "dotplot.R"))

prepare_alignment <- function(path) {
  all_alignments <- read_paf(path)
  alignments <- all_alignments[all_alignments$alen >= 1000, ]
  if (!nrow(alignments)) stop("No alignments remain in ", path, call. = FALSE)
  alignments <- classify(alignments, min_sv = 50)
  list(
    all = all_alignments,
    alignments = alignments,
    segments = make_segments(alignments, min_sv = 50)
  )
}

control <- prepare_alignment(inputs[[1]])
simulated <- prepare_alignment(inputs[[2]])

# Use one range for both panels. This makes changes in position and genome size
# directly comparable instead of allowing each panel to stretch independently.
axis_total <- function(paf, name, length_name) {
  sum(paf[[length_name]][!duplicated(paf[[name]])])
}
x_limit <- max(axis_total(control$all, "t", "tlen"),
               axis_total(simulated$all, "t", "tlen"))
y_limit <- max(axis_total(control$all, "q", "qlen"),
               axis_total(simulated$all, "q", "qlen"))

plot_options <- list(
  display_contigs = TRUE,
  title = NULL,
  format = character(),
  output = "",
  width = options$width,
  dpi = options$dpi
)

control_plot <- plot_dotplot(
  control$alignments, control$segments, plot_options, control$all,
  x_limit = x_limit, y_limit = y_limit
) + labs(title = "Control") + guides(colour = "none")

simulated_plot <- plot_dotplot(
  simulated$alignments, simulated$segments, plot_options, simulated$all,
  x_limit = x_limit, y_limit = y_limit
) + labs(title = "Simulated")

comparison <- (control_plot + simulated_plot) +
  plot_layout(ncol = 2, guides = "collect") &
  theme(legend.position = "bottom")

dir.create(dirname(options$output), recursive = TRUE, showWarnings = FALSE)
ggsave(
  options$output,
  comparison,
  width = options$width,
  height = options$height,
  units = "in",
  dpi = options$dpi,
  bg = "white"
)
message("Wrote ", options$output)

#!/usr/bin/env Rscript

# Plot an inSVert simulated VCF as a Circos-style structural-variant overview.
# Usage:
#   R_LIBS_USER=viz/Rlib Rscript viz/plot_circlize.R \
#     viz/simulated.vcf data/human/Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai \
#     viz/insvert_sv_circos [--title "Plot title"] [--format pdf png svg]

suppressPackageStartupMessages(library(circlize))

args <- commandArgs(trailingOnly = TRUE)
usage <- "Usage: plot_circlize.R <simulated.vcf> <reference.fa.fai> <output-prefix> [--title <title>] [--format <pdf|png|svg> ...]"
if (length(args) < 3) stop(usage)
vcf_path <- args[[1]]
fai_path <- args[[2]]
output_prefix <- args[[3]]
plot_title <- ""
formats <- "png"
i <- 4
while (i <= length(args)) {
  if (args[[i]] == "--title" && i < length(args) && !startsWith(args[[i + 1]], "--")) {
    plot_title <- args[[i + 1]]
    i <- i + 2
  } else if (args[[i]] == "--format") {
    i <- i + 1
    start <- i
    while (i <= length(args) && !startsWith(args[[i]], "--")) i <- i + 1
    if (i == start) stop(usage)
    formats <- args[start:(i - 1)]
    if (any(!formats %in% c("pdf", "png", "svg"))) stop(usage)
  } else {
    stop(usage)
  }
}

read_info <- function(info, key) {
  pattern <- paste0("(?:^|;)", key, "=([^;]+)")
  match <- regexec(pattern, info, perl = TRUE)
  values <- regmatches(info, match)
  vapply(values, function(x) if (length(x) >= 2) x[[2]] else NA_character_, character(1))
}

vcf <- read.delim(
  vcf_path, comment.char = "#", header = FALSE, sep = "\t",
  quote = "", stringsAsFactors = FALSE
)
if (ncol(vcf) < 8) stop("VCF has fewer than eight required columns")
names(vcf)[1:8] <- c("chrom", "pos", "id", "ref", "alt", "qual", "filter", "info")
vcf$pos <- as.numeric(vcf$pos)
vcf$svtype <- read_info(vcf$info, "SVTYPE")
vcf$end <- suppressWarnings(as.numeric(read_info(vcf$info, "END")))
vcf$svlen <- suppressWarnings(abs(as.numeric(read_info(vcf$info, "SVLEN"))))
vcf$event <- read_info(vcf$info, "EVENT")

fai <- read.delim(fai_path, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
names(fai)[1:2] <- c("chrom", "length")
canonical <- c(as.character(1:22), "X", "Y")
fai <- fai[match(canonical, fai$chrom, nomatch = 0), c("chrom", "length")]
if (nrow(fai) != length(canonical)) stop("Reference FAI does not contain all chromosomes 1-22, X and Y")

sv_colors <- c(
  DEL = "#9e2a2b", INS = "#e09f3e", DUP = "#009E73", INV = "#7A5AA6",
  TRA_COPY = "#0072B2", TRA_CUT = "#C44E52"
)

# Convert ordinary symbolic SV records into genomic intervals.
ordinary <- vcf[vcf$svtype %in% c("DEL", "INS", "DUP", "INV") & vcf$chrom %in% fai$chrom, ]
ordinary$plot_end <- ordinary$end
ordinary$plot_end[is.na(ordinary$plot_end)] <- ordinary$pos[is.na(ordinary$plot_end)] + ordinary$svlen[is.na(ordinary$plot_end)]
ordinary$plot_end[is.na(ordinary$plot_end) | ordinary$plot_end <= ordinary$pos] <-
  ordinary$pos[is.na(ordinary$plot_end) | ordinary$plot_end <= ordinary$pos] + 1

# inSVert translocations are BND groups. P1 is always at the destination;
# the other P records identify the source interval. Reverse events swap P3/P4.
bnd <- vcf[vcf$svtype == "BND" & !is.na(vcf$event), ]
event_ids <- unique(bnd$event)
links <- data.frame()
for (event_id in event_ids) {
  rows <- bnd[bnd$event == event_id, ]
  p1 <- rows[grepl("\\.P1$", rows$id), ]
  if (nrow(p1) == 0) next
  p1 <- p1[1, ]
  is_reverse <- grepl("]", p1$alt, fixed = TRUE)
  source_suffixes <- if (is_reverse) c(".P2", ".P4") else c(".P2", ".P3")
  source_rows <- rows[vapply(rows$id, function(x) any(endsWith(x, source_suffixes)), logical(1)), ]
  if (nrow(source_rows) < 2) next
  source_rows <- source_rows[source_rows$chrom == source_rows$chrom[[1]], ]
  if (nrow(source_rows) < 2) next
  source_start <- min(source_rows$pos)
  source_end <- max(source_rows$pos)
  destination_width <- max(250000, min(source_end - source_start, 1000000))
  event_type <- if (grepl("TRA_CUT", event_id, fixed = TRUE)) "TRA_CUT" else "TRA_COPY"
  links <- rbind(links, data.frame(
    event = event_id, type = event_type,
    source_chrom = source_rows$chrom[[1]], source_start = source_start, source_end = source_end,
    destination_chrom = p1$chrom, destination_start = p1$pos,
    destination_end = p1$pos + destination_width,
    stringsAsFactors = FALSE
  ))
}
links <- links[links$source_chrom %in% fai$chrom & links$destination_chrom %in% fai$chrom, ]
links$source_end <- pmin(links$source_end, fai$length[match(links$source_chrom, fai$chrom)])
links$destination_end <- pmin(
  links$destination_end,
  fai$length[match(links$destination_chrom, fai$chrom)]
)

draw_plot <- function() {
  # Keep the circular plot square and reserve the extra device height for the legend.
  layout(matrix(c(1, 2), ncol = 1), heights = c(10, 1))
  par(mar = c(0.5, 0.5, 2, 0.5))
  circos.clear()
  circos.par(
    start.degree = 90, gap.after = c(rep(1.2, 21), 3, 1.2, 4),
    track.margin = c(0.003, 0.003), cell.padding = c(0, 0, 0, 0),
    points.overflow.warning = FALSE
  )
  circos.initialize(
    factors = factor(fai$chrom, levels = fai$chrom),
    xlim = cbind(rep(0, nrow(fai)), fai$length)
  )

  # Outer chromosome ideogram and labels.
  circos.trackPlotRegion(ylim = c(0, 1), track.height = 0.065, bg.border = "#555555", panel.fun = function(x, y) {
    chr <- CELL_META$sector.index
    circos.text(CELL_META$xcenter, 1.7, paste0("chr", chr), cex = 0.48,
                facing = "bending.outside", niceFacing = TRUE)
  })

  # Four inner rings: one for each non-translocation SV type.
  ring_types <- c("INV", "DUP", "DEL", "INS")
  for (sv in ring_types) {
    circos.trackPlotRegion(ylim = c(0, 1), track.height = 0.045, bg.col = "#F5F5F5", bg.border = "#DDDDDD",
      panel.fun = function(x, y) {
        chr <- CELL_META$sector.index
        rows <- ordinary[ordinary$svtype == sv & ordinary$chrom == chr, ]
        if (nrow(rows) == 0) return()
        for (i in seq_len(nrow(rows))) {
          start <- max(rows$pos[[i]], CELL_META$xlim[[1]])
          end <- min(rows$plot_end[[i]], CELL_META$xlim[[2]])
          if (sv == "INS") {
            circos.segments(start, 0.12, start, 0.88, col = sv_colors[[sv]], lwd = 0.7)
          } else if (end > start) {
            circos.rect(start, 0.12, end, 0.88, col = sv_colors[[sv]], border = NA)
          }
        }
      })
  }

  # Draw translocations last so their bands remain visible inside the rings.
  if (nrow(links) > 0) {
    for (i in seq_len(nrow(links))) {
      link <- links[i, ]
      circos.link(
        link$source_chrom, c(link$source_start, link$source_end),
        link$destination_chrom, c(link$destination_start, link$destination_end),
        col = adjustcolor(sv_colors[[link$type]], alpha.f = 0.34),
        border = adjustcolor(sv_colors[[link$type]], alpha.f = 0.65),
        lwd = 0.35
      )
    }
  }

  if (nzchar(plot_title)) title(plot_title, line = -1.2, cex.main = 1.05)
  circos.clear()

  par(mar = c(0, 0, 0, 0))
  plot.new()
  legend(
    "center", legend = c("Cut translocation", "Copy translocation", "Inversion", "Duplication", "Deletion", "Insertion"),
    fill = sv_colors[c("TRA_CUT", "TRA_COPY", "INV", "DUP", "DEL", "INS")],
    ncol = 3, bty = "n", cex = 0.75
  )
}

for (format in unique(formats)) {
  if (format == "svg") svg(paste0(output_prefix, ".svg"), width = 10, height = 11)
  if (format == "pdf") pdf(paste0(output_prefix, ".pdf"), width = 10, height = 11, useDingbats = FALSE)
  if (format == "png") png(paste0(output_prefix, ".png"), width = 3000, height = 3300, res = 300)
  draw_plot()
  dev.off()
}

message(sprintf(
  "Plotted %d translocations and %d other SVs (%s)",
  nrow(links), nrow(ordinary), paste(names(table(ordinary$svtype)), table(ordinary$svtype), collapse = ", ")
))

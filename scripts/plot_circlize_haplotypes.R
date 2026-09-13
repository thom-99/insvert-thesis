#!/usr/bin/env Rscript

# Plot each haplotype from an inSVert simulated VCF as horizontally paired
# Circos-style structural-variant overviews.
# Usage:
#   Rscript scripts/plot_circlize_haplotypes.R \
#     simulated.vcf reference.fa.fai output-prefix \
#     [--ploidy 2] [--format pdf png svg]

suppressPackageStartupMessages({
  library(circlize)
  library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)
usage <- "Usage: plot_circlize_haplotypes.R <simulated.vcf> <reference.fa.fai> <output-prefix> [--ploidy <N>] [--format <pdf|png|svg> ...]"
if (length(args) < 3) stop(usage)
vcf_path <- args[[1]]
fai_path <- args[[2]]
output_prefix <- args[[3]]
ploidy <- NULL
formats <- "png"
i <- 4
while (i <= length(args)) {
  if (args[[i]] == "--ploidy") {
    if (i == length(args) || !grepl("^[1-9][0-9]*$", args[[i + 1]])) {
      stop("--ploidy must be a positive integer")
    }
    ploidy <- as.integer(args[[i + 1]])
    i <- i + 2
  } else if (args[[i]] == "--format") {
    i <- i + 1
    start <- i
    while (i <= length(args) && !startsWith(args[[i]], "--")) i <- i + 1
    if (i == start) stop("--format requires at least one of: pdf, png, svg")
    formats <- args[start:(i - 1)]
    if (any(!formats %in% c("pdf", "png", "svg"))) {
      stop("Unsupported format; choose one or more of: pdf, png, svg")
    }
  } else {
    stop(paste("Unknown argument:", args[[i]]))
  }
}

if (!file.exists(vcf_path)) stop(paste("VCF file not found:", vcf_path))
if (!file.exists(fai_path)) stop(paste("FAI file not found:", fai_path))
if (!dir.exists(dirname(output_prefix))) {
  stop(paste("Output directory does not exist:", dirname(output_prefix)))
}

vcf_first <- readLines(vcf_path, n = 1, warn = FALSE)
if (length(vcf_first) == 0 || !startsWith(vcf_first, "##fileformat=VCF")) {
  if (length(vcf_first) > 0 && startsWith(vcf_first, ">")) {
    stop("The VCF argument appears to be a FASTA file")
  }
  stop("The VCF argument does not appear to be a VCF file")
}

fai_first <- readLines(fai_path, n = 1, warn = FALSE)
if (length(fai_first) == 0) stop("The FAI file is empty")
if (startsWith(fai_first, ">")) {
  stop("The FAI argument appears to be a FASTA file; provide its .fai index")
}
fai_fields <- strsplit(fai_first, "\t", fixed = TRUE)[[1]]
if (length(fai_fields) < 2 || is.na(suppressWarnings(as.numeric(fai_fields[[2]])))) {
  stop("The FAI argument does not appear to be a valid FASTA index")
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
if (ncol(vcf) < 10) stop("VCF must contain FORMAT and one sample column")
names(vcf)[1:10] <- c(
  "chrom", "pos", "id", "ref", "alt", "qual", "filter", "info", "format", "sample"
)
vcf$pos <- as.numeric(vcf$pos)
vcf$svtype <- read_info(vcf$info, "SVTYPE")
vcf$end <- suppressWarnings(as.numeric(read_info(vcf$info, "END")))
vcf$svlen <- suppressWarnings(abs(as.numeric(read_info(vcf$info, "SVLEN"))))
vcf$event <- read_info(vcf$info, "EVENT")

read_gt <- function(format, sample) {
  keys <- strsplit(format, ":", fixed = TRUE)[[1]]
  values <- strsplit(sample, ":", fixed = TRUE)[[1]]
  gt_index <- match("GT", keys)
  if (is.na(gt_index) || gt_index > length(values)) return(character())
  strsplit(values[[gt_index]], "[|/]", perl = TRUE)[[1]]
}
genotypes <- Map(read_gt, vcf$format, vcf$sample)
if (is.null(ploidy)) {
  if (length(genotypes) == 0 || length(genotypes[[1]]) == 0) {
    stop("Cannot infer ploidy: the first VCF record has no GT")
  }
  ploidy <- length(genotypes[[1]])
}
if (any(lengths(genotypes) != ploidy)) {
  stop("Every VCF genotype must have the same allele count as --ploidy")
}
carrier <- t(vapply(genotypes, function(gt) gt == "1", logical(ploidy)))

fai <- read.delim(fai_path, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
names(fai)[1:2] <- c("chrom", "length")
canonical <- c(as.character(1:22), "X", "Y")
fai <- fai[match(canonical, fai$chrom, nomatch = 0), c("chrom", "length")]
if (nrow(fai) != length(canonical)) stop("Reference FAI does not contain all chromosomes 1-22, X and Y")

sv_colors <- c(
  DEL = "#9e2a2b", INS = "#e09f3e", DUP = "#009E73", INV = "#7A5AA6",
  TRA_COPY = "#0072B2", TRA_CUT = "#C44E52"
)

build_haplotype_data <- function(haplotype) {
  hap_vcf <- vcf[carrier[, haplotype], ]

  # Convert ordinary symbolic SV records into genomic intervals.
  ordinary <- hap_vcf[
    hap_vcf$svtype %in% c("DEL", "INS", "DUP", "INV") & hap_vcf$chrom %in% fai$chrom,
  ]
  ordinary$plot_end <- ordinary$end
  ordinary$plot_end[is.na(ordinary$plot_end)] <-
    ordinary$pos[is.na(ordinary$plot_end)] + ordinary$svlen[is.na(ordinary$plot_end)]
  ordinary$plot_end[is.na(ordinary$plot_end) | ordinary$plot_end <= ordinary$pos] <-
    ordinary$pos[is.na(ordinary$plot_end) | ordinary$plot_end <= ordinary$pos] + 1

  # inSVert translocations are BND groups. P1 is always at the destination;
  # the other P records identify the source interval. Reverse events swap P3/P4.
  bnd <- hap_vcf[hap_vcf$svtype == "BND" & !is.na(hap_vcf$event), ]
  event_ids <- unique(bnd$event)
  links <- data.frame(
    event = character(), type = character(),
    source_chrom = character(), source_start = numeric(), source_end = numeric(),
    destination_chrom = character(), destination_start = numeric(), destination_end = numeric(),
    stringsAsFactors = FALSE
  )
  for (event_id in event_ids) {
    rows <- bnd[bnd$event == event_id, ]
    p1 <- rows[grepl("\\.P1$", rows$id), ]
    if (nrow(p1) == 0) next
    p1 <- p1[1, ]
    is_reverse <- grepl("]", p1$alt, fixed = TRUE)
    source_suffixes <- if (is_reverse) c(".P2", ".P4") else c(".P2", ".P3")
    source_rows <- rows[
      vapply(rows$id, function(x) any(endsWith(x, source_suffixes)), logical(1)),
    ]
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
  links <- links[
    links$source_chrom %in% fai$chrom & links$destination_chrom %in% fai$chrom,
  ]
  if (nrow(links) > 0) {
    links$source_end <- pmin(
      links$source_end, fai$length[match(links$source_chrom, fai$chrom)]
    )
    links$destination_end <- pmin(
      links$destination_end,
      fai$length[match(links$destination_chrom, fai$chrom)]
    )
  }

  list(ordinary = ordinary, links = links)
}

haplotype_data <- lapply(seq_len(ploidy), build_haplotype_data)

draw_plot <- function(haplotype) {
  ordinary <- haplotype_data[[haplotype]]$ordinary
  links <- haplotype_data[[haplotype]]$links

  # Keep the circular plot square and reserve enough top margin for its title.
  par(mar = c(0.5, 0.5, 3, 0.5))
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
        lwd = 0.35, directional = 1,
        arr.type = "triangle", arr.length = 0.12, arr.width = 0.06
      )
    }
  }

  title(paste("haplotype", haplotype), line = 1, cex.main = 1.05)
  circos.clear()
}

panels <- lapply(seq_len(ploidy), function(haplotype) {
  local_haplotype <- haplotype
  wrap_elements(full = ~draw_plot(local_haplotype))
})
panel_row <- wrap_plots(panels, nrow = 1)
legend_labels <- c(
  "Cut translocation", "Copy translocation", "Inversion",
  "Duplication", "Deletion", "Insertion"
)
legend_colors <- unname(sv_colors[c("TRA_CUT", "TRA_COPY", "INV", "DUP", "DEL", "INS")])
legend_items <- lapply(seq_along(legend_labels), function(i) {
  grid::grobTree(
    grid::rectGrob(
      x = 0.38, width = 0.07, height = 0.22,
      gp = grid::gpar(fill = legend_colors[[i]], col = NA)
    ),
    grid::textGrob(
      legend_labels[[i]], x = 0.45, just = "left",
      gp = grid::gpar(fontsize = 9)
    ),
    vp = grid::viewport(layout.pos.col = i)
  )
})
legend_grob <- grid::gTree(
  children = do.call(grid::gList, legend_items),
  vp = grid::viewport(
    x = 0.5, width = 0.50,
    layout = grid::grid.layout(nrow = 1, ncol = length(legend_labels))
  )
)
legend_panel <- wrap_elements(full = legend_grob)
combined_plot <- wrap_plots(panel_row, legend_panel, ncol = 1, heights = c(10, 0.55))

for (format in unique(formats)) {
  if (format == "svg") svg(paste0(output_prefix, ".svg"), width = 10 * ploidy, height = 11)
  if (format == "pdf") pdf(paste0(output_prefix, ".pdf"), width = 10 * ploidy, height = 11, useDingbats = FALSE)
  if (format == "png") png(paste0(output_prefix, ".png"), width = 3000 * ploidy, height = 3300, res = 300)
  print(combined_plot)
  dev.off()
}

for (haplotype in seq_len(ploidy)) {
  ordinary <- haplotype_data[[haplotype]]$ordinary
  links <- haplotype_data[[haplotype]]$links
  counts <- table(ordinary$svtype)
  message(sprintf(
    "Haplotype %d: %d translocations and %d other SVs (%s)",
    haplotype, nrow(links), nrow(ordinary),
    paste(names(counts), counts, collapse = ", ")
  ))
}

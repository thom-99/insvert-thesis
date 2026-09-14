#!/usr/bin/env Rscript
# Run from the main project folder (uses ggplot2 from the pixi environment):
# pixi run Rscript fig3-dotplot/dotplot.R fig3-dotplot/simulated-vs-reference.paf -o fig3-dotplot/simulated-dotplot --format svg --display-contigs
# Add --title "Genome comparison" for an optional title.
# Positions and possible structural variants (SVs) are measured in the query genome
# relative to the reference genome. Colors suggest events; they do not prove them.
usage <- function() cat('Usage: Rscript dotplot.R alignment.paf[.gz] [options]
  -o, --output PREFIX    Output prefix (default: input path without .paf[.gz])
  --format png pdf svg   One or more formats (default: png)
  --display-contigs    Show contig labels and boundaries (default: off)
  --title TEXT         Optional plot title (default: none)
  --min-sv N            Minimum candidate event size in bp (default: 50)
  --min-aln N           Minimum alignment block length (default: 1000)
  --min-mapq N          Minimum MAPQ; 255 is unknown (default: 0)
  --width N            Figure width and height in inches (default: 10)
  --dpi N              PNG resolution (default: 200)
  -h, --help           Show this help
Use minimap2 -cx asm5 reference.fa query.fa > alignment.paf for CIGAR gaps.
Colors indicate candidate patterns, not validated SV calls. SAM/BAM unsupported.
')

# Read the input filename and command options, and check their values before plotting.
parse_args <- function(args) {
  opt <- list(input = NULL, output = NULL, format = 'png', min_sv = 50,
              min_aln = 1000, min_mapq = 0, width = 10, dpi = 200,
              display_contigs = FALSE, title = NULL)
  i <- 1L
  while (i <= length(args)) {
    a <- args[i]
    if (a %in% c('-h', '--help')) { usage(); return(NULL) }
    if (a == '--display-contigs') { opt$display_contigs <- TRUE; i <- i + 1L; next }
    # Collect all formats up to the next option, allowing --format svg pdf.
    if (a == '--format') {
      j <- i + 1L
      while (j <= length(args) && !startsWith(args[j], '-')) j <- j + 1L
      if (j == i + 1L) stop('--format needs png, pdf and/or svg.')
      opt$format <- unique(args[seq.int(i + 1L, j - 1L)])
      if (any(!opt$format %in% c('png', 'pdf', 'svg'))) stop('Unknown output format.')
      i <- j; next
    }
    keys <- c('-o' = 'output', '--output' = 'output', '--min-sv' = 'min_sv',
              '--min-aln' = 'min_aln', '--min-mapq' = 'min_mapq',
              '--width' = 'width', '--dpi' = 'dpi', '--title' = 'title')
    if (a %in% names(keys)) {
      if (i == length(args)) stop('Missing value for ', a)
      key <- unname(keys[a]); value <- args[i + 1L]
      if (!key %in% c('output', 'title')) {
        value <- suppressWarnings(as.numeric(value))
        if (!is.finite(value) || value < 0 ||
            (key %in% c('min_sv', 'width', 'dpi') && value == 0)) stop('Invalid value for ', a)
        if (key == 'min_mapq' && value > 254) stop('--min-mapq must be between 0 and 254.')
      }
      opt[[key]] <- value; i <- i + 2L; next
    }
    if (startsWith(a, '-') || !is.null(opt$input)) stop('Unexpected argument: ', a)
    opt$input <- a; i <- i + 1L
  }
  if (is.null(opt$input)) stop('An input PAF file is required. Use --help.')
  if (is.null(opt$output)) opt$output <- sub('\\.paf(\\.gz)?$', '', opt$input, ignore.case = TRUE)
  opt$output <- sub('\\.(png|pdf|svg)$', '', opt$output, ignore.case = TRUE)
  if (!nzchar(opt$output)) stop('Output prefix must not be empty.')
  opt
}

# PAF stores one aligned region per row. Read plain or compressed files identically.
read_paf <- function(path) {
  con <- if (grepl('\\.gz$', path)) gzfile(path, 'rt') else file(path, 'rt')
  on.exit(close(con))
  lines <- readLines(con, warn = FALSE)
  lines <- lines[nzchar(lines) & !startsWith(lines, '#')]
  if (!length(lines)) stop('The PAF file is empty.')
  fields <- strsplit(lines, '\t', fixed = TRUE)
  if (any(lengths(fields) < 12L)) stop('PAF requires at least 12 tab-separated columns per row.')
  p <- as.data.frame(do.call(rbind, lapply(fields, head, 12L)), stringsAsFactors = FALSE)
  # q = query, t = reference; s/e = start/end; len = full sequence length.
  # PAF positions start at zero, and the end position is excluded from the region.
  names(p) <- c('q', 'qlen', 'qs', 'qe', 'strand', 't', 'tlen', 'ts', 'te', 'matches', 'alen', 'mapq')
  for (n in c('qlen', 'qs', 'qe', 'tlen', 'ts', 'te', 'matches', 'alen', 'mapq')) {
    p[[n]] <- suppressWarnings(as.numeric(p[[n]]))
    if (any(!is.finite(p[[n]]) | p[[n]] < 0 | p[[n]] != floor(p[[n]]))) stop('Invalid PAF column: ', n)
  }
  if (any(!p$strand %in% c('+', '-') | p$qs >= p$qe | p$qe > p$qlen |
          p$ts >= p$te | p$te > p$tlen | p$matches > p$alen | p$mapq > 255)) stop('Invalid PAF coordinates or values.')
  # A sequence must have the same total length in every row, or the axes are wrong.
  for (pair in list(c('q', 'qlen'), c('t', 'tlen'))) {
    if (any(vapply(split(p[[pair[2]]], p[[pair[1]]]), function(x) length(unique(x)) != 1L, logical(1))))
      stop('Inconsistent sequence lengths in PAF.')
  }
  tag <- function(f, prefix) {
    x <- f[startsWith(f, prefix)]
    if (length(x)) substring(x[1], nchar(prefix) + 1L) else ''
  }
  # cg describes matches and gaps within an alignment (the CIGAR string).
  # tp marks secondary matches: alternative placements of the same query sequence.
  p$cigar <- vapply(fields, tag, character(1), prefix = 'cg:Z:')
  p$secondary <- vapply(fields, tag, character(1), prefix = 'tp:A:') == 'S'
  p
}

# Assign possible SV types from the arrangement of aligned regions. These rules
# cannot distinguish every real SV from repeats or differences in assembly orientation.
classify <- function(p, min_sv) {
  p$type <- 'COLLINEAR'
  # An anchor is an alignment selected for comparing neighboring query regions.
  p$anchor <- FALSE
  homes <- setNames(character(length(unique(p$q))), unique(p$q))
  for (q in unique(p$q)) {
    ids <- which(p$q == q)
    primary <- ids[!p$secondary[ids]]
    votes <- if (length(primary)) primary else ids
    support <- tapply(p$matches[votes], p$t[votes], sum)
    # Prefer an identically named reference sequence; otherwise use the one with
    # the most matching bases. Different naming can make this assignment uncertain.
    homes[q] <- if (q %in% p$t) q else names(which.max(support))
    # Consider the strongest matches first. Keep a region only if its overlap
    # with each selected query region is less than half the shorter region.
    # Alternative matches stay visible in gray instead of being counted as SVs.
    accepted <- integer()
    for (i in primary[order(-p$matches[primary])]) {
      overlap <- pmax(0, pmin(p$qe[i], p$qe[accepted]) - pmax(p$qs[i], p$qs[accepted]))
      if (!length(accepted) || all(overlap < 0.5 * pmin(p$qe[i] - p$qs[i], p$qe[accepted] - p$qs[accepted])))
        accepted <- c(accepted, i)
    }
    p$anchor[accepted] <- TRUE
    p$type[setdiff(ids, accepted)] <- 'OTHER'
    # Reverse matches suggest INV. A match to another reference sequence is
    # labeled TRA even when reversed, so each aligned region has one main color.
    p$type[accepted[p$strand[accepted] == '-']] <- 'INV'
    p$type[accepted[p$t[accepted] != homes[q]]] <- 'TRA'
    ordered <- accepted[order(p$qs[accepted])]
    for (k in seq_along(ordered)) {
      i <- ordered[k]
      if (p$type[i] != 'COLLINEAR' || k == 1L) next
      earlier <- ordered[seq_len(k - 1L)]
      earlier <- earlier[p$t[earlier] == p$t[i] & p$strand[earlier] == '+' & p$qe[earlier] <= p$qs[i] + min_sv]
      # Separate query regions matching much of the same reference interval suggest
      # DUP. A backward jump without overlap suggests movement within a chromosome.
      ov <- pmax(0, pmin(p$te[i], p$te[earlier]) - pmax(p$ts[i], p$ts[earlier]))
      if (any(ov >= min_sv & ov >= 0.5 * pmin(p$te[i] - p$ts[i], p$te[earlier] - p$ts[earlier]))) {
        p$type[i] <- 'DUP'
      } else if (length(earlier) && p$te[i] < p$ts[tail(earlier, 1)] - min_sv) p$type[i] <- 'TRA'
    }
  }
  p
}

# Turn aligned regions into line endpoints, splitting lines at sufficiently large gaps.
make_segments <- function(p, min_sv) {
  out <- list(); n <- 0L
  add <- function(i, x0, y0, x1, y1, type, gap = FALSE) {
    n <<- n + 1L
    out[[n]] <<- data.frame(t = p$t[i], q = p$q[i], x0, y0, x1, y1, type, gap)
  }
  for (i in seq_len(nrow(p))) {
    # Reference positions always increase. For reverse matches, query positions
    # decrease from the query end; this produces the downward slope of an inversion.
    direction <- if (p$strand[i] == '+') 1 else -1
    x <- p$ts[i]; y <- if (direction == 1) p$qs[i] else p$qe[i]
    if (!nzchar(p$cigar[i])) {
      add(i, x, y, p$te[i], if (direction == 1) p$qe[i] else p$qs[i], p$type[i]); next
    }
    # CIGAR examples: 500M = 500 aligned bases, 100I = 100 inserted query bases,
    # 100D = 100 deleted reference bases. = and X also advance both positions;
    # N skips reference bases and is left ambiguous rather than called a deletion.
    tokens <- regmatches(p$cigar[i], gregexpr('[0-9]+[MIDN=X]', p$cigar[i]))[[1]]
    if (!length(tokens) || paste0(tokens, collapse = '') != p$cigar[i]) stop('Unsupported CIGAR at PAF row ', i)
    lens <- as.numeric(sub('[A-Z=]$', '', tokens)); ops <- sub('^[0-9]+', '', tokens)
    # Stop if these lengths disagree with the reported endpoints; continuing
    # would place the colored gaps at incorrect positions.
    if (any(lens <= 0) || sum(lens[ops %in% c('M', '=', 'X', 'D', 'N')]) != p$te[i] - p$ts[i] ||
        sum(lens[ops %in% c('M', '=', 'X', 'I')]) != p$qe[i] - p$qs[i]) stop('CIGAR span disagrees with PAF row ', i)
    sx <- x; sy <- y
    for (j in seq_along(ops)) {
      op <- ops[j]; len <- lens[j]
      nx <- x + if (op != 'I') len else 0
      ny <- y + if (!op %in% c('D', 'N')) direction * len else 0
      # Small gaps remain within the overall line. For a large gap, end the
      # preceding line and draw a separate dashed segment in the event color.
      if (op %in% c('I', 'D', 'N') && len >= min_sv) {
        if (x != sx || y != sy) add(i, sx, sy, x, y, p$type[i])
        type <- if (p$anchor[i]) switch(op, I = 'INS', D = 'DEL', N = 'OTHER') else 'OTHER'
        add(i, x, y, nx, ny, type, TRUE)
        sx <- nx; sy <- ny
      }
      x <- nx; y <- ny
    }
    if (x != sx || y != sy) add(i, sx, sy, x, y, p$type[i])
  }
  # Also look for gaps between neighboring selected regions. Compare only forward
  # regions in the expected order on the same reference sequence, so inversions
  # and translocations are not accidentally connected and labeled as INS or DEL.
  for (q in unique(p$q)) {
    ids <- which(p$q == q & p$anchor)
    ids <- ids[order(p$qs[ids])]
    if (length(ids) < 2L) next
    for (k in seq_len(length(ids) - 1L)) {
      a <- ids[k]; b <- ids[k + 1L]
      if (p$t[a] != p$t[b] || any(p$type[c(a, b)] != 'COLLINEAR')) next
      qgap <- p$qs[b] - p$qe[a]; tgap <- p$ts[b] - p$te[a]
      if (min(qgap, tgap) < 0 || max(qgap, tgap) < min_sv) next
      # More unaligned query bases suggest INS; more reference bases suggest DEL.
      # Similar gaps on both sides remain ambiguous.
      type <- if (qgap - tgap >= min_sv) 'INS' else if (tgap - qgap >= min_sv) 'DEL' else 'OTHER'
      add(a, p$te[a], p$qe[a], p$ts[b], p$qs[b], type, TRUE)
    }
  }
  do.call(rbind, out)
}

# Draw the same classified regions in each requested format. Styling below does
# not change which regions were selected or how their possible SV types were assigned.
plot_dotplot <- function(p, s, opt, all_p = p) {
  if (!requireNamespace('ggplot2', quietly = TRUE))
    stop('ggplot2 is required. Run with: pixi run Rscript fig3-dotplot/dotplot.R ...')
  library(ggplot2)
  palette <- c(COLLINEAR = '#000000', INS = '#009E73', DEL = '#E41A1C',
               DUP = '#E69F00', INV = '#0072B2', TRA = '#8E44AD', OTHER = '#999999')
  # Place sequences end to end on each axis, keeping shared names in the same order.
  # Use lengths from before filtering so removing weak matches does not shift positions.
  targets <- sort(unique(all_p$t))
  queries <- c(targets[targets %in% all_p$q], sort(setdiff(unique(all_p$q), targets)))
  lengths_for <- function(ids, name, len) setNames(all_p[[len]][match(ids, all_p[[name]])], ids)
  tl <- lengths_for(targets, 't', 'tlen'); ql <- lengths_for(queries, 'q', 'qlen')
  to <- setNames(c(0, head(cumsum(tl), -1)), targets)
  qo <- setNames(c(0, head(cumsum(ql), -1)), queries)
  # Add the preceding sequence lengths, then convert base positions to megabases.
  s$x0 <- (s$x0 + to[s$t]) / 1e6; s$x1 <- (s$x1 + to[s$t]) / 1e6
  s$y0 <- (s$y0 + qo[s$q]) / 1e6; s$y1 <- (s$y1 + qo[s$q]) / 1e6
  # Draw ambiguous and normal matches first so candidate events remain visible.
  s <- s[order(match(s$type, c('OTHER', 'COLLINEAR', 'INS', 'DEL', 'DUP', 'INV', 'TRA'))), ]
  g <- ggplot(s, aes(x = x0, y = y0, colour = type))
  if (opt$display_contigs) {
    g <- g + geom_vline(xintercept = c(to, sum(tl)) / 1e6, colour = 'grey90', linewidth = 0.3) +
      geom_hline(yintercept = c(qo, sum(ql)) / 1e6, colour = 'grey90', linewidth = 0.3)
  }
  g <- g +
    geom_segment(aes(xend = x1, yend = y1, linetype = gap, linewidth = type != 'COLLINEAR')) +
    geom_point(data = s[s$gap, ], aes(x = (x0 + x1) / 2, y = (y0 + y1) / 2), size = 0.7, show.legend = FALSE) +
    scale_colour_manual(values = palette, limits = names(palette),
                        labels = c('Collinear', 'INS', 'DEL', 'DUP', 'INV', 'TRA', 'Other / ambiguous'), drop = FALSE) +
    scale_linetype_manual(values = c('FALSE' = 'solid', 'TRUE' = 'dashed'), guide = 'none') +
    # Line widths: 0.6 for black matches and 1.0 for colored or ambiguous regions.
    scale_linewidth_manual(values = c('FALSE' = 0.6, 'TRUE' = 1.0), guide = 'none') +
    scale_x_continuous(limits = c(0, sum(tl)) / 1e6, expand = expansion(mult = 0.01),
      sec.axis = if (opt$display_contigs) dup_axis(breaks = (to + tl / 2) / 1e6, labels = targets, name = NULL) else waiver()) +
    scale_y_continuous(limits = c(0, sum(ql)) / 1e6, expand = expansion(mult = 0.01),
      sec.axis = if (opt$display_contigs) dup_axis(breaks = (qo + ql / 2) / 1e6, labels = queries, name = NULL) else waiver()) +
    labs(x = 'Reference position (Mb)', y = 'Query position (Mb)', title = opt$title, colour = NULL) +
    theme_classic(base_size = 12) +
    theme(legend.position = 'bottom', aspect.ratio = 1,
          axis.text.x.top = element_text(angle = 90, hjust = 0, size = 8),
          axis.text.y.right = element_text(size = 8),
          axis.ticks.x.top = element_blank(), axis.ticks.y.right = element_blank(),
          plot.title = element_text(hjust = 0.5), plot.margin = margin(12, 12, 12, 12)) +
    guides(colour = guide_legend(nrow = 2, override.aes = list(linewidth = 1.0, linetype = 'solid')))
  for (fmt in opt$format) {
    path <- paste0(opt$output, '.', fmt)
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    # Use R's built-in image writers, including SVG, without needing extra packages.
    device <- switch(fmt, png = grDevices::png, pdf = grDevices::pdf, svg = grDevices::svg)
    ggsave(path, plot = g, device = device, width = opt$width, height = opt$width,
           units = 'in', dpi = opt$dpi, bg = 'white')
    message('Wrote ', path)
  }
  invisible(g)
}

# Run the steps in order: read options, read alignments, filter, classify, then plot.
main <- function(args = commandArgs(trailingOnly = TRUE)) {
  opt <- parse_args(args)
  if (is.null(opt)) return(invisible(NULL))
  all_p <- read_paf(opt$input)
  # MAPQ is alignment confidence. A value of 255 means unknown, not high confidence;
  # exclude it when a positive minimum confidence is requested.
  p <- all_p[all_p$alen >= opt$min_aln &
             (if (opt$min_mapq == 0) rep(TRUE, nrow(all_p)) else all_p$mapq != 255 & all_p$mapq >= opt$min_mapq), ]
  if (!nrow(p)) stop('No alignments remain after filtering.')
  if (any(!nzchar(p$cigar))) message('Some rows lack cg:Z CIGAR: gaps inside those alignments cannot be shown.')
  p <- classify(p, opt$min_sv)
  s <- make_segments(p, opt$min_sv)
  plot_dotplot(p, s, opt, all_p)
  message('Plotted ', nrow(p), ' alignment blocks. Colors show possible SVs relative to the reference, not confirmed events.')
  invisible(list(alignments = p, segments = s))
}
# Run automatically from the command line, but not when the checks load this file.
# Report failures clearly and return a failure status to the calling command.
if (sys.nframe() == 0L) tryCatch(main(), error = function(e) {
  message('Error: ', conditionMessage(e)); quit(status = 1L)
})

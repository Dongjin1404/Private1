#!/usr/bin/env Rscript
# plot_snp_summary.R — Confident SNPs vs Ambiguous positions from 18_variant_calling
#
# USAGE:  Rscript scripts/plot_snp_summary.R <vardir> [outdir]
# EXAMPLE (BREADTH_MIN=0.30 run):
#   Rscript scripts/plot_snp_summary.R 18_variant_calling_0.3 18_variant_calling_0.3/plots
#
# OUTPUTS (in outdir, default = vardir/plots/):
#   snp_bar.pdf           — stacked bar: Confident vs Ambiguous per sample, coloured by pass/fail
#   snp_scatter.pdf       — scatter: Confident vs Ambiguous, sized by coverage breadth
#   snp_ratio.pdf         — bar: ambiguous fraction (%) per sample — best single QC indicator
#   snp_heterogeneity.pdf — bar: strain-heterogeneous SNPs (ALT fraction band, default 40-80%) per
#                           sample — only written if the Het_SNPs_* column is present in the TSV
#
# REQUIRES: ggplot2, dplyr, tidyr, scales   (install.packages(c("ggplot2","dplyr","tidyr","scales")))

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(scales)
})

# ── Arguments ────────────────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript plot_snp_summary.R <vardir> [outdir]\n",
       "  vardir  = output folder of 18_variant_calling (contains <GENOME>/...variant_summary.tsv)")
}
vardir <- sub("/$", "", args[1])
outdir <- if (length(args) >= 2) args[2] else file.path(vardir, "plots")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ── 1. Read all variant_summary.tsv files ────────────────────────────────────
tsv_files <- list.files(vardir, pattern = "variant_summary\\.tsv$",
                        recursive = TRUE, full.names = TRUE)
if (length(tsv_files) == 0)
  stop("No variant_summary.tsv files found under: ", vardir)

message("Found ", length(tsv_files), " summary file(s):")
message(paste(" ", tsv_files, collapse = "\n"))

dat <- lapply(tsv_files, function(f) {
  d <- read.table(f, header = TRUE, sep = "\t",
                  stringsAsFactors = FALSE, check.names = FALSE)
  # The pass column is named Pass<breadth_pct> (e.g. Pass30, Pass70) — normalise it.
  pass_col <- grep("^Pass[0-9]", names(d), value = TRUE)[1]
  if (!is.na(pass_col)) names(d)[names(d) == pass_col] <- "Pass"
  d
})
dat <- do.call(rbind, dat)
stopifnot("Pass" %in% names(dat))

# ── 2. Derived columns + ordering ────────────────────────────────────────────
dat <- dat %>%
  mutate(
    Pass       = factor(Pass, levels = c("yes", "no")),
    Total_SNP  = Confident_SNPs + Ambiguous_Masked,
    Ambig_frac = ifelse(Total_SNP > 0, 100 * Ambiguous_Masked / Total_SNP, 0),
    # Short label: Dataset prefix + Sample
    Tag = paste0(sub("_genomemapp.*|_SRR|_ERR|_[0-9]+$", "", Dataset), "\n", Sample)
  ) %>%
  # Passing samples first, then sort by Confident_SNPs descending within each group
  arrange(Genome, Pass, desc(Confident_SNPs)) %>%
  mutate(Tag = factor(Tag, levels = unique(Tag)))

n_genomes <- length(unique(dat$Genome))
n_samp    <- nrow(dat)

# Colour palette — blue = confident/pass, red = ambiguous/fail, grey = fail fill
COL_CONF  <- "#2E86AB"   # confident SNPs
COL_AMBIG <- "#E84855"   # ambiguous masked
COL_PASS  <- "#2E86AB"
COL_FAIL  <- "#E84855"

# ── 3. PLOT A — Stacked bar: Confident + Ambiguous per sample ────────────────
# Best for: seeing absolute numbers and which component dominates.
# Faded bars = samples that FAILED the breadth filter (excluded from alignment).

dat_long <- dat %>%
  select(Genome, Tag, Pass, Confident_SNPs, Ambiguous_Masked) %>%
  pivot_longer(cols = c(Confident_SNPs, Ambiguous_Masked),
               names_to = "SNP_type", values_to = "Count") %>%
  mutate(SNP_type = factor(SNP_type,
    levels = c("Confident_SNPs", "Ambiguous_Masked"),
    labels = c("Confident (→ called base)", "Ambiguous (→ masked N)")))

bar_width  <- max(9, n_samp * 0.32)
bar_height <- 4.5 * ceiling(n_genomes / 2)

pA <- ggplot(dat_long,
             aes(x = Tag, y = Count, fill = SNP_type, alpha = Pass)) +
  geom_col(width = 0.78, colour = NA) +
  scale_fill_manual(values = c("Confident (→ called base)" = COL_CONF,
                                "Ambiguous (→ masked N)"   = COL_AMBIG),
                    name = NULL) +
  scale_alpha_manual(values = c("yes" = 1.0, "no" = 0.28),
                     labels  = c("yes" = "PASS breadth", "no" = "FAIL breadth"),
                     name = NULL) +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.06))) +
  facet_wrap(~Genome, scales = "free", ncol = min(n_genomes, 2)) +
  labs(
    title    = "Confident SNPs vs Ambiguous positions per sample",
    subtitle = paste0("Source: ", vardir,
                      "   |   Faded = failed breadth coverage filter (excluded from alignment)"),
    x = NULL,
    y = "Number of positions"
  ) +
  theme_bw(base_size = 10) +
  theme(
    axis.text.x      = element_text(angle = 55, hjust = 1, size = 6.5),
    strip.text       = element_text(face = "bold", size = 9),
    legend.position  = "top",
    legend.text      = element_text(size = 9),
    panel.grid.major.x = element_blank(),
    plot.title       = element_text(face = "bold")
  ) +
  guides(fill  = guide_legend(order = 1),
         alpha = guide_legend(order = 2,
                              override.aes = list(fill = "grey40")))

ggsave(file.path(outdir, "snp_bar.pdf"), pA,
       width = bar_width, height = bar_height, units = "in", limitsize = FALSE)
message("Saved: ", file.path(outdir, "snp_bar.pdf"))

# ── 4. PLOT B — Scatter: Confident vs Ambiguous ──────────────────────────────
# Best for: spotting outlier samples (e.g. one with far more ambiguous than confident → contamination).
# The 1:1 dashed line is key: points above it have MORE ambiguous than confident positions.
# Point size encodes coverage breadth (larger = more of the genome covered).

pB <- ggplot(dat,
             aes(x = Confident_SNPs, y = Ambiguous_Masked,
                 colour = Pass, shape = substr(Dataset, 1, 20),
                 size = CovBreadth_Pct)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              colour = "grey55", linewidth = 0.5) +
  geom_point(alpha = 0.82) +
  scale_colour_manual(values = c("yes" = COL_PASS, "no" = COL_FAIL),
                      labels  = c("yes" = "PASS", "no" = "FAIL"),
                      name = "Breadth filter") +
  scale_size_continuous(name = "Breadth\ncoverage (%)", range = c(1.5, 7)) +
  scale_shape_manual(values = 0:15, name = "Dataset") +
  scale_x_continuous(labels = comma) +
  scale_y_continuous(labels = comma) +
  facet_wrap(~Genome, scales = "free", ncol = min(n_genomes, 2)) +
  labs(
    title    = "Confident SNPs vs Ambiguous positions (scatter)",
    subtitle = "Dashed 1:1 line — points above it have more ambiguous than confident positions",
    x = "Confident SNPs (called base)",
    y = "Ambiguous positions (masked → N)"
  ) +
  theme_bw(base_size = 10) +
  theme(
    strip.text    = element_text(face = "bold", size = 9),
    legend.position = "right",
    plot.title    = element_text(face = "bold")
  )

ggsave(file.path(outdir, "snp_scatter.pdf"), pB,
       width = 7 + 3 * (n_genomes > 2), height = 4.5 * ceiling(n_genomes / 2),
       units = "in", limitsize = FALSE)
message("Saved: ", file.path(outdir, "snp_scatter.pdf"))

# ── 5. PLOT C — Ambiguous fraction bar ───────────────────────────────────────
# The single most useful QC chart: what % of "called" positions were flagged uncertain?
# Samples with high ambiguous fraction (>30-50%) are suspect (damage, contamination, mixed strains).
# A dashed line marks a typical warning threshold.

AMBIG_WARN <- 30   # % — adjust if needed

pC <- ggplot(dat, aes(x = Tag, y = Ambig_frac, fill = Pass)) +
  geom_col(width = 0.78, colour = NA, alpha = 0.9) +
  geom_hline(yintercept = AMBIG_WARN, linetype = "dashed",
             colour = "grey30", linewidth = 0.6) +
  annotate("text", x = -Inf, y = AMBIG_WARN + 1.5,
           label = paste0(AMBIG_WARN, "% warning"), hjust = -0.05,
           size = 2.8, colour = "grey30") +
  scale_fill_manual(values = c("yes" = COL_CONF, "no" = COL_FAIL),
                    labels  = c("yes" = "PASS", "no" = "FAIL"),
                    name = "Breadth filter") +
  scale_y_continuous(labels = function(x) paste0(x, "%"),
                     expand = expansion(mult = c(0, 0.08)),
                     limits = c(0, max(dat$Ambig_frac, AMBIG_WARN) * 1.12)) +
  facet_wrap(~Genome, scales = "free_x", ncol = min(n_genomes, 2)) +
  labs(
    title    = "Ambiguous fraction per sample  [= Ambiguous / (Confident + Ambiguous)]",
    subtitle = paste0("High % → mixed signal (DNA damage, contamination, or strain mixture)  ",
                      "  |  Dashed line = ", AMBIG_WARN, "% warning threshold"),
    x = NULL,
    y = "Ambiguous fraction (%)"
  ) +
  theme_bw(base_size = 10) +
  theme(
    axis.text.x      = element_text(angle = 55, hjust = 1, size = 6.5),
    strip.text       = element_text(face = "bold", size = 9),
    legend.position  = "top",
    panel.grid.major.x = element_blank(),
    plot.title       = element_text(face = "bold")
  )

ggsave(file.path(outdir, "snp_ratio.pdf"), pC,
       width = bar_width, height = bar_height, units = "in", limitsize = FALSE)
message("Saved: ", file.path(outdir, "snp_ratio.pdf"))

# ── 5b. PLOT D — Strain-heterogeneous SNPs (allele-fraction band) ─────────────
# Absolute count of SNPs whose ALT fraction fell in the configured band (default 40-80%) — the
# positions most diagnostic of within-sample strain mixture, removed (masked → N) from the alignment.
# Only drawn if 18_variant_calling.sh wrote a Het_SNPs_<lo>_<hi> column (older TSVs lack it).
het_col <- grep("^Het_SNPs_", names(dat), value = TRUE)[1]
if (!is.na(het_col)) {
  band_lbl <- gsub("_", "-", sub("^Het_SNPs_", "", het_col))   # "40_80" -> "40-80"
  dat$Het  <- suppressWarnings(as.numeric(dat[[het_col]]))
  pD <- ggplot(dat, aes(x = Tag, y = Het, fill = Pass)) +
    geom_col(width = 0.78, colour = NA, alpha = 0.9) +
    scale_fill_manual(values = c("yes" = COL_CONF, "no" = COL_FAIL),
                      labels  = c("yes" = "PASS", "no" = "FAIL"),
                      name = "Breadth filter") +
    scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.06))) +
    facet_wrap(~Genome, scales = "free", ncol = min(n_genomes, 2)) +
    labs(
      title    = paste0("Strain-heterogeneous SNPs per sample  (ALT fraction ", band_lbl, "%)"),
      subtitle = "Mixed-allele positions removed (masked \u2192 N) as strain heterogeneity \u2014 higher = more within-sample strain mixture",
      x = NULL,
      y = "Number of heterogeneous positions"
    ) +
    theme_bw(base_size = 10) +
    theme(
      axis.text.x        = element_text(angle = 55, hjust = 1, size = 6.5),
      strip.text         = element_text(face = "bold", size = 9),
      legend.position    = "top",
      panel.grid.major.x = element_blank(),
      plot.title         = element_text(face = "bold")
    )
  ggsave(file.path(outdir, "snp_heterogeneity.pdf"), pD,
         width = bar_width, height = bar_height, units = "in", limitsize = FALSE)
  message("Saved: ", file.path(outdir, "snp_heterogeneity.pdf"))
} else {
  message("Note: no Het_SNPs_* column in the summary \u2014 skipping heterogeneity plot ",
          "(re-run 18_variant_calling.sh to populate it).")
}

# ── 6. Console summary table ──────────────────────────────────────────────────
cat("\n=== SNP Summary Table ===\n")
het_col <- grep("^Het_SNPs_", names(dat), value = TRUE)[1]
base_cols <- c("Genome", "Dataset", "Sample", "Pass",
               "CovBreadth_Pct", "Confident_SNPs", "Ambiguous_Masked",
               if (!is.na(het_col)) het_col, "Ambig_frac")
print(
  dat %>%
    select(all_of(base_cols)) %>%
    rename(`Breadth%` = CovBreadth_Pct,
           `Confident` = Confident_SNPs,
           `Ambiguous` = Ambiguous_Masked,
           `Ambig%` = Ambig_frac) %>%
    mutate(`Ambig%` = round(`Ambig%`, 1)) %>%
    arrange(Genome, Pass, desc(Confident)),
  n = Inf
)
cat("\nOutputs written to: ", outdir, "\n")

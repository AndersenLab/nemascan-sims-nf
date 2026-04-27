#!/usr/bin/env Rscript
# summarize_workdir.R — Analyze Nextflow work directory metrics
#
# Usage:
#   Rscript summarize_workdir.R <workdir_metrics.tsv> [output_dir]
#
# Reads the TSV produced by scan_workdir.sh and produces:
#   - process_summary.tsv       per-process size + symlink statistics
#   - panel_size_summary.tsv    size breakdown by process × panel size
#   - plot_total_usage.pdf      total GB per process (bar chart)
#   - plot_size_distribution.pdf per-task size distributions (violin + box)
#   - plot_panel_scaling.pdf    mean task size vs n_strains for GWAS processes
#   - plot_symlink_fraction.pdf fraction of files that are symlinks per process
#   - plot_largest_file_type.pdf dominant file extension per process
#   - plot_total_by_panel.pdf   total GB stacked by process + panel size

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
  library(readr)
  library(stringr)
  library(forcats)
})

args    <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: summarize_workdir.R <workdir_metrics.tsv> [output_dir]")
tsv_in  <- args[1]
out_dir <- if (length(args) >= 2) args[2] else "."
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cat("Reading:", tsv_in, "\n")

# ─── Pipeline phase lookup ──────────────────────────────────────────────────
PHASES <- tribble(
  ~process,                              ~phase,
  "LOCAL_GET_CONTIG_INFO",               "1_PrepMarkers",
  "BCFTOOLS_EXTRACT_STRAINS",            "1_PrepMarkers",
  "BCFTOOLS_RENAME_CHROMS",              "1_PrepMarkers",
  "PLINK_RECODE_MS_VCF",                 "1_PrepMarkers",
  "PLINK_RECODE_CV_VCF",                 "1_PrepMarkers",
  "BCFTOOLS_CREATE_GENOTYPE_MATRIX",     "1_PrepMarkers",
  "R_FIND_GENOTYPE_MATRIX_EIGEN",        "1_PrepMarkers",
  "LOCAL_COMPILE_EIGENS",                "1_PrepMarkers",
  "PYTHON_SIMULATE_EFFECTS_GLOBAL",      "2_SimPheno",
  "GCTA_SIMULATE_PHENOTYPES",            "2_SimPheno",
  "PLINK_UPDATE_BY_H2",                  "2_SimPheno",
  "GCTA_MAKE_GRM",                       "3_GWAS",
  "GCTA_PERFORM_GWA",                    "3_GWAS",
  "DB_MIGRATION_WRITE_MARKER_SET",       "4_DB",
  "DB_MIGRATION_WRITE_GENOTYPE_MATRIX",  "4_DB",
  "DB_MIGRATION_WRITE_TRAIT_DATA",       "4_DB",
  "DB_MIGRATION_WRITE_GWA_TO_DB",        "4_DB",
  "DB_MIGRATION_AGGREGATE_METADATA",     "4_DB",
  "DB_MIGRATION_ANALYZE_QTL",            "4_DB",
  "DB_MIGRATION_ASSESS_SIMS",            "4_DB"
)

PHASE_COLORS <- c(
  "1_PrepMarkers" = "#4E79A7",
  "2_SimPheno"    = "#F28E2B",
  "3_GWAS"        = "#E15759",
  "4_DB"          = "#76B7B2",
  "5_Unknown"     = "#B07AA1"
)

# ─── Read + clean ────────────────────────────────────────────────────────────
raw <- read_tsv(tsv_in, col_types = cols(
  work_hash          = col_character(),
  process            = col_character(),
  tag                = col_character(),
  panel_size         = col_character(),
  total_bytes        = col_double(),
  n_files            = col_integer(),
  n_symlinks         = col_integer(),
  n_real_files       = col_integer(),
  largest_file_bytes = col_double(),
  largest_file_name  = col_character()
), na = c("", "NA"))

df <- raw |>
  mutate(
    process      = str_remove(process, "^nf-"),
    process      = str_trim(process),
    process      = if_else(is.na(process) | process == "", "UNKNOWN", process),
    panel_size_n = suppressWarnings(as.integer(panel_size)),
    symlink_frac = n_symlinks / pmax(n_files, 1L),
    total_gb     = total_bytes / 1e9,
    largest_ext  = str_extract(largest_file_name, "\\.[^./]+$"),
    largest_ext  = replace_na(largest_ext, "(none)")
  ) |>
  left_join(PHASES, by = "process") |>
  mutate(phase = replace_na(phase, "5_Unknown"))

cat("Task dirs loaded:", nrow(df), "\n")
cat("Processes found:", n_distinct(df$process), "\n")
cat("Panel sizes found:", paste(sort(na.omit(unique(df$panel_size_n))), collapse = ", "), "\n")
cat("Total data volume:", round(sum(df$total_gb, na.rm = TRUE), 1), "GB\n\n")

# ─── Summary table 1: per-process ────────────────────────────────────────────
process_summary <- df |>
  group_by(process, phase) |>
  summarise(
    n_tasks          = n(),
    median_gb        = median(total_gb, na.rm = TRUE),
    mean_gb          = mean(total_gb, na.rm = TRUE),
    p95_gb           = quantile(total_gb, 0.95, na.rm = TRUE),
    total_gb_sum     = sum(total_gb, na.rm = TRUE),
    median_symlink_f = median(symlink_frac, na.rm = TRUE),
    top_ext          = names(sort(table(largest_ext), decreasing = TRUE))[1],
    .groups          = "drop"
  ) |>
  arrange(desc(total_gb_sum))

write_tsv(process_summary, file.path(out_dir, "process_summary.tsv"))
cat("process_summary.tsv written\n")

# ─── Summary table 2: panel-size breakdown ───────────────────────────────────
key_procs <- c(
  "PLINK_RECODE_MS_VCF", "PLINK_RECODE_CV_VCF",
  "BCFTOOLS_EXTRACT_STRAINS", "BCFTOOLS_RENAME_CHROMS",
  "BCFTOOLS_CREATE_GENOTYPE_MATRIX", "R_FIND_GENOTYPE_MATRIX_EIGEN",
  "GCTA_MAKE_GRM", "GCTA_PERFORM_GWA"
)

panel_summary <- df |>
  filter(process %in% key_procs, !is.na(panel_size_n)) |>
  group_by(process, panel_size_n) |>
  summarise(
    n_tasks   = n(),
    median_gb = median(total_gb, na.rm = TRUE),
    mean_gb   = mean(total_gb, na.rm = TRUE),
    p95_gb    = quantile(total_gb, 0.95, na.rm = TRUE),
    total_gb  = sum(total_gb, na.rm = TRUE),
    .groups   = "drop"
  ) |>
  arrange(process, panel_size_n)

write_tsv(panel_summary, file.path(out_dir, "panel_size_summary.tsv"))
cat("panel_size_summary.tsv written\n")

# ─── Plotting helpers ────────────────────────────────────────────────────────
save_plot <- function(p, name, width = 10, height = 7) {
  path_pdf <- file.path(out_dir, paste0(name, ".pdf"))
  path_png <- file.path(out_dir, paste0(name, ".png"))
  ggsave(path_pdf, p, width = width, height = height)
  ggsave(path_png, p, width = width, height = height, dpi = 150)
  cat(basename(path_pdf), "written\n")
}

top10 <- process_summary |> slice_max(total_gb_sum, n = 10) |> pull(process)

# ─── Plot 1: Total GB per process (bar chart) ────────────────────────────────
p1 <- process_summary |>
  filter(total_gb_sum > 0) |>
  mutate(process = fct_reorder(process, total_gb_sum)) |>
  ggplot(aes(x = process, y = total_gb_sum, fill = phase)) +
  geom_col() +
  geom_text(aes(label = paste0(round(total_gb_sum, 1), " GB")),
            hjust = -0.1, size = 3) +
  coord_flip() +
  scale_fill_manual(values = PHASE_COLORS) +
  scale_x_discrete(expand = expansion(mult = c(0.05, 0.05))) +
  scale_y_continuous(
    labels = label_comma(suffix = " GB"),
    expand = expansion(mult = c(0, 0.2))
  ) +
  labs(
    title = "Total work directory disk usage by process",
    x = NULL, y = "Total GB", fill = "Phase"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")
save_plot(p1, "plot_total_usage", width = 10, height = 8)

# ─── Plot 2: Size distribution per process (violin + boxplot) ────────────────
p2 <- df |>
  filter(process %in% top10, total_gb > 0) |>
  mutate(process = fct_reorder(process, total_gb, .fun = median)) |>
  ggplot(aes(x = process, y = total_gb, fill = phase)) +
  geom_violin(scale = "width", trim = TRUE, alpha = 0.7, linewidth = 0.3) +
  geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white",
               alpha = 0.85, linewidth = 0.4) +
  coord_flip() +
  scale_y_log10(
    labels = label_comma(suffix = " GB"),
    breaks = 10^seq(-3, 3)
  ) +
  scale_fill_manual(values = PHASE_COLORS) +
  labs(
    title = "Work directory size distribution (top 10 processes by total usage)",
    subtitle = "Log10 scale; white boxes show median and IQR",
    x = NULL, y = "Task dir size (GB, log10)"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")
save_plot(p2, "plot_size_distribution", width = 11, height = 7)

# ─── Plot 3: Panel size scaling ───────────────────────────────────────────────
gwas_procs <- c("PLINK_RECODE_MS_VCF", "PLINK_RECODE_CV_VCF",
                "GCTA_MAKE_GRM", "GCTA_PERFORM_GWA",
                "BCFTOOLS_EXTRACT_STRAINS", "BCFTOOLS_RENAME_CHROMS")

p3_data <- panel_summary |>
  filter(process %in% gwas_procs, !is.na(panel_size_n))

if (nrow(p3_data) > 0) {
  p3 <- p3_data |>
    ggplot(aes(x = panel_size_n, y = mean_gb,
               color = process, group = process)) +
    geom_ribbon(aes(ymin = median_gb, ymax = p95_gb, fill = process),
                alpha = 0.15, color = NA) +
    geom_line(linewidth = 0.8) +
    geom_point(aes(size = n_tasks)) +
    scale_x_continuous(breaks = c(100, 200, 300, 400, 500),
                       labels = function(x) paste0(x, " strains")) +
    scale_size_continuous(range = c(2, 6), name = "N tasks") +
    scale_color_brewer(palette = "Dark2", name = "Process") +
    scale_fill_brewer(palette = "Dark2", guide = "none") +
    labs(
      title = "Work directory size scaling with panel size",
      subtitle = "Line = mean; ribbon = median–p95",
      x = "Panel size (n strains)", y = "Task dir size (GB)"
    ) +
    theme_bw(base_size = 11) +
    theme(legend.position = "right")
  save_plot(p3, "plot_panel_scaling", width = 10, height = 6)
}

# ─── Plot 4: Symlink fraction per process ─────────────────────────────────────
p4 <- df |>
  filter(process %in% top10, n_files > 0) |>
  mutate(process = fct_reorder(process, symlink_frac, .fun = median, .desc = TRUE)) |>
  ggplot(aes(x = process, y = symlink_frac, fill = phase)) +
  geom_violin(scale = "width", trim = TRUE, alpha = 0.7, linewidth = 0.3) +
  geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white",
               alpha = 0.85, linewidth = 0.4) +
  coord_flip() +
  scale_y_continuous(labels = label_percent(), limits = c(0, 1)) +
  scale_fill_manual(values = PHASE_COLORS) +
  labs(
    title = "Symlink fraction per process",
    subtitle = "Higher = more staged inputs (symlinks), less redundant copying",
    x = NULL, y = "Fraction of files that are symbolic links"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")
save_plot(p4, "plot_symlink_fraction", width = 11, height = 7)

# ─── Plot 5: Largest file type distribution per process ───────────────────────
p5_data <- df |>
  filter(process %in% top10, !is.na(largest_ext)) |>
  count(process, phase, largest_ext) |>
  group_by(process) |>
  mutate(pct = n / sum(n)) |>
  ungroup() |>
  mutate(process = fct_reorder(process, pct, .fun = first, .desc = TRUE))

if (nrow(p5_data) > 0) {
  top_exts <- p5_data |> count(largest_ext, wt = n, sort = TRUE) |> slice_head(n = 10) |> pull(largest_ext)
  p5_data <- p5_data |>
    mutate(largest_ext = if_else(largest_ext %in% top_exts, largest_ext, "other"))

  p5 <- p5_data |>
    ggplot(aes(x = process, y = pct, fill = largest_ext)) +
    geom_col(position = "stack") +
    coord_flip() +
    scale_y_continuous(labels = label_percent()) +
    scale_fill_brewer(palette = "Paired", name = "File extension") +
    labs(
      title = "Dominant file type (by extension) per process",
      subtitle = "Each bar = fraction of tasks where this extension was the largest file",
      x = NULL, y = "Fraction of tasks"
    ) +
    theme_bw(base_size = 11) +
    theme(legend.position = "right")
  save_plot(p5, "plot_largest_file_type", width = 11, height = 7)
}

# ─── Plot 6: Total GB stacked by process + panel size ─────────────────────────
p6_data <- panel_summary |>
  filter(process %in% c("PLINK_RECODE_MS_VCF", "PLINK_RECODE_CV_VCF",
                        "BCFTOOLS_EXTRACT_STRAINS", "BCFTOOLS_RENAME_CHROMS",
                        "BCFTOOLS_CREATE_GENOTYPE_MATRIX",
                        "R_FIND_GENOTYPE_MATRIX_EIGEN",
                        "GCTA_MAKE_GRM", "GCTA_PERFORM_GWA")) |>
  mutate(panel_size_f = factor(panel_size_n))

if (nrow(p6_data) > 0) {
  p6 <- p6_data |>
    mutate(process = fct_reorder(process, total_gb, .fun = sum, .desc = FALSE)) |>
    ggplot(aes(x = process, y = total_gb, fill = panel_size_f)) +
    geom_col(position = "stack") +
    coord_flip() +
    scale_fill_viridis_d(option = "C", direction = -1, name = "Panel size\n(n strains)") +
    scale_y_continuous(labels = label_comma(suffix = " GB")) +
    labs(
      title = "Total storage contribution by process and panel size",
      x = NULL, y = "Total GB"
    ) +
    theme_bw(base_size = 11) +
    theme(legend.position = "right")
  save_plot(p6, "plot_total_by_panel", width = 11, height = 6)
}

cat("\n──────────────────────────────────────────────────────\n")
cat("Done. All outputs written to:", out_dir, "\n")
cat("\nTop 10 processes by total disk usage:\n")
process_summary |>
  slice_head(n = 10) |>
  mutate(total_gb_sum = round(total_gb_sum, 2),
         median_gb    = round(median_gb, 3)) |>
  select(process, phase, n_tasks, median_gb, total_gb_sum, top_ext) |>
  as.data.frame() |>
  print(row.names = FALSE)

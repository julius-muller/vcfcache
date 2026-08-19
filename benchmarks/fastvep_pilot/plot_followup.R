#!/usr/bin/env Rscript

local_library <- Sys.getenv(
  "VCFCACHE_R_LIBRARY",
  unset = file.path(
    path.expand("~"), ".local", "lib", "R",
    paste(R.version$major, strsplit(R.version$minor, "\\.")[[1]][1], sep = ".")
  )
)
.libPaths(c(local_library, .libPaths()))

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: plot_followup.R CORE_SCALING_TSV HEAVY_TSV OUTPUT_DIR")
}

scaling_path <- args[[1]]
heavy_path <- args[[2]]
output_dir <- args[[3]]
panel_dir <- file.path(output_dir, "panels")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(panel_dir, recursive = TRUE, showWarnings = FALSE)

theme_vcfcache <- function() {
  theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", colour = "#263238"),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      legend.position = "bottom",
      legend.title = element_blank()
    )
}

save_formats <- function(plot, directory, stem, width, height) {
  prefix <- file.path(directory, stem)
  ggsave(
    paste0(prefix, ".png"), plot,
    width = width, height = height, dpi = 320, bg = "white"
  )
  ggsave(
    paste0(prefix, ".pdf"), plot,
    width = width, height = height, device = cairo_pdf, bg = "white"
  )
  ggsave(
    paste0(prefix, ".svg"), plot,
    width = width, height = height, device = grDevices::svg, bg = "white"
  )
}

method_values <- c("Direct fastVEP" = "#546E7A", "VCFcache (90% hits)" = "#00897B")

scaling <- read.delim(scaling_path, stringsAsFactors = FALSE, check.names = FALSE)
if (!("configured_thread_limit" %in% names(scaling))) {
  scaling$configured_thread_limit <- scaling$cores
}
if (!("speedup_vs_direct_same_threads" %in% names(scaling))) {
  scaling$speedup_vs_direct_same_threads <- scaling$speedup_vs_direct_same_cores
}
scaling$method <- ifelse(
  scaling$mode == "direct", "Direct fastVEP", "VCFcache (90% hits)"
)
scaling$method <- factor(scaling$method, levels = names(method_values))
scaling$minutes <- scaling$wall_seconds / 60
scaling$core_label <- factor(
  scaling$configured_thread_limit,
  levels = sort(unique(scaling$configured_thread_limit))
)

runtime_mean <- aggregate(minutes ~ core_label + method, scaling, mean)
runtime_min <- aggregate(minutes ~ core_label + method, scaling, min)
runtime_max <- aggregate(minutes ~ core_label + method, scaling, max)
names(runtime_min)[3] <- "ymin"
names(runtime_max)[3] <- "ymax"
runtime_summary <- merge(merge(runtime_mean, runtime_min), runtime_max)

cache_rows <- subset(scaling, mode == "cached")
benefit_mean <- aggregate(
  speedup_vs_direct_same_threads ~ core_label, cache_rows, mean
)
benefit_min <- aggregate(
  speedup_vs_direct_same_threads ~ core_label, cache_rows, min
)
benefit_max <- aggregate(
  speedup_vs_direct_same_threads ~ core_label, cache_rows, max
)
names(benefit_min)[2] <- "ymin"
names(benefit_max)[2] <- "ymax"
benefit_summary <- merge(merge(benefit_mean, benefit_min), benefit_max)

scaling_runtime <- ggplot(
  runtime_summary,
  aes(core_label, minutes, colour = method, group = method)
) +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0.10, linewidth = 0.5) +
  geom_line(linewidth = 1.05) +
  geom_point(size = 3.4) +
  geom_text(
    aes(label = sprintf("%.1f min", minutes)),
    vjust = -0.75,
    size = 3.4,
    show.legend = FALSE
  ) +
  scale_colour_manual(values = method_values) +
  scale_y_continuous(
    limits = c(0, max(scaling$minutes) * 1.12),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(
    x = "Configured thread limit (enforced CPU affinity)",
    y = "End-to-end wall time (minutes)",
    title = "A  Both workflows benefit from more cores"
  ) +
  theme_vcfcache()

scaling_benefit <- ggplot(
  benefit_summary,
  aes(core_label, speedup_vs_direct_same_threads)
) +
  geom_hline(yintercept = 1, colour = "#90A4AE", linewidth = 0.5) +
  geom_col(width = 0.62, fill = "#00897B", colour = "#263238", linewidth = 0.35) +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0.14, linewidth = 0.55) +
  geom_text(
    aes(label = sprintf("%.2fx", speedup_vs_direct_same_threads)),
    vjust = -0.55,
    size = 4,
    fontface = "bold"
  ) +
  scale_y_continuous(
    limits = c(0, max(cache_rows$speedup_vs_direct_same_threads) * 1.18),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(
    x = "Configured thread limit (enforced CPU affinity)",
    y = "Speedup over direct fastVEP",
    title = "B  VCFcache benefit at 90% hits",
    subtitle = "Direct and cached cells use the same process-wide CPU affinity"
  ) +
  theme_vcfcache() +
  theme(legend.position = "none")

scaling_combined <- scaling_runtime + scaling_benefit +
  plot_layout(widths = c(1.18, 0.82)) +
  plot_annotation(
    title = "VCFcache remains faster as fastVEP receives more cores",
    subtitle = "One 4.33M-variant WGS per cell; 90% cache hits; taskset-enforced affinity; range shown for repeated control",
    theme = theme(plot.title = element_text(face = "bold", size = 15))
  )

save_formats(scaling_runtime, panel_dir, "fastvep_core_scaling_runtime", 6.4, 4.5)
save_formats(scaling_benefit, panel_dir, "fastvep_core_scaling_speedup", 5.4, 4.5)
save_formats(scaling_combined, output_dir, "fastvep_core_scaling", 12.0, 5.0)

heavy <- read.delim(heavy_path, stringsAsFactors = FALSE, check.names = FALSE)
heavy$condition_label <- c(
  "direct" = "Direct",
  "hit-090" = "90% hits",
  "hit-100" = "100% hits"
)[heavy$condition]
heavy$condition_label <- factor(
  heavy$condition_label,
  levels = c("Direct", "90% hits", "100% hits")
)
heavy$minutes <- heavy$wall_seconds / 60
heavy$fill_group <- ifelse(heavy$condition == "direct", "Direct fastVEP", "VCFcache")

heavy_runtime <- ggplot(heavy, aes(condition_label, minutes, fill = fill_group)) +
  geom_col(width = 0.64, colour = "#263238", linewidth = 0.35) +
  geom_text(
    aes(label = sprintf("%.1f min", minutes)),
    vjust = -0.5,
    size = 3.8
  ) +
  scale_fill_manual(values = c("Direct fastVEP" = "#546E7A", "VCFcache" = "#00897B")) +
  scale_y_continuous(
    limits = c(0, max(heavy$minutes) * 1.15),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(
    x = NULL,
    y = "End-to-end wall time (minutes)",
    title = "A  Dense supplementary pipeline"
  ) +
  theme_vcfcache() +
  theme(legend.position = "none")

heavy_cached <- subset(heavy, condition != "direct")
heavy_speedup <- ggplot(heavy_cached, aes(condition_label, speedup, fill = condition_label)) +
  geom_hline(yintercept = 1, colour = "#90A4AE", linewidth = 0.5) +
  geom_col(width = 0.60, colour = "#263238", linewidth = 0.35) +
  geom_text(
    aes(label = sprintf("%.2fx", speedup)),
    vjust = -0.55,
    size = 4.1,
    fontface = "bold"
  ) +
  scale_fill_manual(values = c("90% hits" = "#26A69A", "100% hits" = "#00796B")) +
  scale_y_continuous(
    limits = c(0, max(heavy_cached$speedup) * 1.16),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(
    x = NULL,
    y = "Speedup over direct fastVEP",
    title = "B  Acceleration from reuse"
  ) +
  theme_vcfcache() +
  theme(legend.position = "none")

heavy_combined <- heavy_runtime + heavy_speedup +
  plot_layout(widths = c(1, 1)) +
  plot_annotation(
    title = "VCFcache substantially accelerates heavily configured fastVEP",
    subtitle = "Ten supplementary databases; one 4.33M-variant WGS; 16 CPU threads; complete output equality",
    theme = theme(plot.title = element_text(face = "bold", size = 15))
  )

save_formats(heavy_runtime, panel_dir, "fastvep_heavy_runtime", 5.8, 4.5)
save_formats(heavy_speedup, panel_dir, "fastvep_heavy_speedup", 5.8, 4.5)
save_formats(heavy_combined, output_dir, "fastvep_heavy_run", 11.6, 5.0)

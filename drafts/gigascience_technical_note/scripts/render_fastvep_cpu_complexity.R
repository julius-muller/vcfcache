arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) != 1) {
  stop("Usage: Rscript --vanilla scripts/render_fastvep_cpu_complexity.R <draft-directory>")
}

draft_dir <- normalizePath(arguments[[1]], mustWork = TRUE)
repo_dir <- normalizePath(file.path(draft_dir, "..", ".."), mustWork = TRUE)
source(file.path(repo_dir, "benchmarks", "figures", "R", "common.R"))

core_path <- file.path(
  repo_dir, "benchmarks", "fastvep_pilot", "source_data",
  "core_scaling_kpgp00319.tsv"
)
heavy_path <- file.path(
  repo_dir, "benchmarks", "fastvep_pilot", "source_data",
  "heavy_core_scaling_kpgp00319.tsv"
)
if (!file.exists(heavy_path)) {
  stop("The enforced-CPU heavy fastVEP scaling result has not been collected")
}

core <- read_tsv(core_path)
heavy <- read_tsv(heavy_path)
core <- core[core[["repeat"]] == 1, ]
heavy <- heavy[heavy[["repeat"]] == 1, ]
core$workload <- "Core consequences"
heavy$workload <- "Ten supplementary databases"
rows <- rbind(core, heavy)

if (nrow(rows) != 10 || !all(as_logical_semantic(rows$semantic_pass))) {
  stop("Expected 10 semantically validated fastVEP CPU-scaling observations")
}

rows$workload <- factor(
  rows$workload,
  levels = c("Core consequences", "Ten supplementary databases")
)
rows$mode_label <- factor(
  ifelse(rows$mode == "direct", "Direct fastVEP", "VCFcache, 90% hits"),
  levels = c("Direct fastVEP", "VCFcache, 90% hits")
)

ink <- "#172A34"
blue <- "#2B6078"
orange <- "#D97824"
green <- "#3D8E55"
purple <- "#745A9C"
line_grey <- "#D7DEE1"

journal_theme <- function(base_size = 9.5) {
  theme_classic(base_size = base_size, base_family = "Lato") +
    theme(
      text = element_text(color = ink),
      axis.text = element_text(color = ink),
      axis.title = element_text(color = ink),
      axis.line = element_line(color = ink, linewidth = 0.35),
      axis.ticks = element_line(color = ink, linewidth = 0.35),
      panel.grid.major.y = element_line(color = line_grey, linewidth = 0.3),
      panel.grid.minor = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", color = ink),
      legend.position = "bottom",
      legend.title = element_blank(),
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      plot.caption = element_blank(),
      plot.tag = element_text(face = "bold", size = 12, color = ink),
      plot.tag.position = c(0.005, 0.995),
      plot.margin = margin(8, 8, 8, 8)
    )
}

mode_colours <- c("Direct fastVEP" = ink, "VCFcache, 90% hits" = green)
mode_shapes <- c("Direct fastVEP" = 16, "VCFcache, 90% hits" = 15)

p_runtime <- ggplot(
  rows,
  aes(
    x = configured_thread_limit,
    y = wall_seconds / 60,
    colour = mode_label,
    shape = mode_label,
    group = mode_label
  )
) +
  geom_line(linewidth = 0.72) +
  geom_point(size = 2.1, stroke = 0.3) +
  facet_wrap(~workload, nrow = 1) +
  scale_colour_manual(values = mode_colours) +
  scale_shape_manual(values = mode_shapes) +
  scale_x_continuous(breaks = c(1, 10, 32)) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.08))) +
  labs(
    tag = "A", x = "CPUs available to the process",
    y = "Observed runtime per genome (minutes)"
  ) +
  journal_theme()

speedup <- rows[rows$mode == "cached", ]
workload_colours <- c(
  "Core consequences" = blue,
  "Ten supplementary databases" = purple
)
workload_shapes <- c("Core consequences" = 16, "Ten supplementary databases" = 17)

p_speedup <- ggplot(
  speedup,
  aes(
    x = configured_thread_limit,
    y = speedup_vs_direct_same_threads,
    colour = workload,
    shape = workload,
    group = workload
  )
) +
  geom_hline(yintercept = 1, colour = line_grey, linewidth = 0.5) +
  geom_line(linewidth = 0.72) +
  geom_point(size = 2.1, stroke = 0.3) +
  scale_colour_manual(values = workload_colours) +
  scale_shape_manual(values = workload_shapes) +
  scale_x_continuous(breaks = c(1, 10, 32)) +
  scale_y_continuous(
    breaks = scales::breaks_pretty(n = 5),
    limits = c(0.9, NA),
    expand = expansion(mult = c(0, 0.08))
  ) +
  labs(
    tag = "B", x = "CPUs available to the process",
    y = "Speedup over direct fastVEP"
  ) +
  journal_theme()

figure <- p_runtime + p_speedup + plot_layout(widths = c(1.3, 1))
save_plot(
  figure,
  file.path(
    draft_dir, "supplement", "figures",
    "supplementary_figure2_fastvep_cpu_complexity"
  ),
  width = 7.2, height = 4.5, dpi = 600
)
save_plot(
  p_runtime,
  file.path(
    draft_dir, "supplement", "figures", "panels",
    "supplementary_figure2_A_runtime"
  ),
  width = 5.2, height = 4.1, dpi = 600
)
save_plot(
  p_speedup,
  file.path(
    draft_dir, "supplement", "figures", "panels",
    "supplementary_figure2_B_speedup"
  ),
  width = 4.2, height = 4.1, dpi = 600
)

write_tsv(
  rows,
  file.path(
    draft_dir, "figures", "source_data", "fastvep_cpu_complexity.tsv"
  )
)

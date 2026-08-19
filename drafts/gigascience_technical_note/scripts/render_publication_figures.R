arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) != 1) {
  stop("Usage: Rscript --vanilla scripts/render_publication_figures.R <draft-directory>")
}

draft_dir <- normalizePath(arguments[[1]], mustWork = TRUE)
repo_dir <- normalizePath(file.path(draft_dir, "..", ".."), mustWork = TRUE)
source(file.path(repo_dir, "benchmarks", "figures", "R", "common.R"))

metrics_path <- file.path(
  repo_dir, "benchmarks", "figures", "source_data", "final",
  "2026-08-09T223718Z-light-matched-final", "annotator_external_wgs_metrics.tsv"
)
complexity_path <- file.path(
  repo_dir, "benchmarks", "figures", "output", "final",
  "2026-08-12T223200Z-pipeline-complexity-dual-cache-final-v2", "source",
  "dual_cache_pipeline_spectrum.tsv"
)

metrics <- read_tsv(metrics_path)
complexity <- read_tsv(complexity_path)

if (nrow(metrics) != 312 || !all(as_logical_semantic(metrics$semantic_pass))) {
  stop("Expected 312 semantically validated external-WGS observations")
}
if (nrow(complexity) != 12 || !all(as_logical_semantic(complexity$semantic_pass))) {
  stop("Expected 12 semantically validated pipeline-complexity observations")
}

ink <- "#172A34"
blue <- "#2B6078"
orange <- "#D97824"
green <- "#3D8E55"
purple <- "#745A9C"
light_blue <- "#DCEAF0"
light_orange <- "#F7E6D5"
light_green <- "#DDEEDF"
light_purple <- "#E9E3F1"
line_grey <- "#D7DEE1"

journal_theme <- function(base_size = 10) {
  theme_classic(base_size = base_size, base_family = "Lato") +
    theme(
      text = element_text(color = ink),
      axis.text = element_text(color = ink),
      axis.title = element_text(color = ink, face = "plain"),
      axis.line = element_line(color = ink, linewidth = 0.35),
      axis.ticks = element_line(color = ink, linewidth = 0.35),
      panel.grid.major.y = element_line(color = line_grey, linewidth = 0.3),
      panel.grid.minor = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", color = ink),
      legend.position = "bottom",
      legend.title = element_text(face = "plain"),
      legend.key.height = grid::unit(0.36, "cm"),
      legend.key.width = grid::unit(0.55, "cm"),
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      plot.caption = element_blank(),
      plot.margin = margin(8, 8, 8, 8)
    )
}

panel_tag_theme <- theme(
  plot.tag = element_text(
    family = "Lato", face = "bold", size = 12, color = ink
  ),
  plot.tag.position = c(0.005, 0.995)
)

box_layer <- function(xmin, xmax, ymin, ymax, fill, label, size = 3.5) {
  list(
    annotate(
      "rect", xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
      fill = fill, color = ink, linewidth = 0.55
    ),
    annotate(
      "text", x = (xmin + xmax) / 2, y = (ymin + ymax) / 2,
      label = label, family = "Lato", color = ink, size = size,
      lineheight = 0.95
    )
  )
}

arrow_layer <- function(x, xend, y, yend, colour = ink) {
  annotate(
    "segment", x = x, xend = xend, y = y, yend = yend,
    colour = colour, linewidth = 0.65,
    arrow = grid::arrow(length = grid::unit(0.12, "inches"), type = "closed")
  )
}

# Figure 1: workflow, with cache construction separated from per-sample use.
workflow_build <- ggplot() +
  box_layer(0.3, 2.2, 0.7, 1.7, light_orange, "Public variants\nor cohort variants") +
  arrow_layer(2.25, 3.0, 1.2, 1.2) +
  box_layer(3.05, 4.85, 0.7, 1.7, light_purple, "Annotate once\nwith a fixed recipe") +
  arrow_layer(4.9, 5.65, 1.2, 1.2) +
  box_layer(5.7, 7.7, 0.7, 1.7, light_green, "Recipe-specific\nannotation cache") +
  coord_cartesian(xlim = c(0, 8), ylim = c(0.35, 1.95), clip = "off") +
  theme_void(base_family = "Lato") +
  theme(plot.margin = margin(8, 20, 0, 20))

workflow_sample <- ggplot() +
  box_layer(0.15, 1.55, 2.0, 3.0, light_orange, "New\nsample", 3.0) +
  arrow_layer(1.6, 2.15, 2.5, 2.5) +
  box_layer(2.2, 4.0, 2.0, 3.0, light_blue, "Cache lookup\nand split", 3.0) +
  arrow_layer(4.05, 4.75, 2.55, 3.65, green) +
  arrow_layer(4.05, 4.75, 2.45, 1.35, purple) +
  box_layer(4.8, 7.3, 3.2, 4.2, light_green, "Cache hits\ncopy cached\nannotations", 2.7) +
  box_layer(4.8, 7.3, 0.8, 1.8, light_purple, "Cache misses\nannotate remaining\nvariants", 2.7) +
  arrow_layer(7.35, 8.05, 3.65, 2.6, green) +
  arrow_layer(7.35, 8.05, 1.35, 2.4, purple) +
  box_layer(8.1, 9.45, 2.0, 3.0, "#EEF1F2", "Merge\nrecords", 3.0) +
  arrow_layer(9.5, 10.05, 2.5, 2.5) +
  box_layer(10.1, 12.55, 2.0, 3.0, light_green, "Sorted and indexed\nannotated BCF", 2.9) +
  annotate(
    "text", x = 6.05, y = 4.45, label = "reuse", family = "Lato",
    fontface = "bold", size = 3.0, color = green
  ) +
  annotate(
    "text", x = 6.05, y = 0.55, label = "annotate only what is new",
    family = "Lato", fontface = "bold", size = 3.0, color = purple
  ) +
  coord_cartesian(xlim = c(0, 12.75), ylim = c(0.25, 4.5), clip = "off") +
  theme_void(base_family = "Lato") +
  theme(plot.margin = margin(0, 12, 8, 12))

figure1 <- workflow_build / workflow_sample +
  plot_layout(heights = c(0.42, 1)) +
  plot_annotation(tag_levels = "A") & panel_tag_theme

save_plot(
  figure1,
  file.path(draft_dir, "figures", "main", "figure1_workflow_v2"),
  width = 7.2, height = 5.2, dpi = 600
)
save_plot(
  workflow_build + labs(tag = "A") + panel_tag_theme,
  file.path(draft_dir, "figures", "main", "panels", "figure1_A_cache_construction"),
  width = 7.2, height = 2.0, dpi = 600
)
save_plot(
  workflow_sample + labs(tag = "B") + panel_tag_theme,
  file.path(draft_dir, "figures", "main", "panels", "figure1_B_new_sample"),
  width = 7.2, height = 3.5, dpi = 600
)

# Figure 2: every external genome is visible; no source-overlap samples are used.
strategy_labels <- c(
  gnomad_af_0.1 = "gnomAD\nAF ≥ 10%",
  gnomad_af_0.01 = "gnomAD\nAF ≥ 1%",
  cohort_3_genomes = "Three-genome\ncohort"
)
metrics$strategy_label <- factor(
  strategy_labels[metrics$strategy], levels = unname(strategy_labels)
)
metrics$cohort_label <- factor(
  toupper(metrics$cohort), levels = c("KPGP", "SGDP", "PGP")
)

cohort_colours <- c(KPGP = blue, SGDP = orange, PGP = green)
cohort_shapes <- c(KPGP = 16, SGDP = 17, PGP = 15)

per_sample_panel <- function(tool_name, y_breaks) {
  values <- metrics[metrics$tool == tool_name, ]
  summaries <- do.call(
    rbind,
    lapply(seq_along(levels(values$strategy_label)), function(index) {
      strategy <- levels(values$strategy_label)[[index]]
      group <- values[values$strategy_label == strategy, ]
      interval <- stratified_bootstrap_median_interval(
        group$speedup, group$cohort_label,
        seed = 20260817 + index + ifelse(tool_name == "vep", 0, 100)
      )
      data.frame(
        strategy_label = factor(strategy, levels = levels(values$strategy_label)),
        median = median(group$speedup),
        lower = interval[[1]],
        upper = interval[[2]]
      )
    })
  )
  ggplot(values, aes(x = strategy_label, y = speedup)) +
    geom_jitter(
      aes(colour = cohort_label, shape = cohort_label),
      width = 0.14, height = 0, size = 1.65, alpha = 0.78, stroke = 0.25
    ) +
    geom_errorbar(
      data = summaries,
      aes(x = strategy_label, ymin = lower, ymax = upper),
      inherit.aes = FALSE, width = 0.08, colour = ink, linewidth = 0.6
    ) +
    geom_point(
      data = summaries,
      aes(x = strategy_label, y = median),
      inherit.aes = FALSE, shape = 18, size = 2.6, colour = ink
    ) +
    scale_colour_manual(values = cohort_colours, name = "Cohort") +
    scale_shape_manual(values = cohort_shapes, name = "Cohort") +
    scale_y_continuous(breaks = y_breaks, expand = expansion(mult = c(0.02, 0.08))) +
    labs(x = NULL, y = "Speedup over direct annotation") +
    journal_theme(9.5)
}

p2a <- per_sample_panel("vep", seq(0, 14, 2))
p2b <- per_sample_panel("fastvep", seq(1, 2.6, 0.4))
figure2 <- p2a + p2b +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A") &
  panel_tag_theme &
  theme(legend.position = "bottom")

save_plot(
  figure2,
  file.path(draft_dir, "figures", "main", "figure2_real_world_annotators_v2"),
  width = 7.2, height = 4.4, dpi = 600
)
save_plot(
  p2a + labs(tag = "A") + panel_tag_theme,
  file.path(draft_dir, "figures", "main", "panels", "figure2_A_vep_speedup"),
  width = 4.2, height = 4.2, dpi = 600
)
save_plot(
  p2b + labs(tag = "B") + panel_tag_theme,
  file.path(draft_dir, "figures", "main", "panels", "figure2_B_fastvep_speedup"),
  width = 4.2, height = 4.2, dpi = 600
)

# Figure 3: measured absolute runtime and speedup across controlled VEP loads.
complexity$cache_strategy <- factor(
  complexity$cache_strategy,
  levels = c("AF ≥ 10% cache", "AF ≥ 1% cache")
)
levels(complexity$cache_strategy) <- c("gnomAD AF ≥ 10%", "gnomAD AF ≥ 1%")

direct <- unique(complexity[c("delay_ms", "direct_hours")])
direct$condition <- "Direct annotation"
direct$runtime_hours <- direct$direct_hours
direct$direct_runtime_hours <- direct$direct_hours

cached <- data.frame(
  delay_ms = complexity$delay_ms,
  direct_hours = complexity$direct_hours,
  condition = as.character(complexity$cache_strategy),
  runtime_hours = complexity$cached_hours,
  direct_runtime_hours = complexity$direct_hours
)
runtime <- rbind(
  data.frame(
    delay_ms = direct$delay_ms,
    direct_hours = direct$direct_hours,
    condition = direct$condition,
    runtime_hours = direct$runtime_hours,
    direct_runtime_hours = direct$direct_runtime_hours
  ),
  cached
)
runtime$condition <- factor(
  runtime$condition,
  levels = c("Direct annotation", "gnomAD AF ≥ 10%", "gnomAD AF ≥ 1%")
)

runtime_colours <- c(
  "Direct annotation" = ink,
  "gnomAD AF ≥ 10%" = orange,
  "gnomAD AF ≥ 1%" = green
)
runtime_shapes <- c(
  "Direct annotation" = 16,
  "gnomAD AF ≥ 10%" = 17,
  "gnomAD AF ≥ 1%" = 15
)

p3a <- ggplot(
  runtime,
  aes(x = direct_runtime_hours, y = runtime_hours, colour = condition, shape = condition)
) +
  geom_line(aes(group = condition), linewidth = 0.75) +
  geom_point(size = 2.2, stroke = 0.3) +
  scale_colour_manual(values = runtime_colours, name = NULL) +
  scale_shape_manual(values = runtime_shapes, name = NULL) +
  scale_x_continuous(breaks = c(5, 10, 15, 20), limits = c(3.5, 23.5)) +
  scale_y_continuous(breaks = c(0, 5, 10, 15, 20), limits = c(0, 23.5)) +
  labs(x = "Direct pipeline runtime (hours)", y = "Observed runtime per genome (hours)") +
  journal_theme(9.5)

p3b <- ggplot(
  complexity,
  aes(x = direct_hours, y = speedup, colour = cache_strategy, shape = cache_strategy)
) +
  geom_line(aes(group = cache_strategy), linewidth = 0.75) +
  geom_point(size = 2.2, stroke = 0.3) +
  scale_colour_manual(
    values = runtime_colours, limits = names(runtime_colours), drop = FALSE,
    name = NULL
  ) +
  scale_shape_manual(
    values = runtime_shapes, limits = names(runtime_shapes), drop = FALSE,
    name = NULL
  ) +
  scale_x_continuous(breaks = c(5, 10, 15, 20), limits = c(3.5, 23.5)) +
  scale_y_continuous(breaks = c(2, 4, 6, 8, 10), limits = c(1, 10.8)) +
  labs(x = "Direct pipeline runtime (hours)", y = "Speedup over direct annotation") +
  journal_theme(9.5)

p3a <- p3a + labs(tag = "A") + panel_tag_theme + theme(legend.position = "none")
p3b <- p3b + labs(tag = "B") + panel_tag_theme + theme(legend.position = "none")

legend_values <- data.frame(
  x = c(0.5, 3.9, 7.3),
  label = names(runtime_colours),
  condition = factor(names(runtime_colours), levels = names(runtime_colours))
)
runtime_legend <- ggplot(
  legend_values, aes(x = x, y = 0, colour = condition, shape = condition)
) +
  geom_segment(
    aes(x = x - 0.3, xend = x + 0.3, y = 0, yend = 0), linewidth = 0.7
  ) +
  geom_point(size = 2.2, stroke = 0.3) +
  geom_text(
    aes(x = x + 0.45, label = label), hjust = 0, colour = ink,
    family = "Lato", size = 3.0, show.legend = FALSE
  ) +
  scale_colour_manual(values = runtime_colours, guide = "none") +
  scale_shape_manual(values = runtime_shapes, guide = "none") +
  coord_cartesian(xlim = c(0, 10.5), ylim = c(-0.5, 0.5), clip = "off") +
  theme_void(base_family = "Lato") +
  theme(plot.margin = margin(0, 5, 0, 5))

figure3 <- (p3a + p3b) / runtime_legend +
  plot_layout(heights = c(1, 0.09))

save_plot(
  figure3,
  file.path(draft_dir, "figures", "main", "figure3_pipeline_complexity_v2"),
  width = 7.2, height = 4.4, dpi = 600
)
save_plot(
  p3a,
  file.path(draft_dir, "figures", "main", "panels", "figure3_A_runtime"),
  width = 4.2, height = 4.0, dpi = 600
)
save_plot(
  p3b,
  file.path(draft_dir, "figures", "main", "panels", "figure3_B_speedup"),
  width = 4.2, height = 4.0, dpi = 600
)

# Supplementary Figure 1: detailed cache coverage and time returned per genome.
metrics$time_saved_percent <- 100 * (1 - metrics$relative_runtime)
metrics$tool_label <- factor(
  ifelse(metrics$tool == "vep", "VEP", "fastVEP"),
  levels = c("VEP", "fastVEP")
)

s1a <- ggplot(metrics, aes(x = strategy_label, y = 100 * cache_hit_rate)) +
  geom_jitter(
    aes(colour = cohort_label, shape = cohort_label),
    width = 0.14, height = 0, size = 1.45, alpha = 0.72, stroke = 0.25
  ) +
  stat_summary(
    fun = median, geom = "point", shape = 95, size = 7.5, colour = ink
  ) +
  scale_colour_manual(values = cohort_colours, name = "Cohort") +
  scale_shape_manual(values = cohort_shapes, name = "Cohort") +
  scale_y_continuous(
    breaks = seq(50, 100, 10), limits = c(48, 100),
    labels = function(x) paste0(x, "%")
  ) +
  labs(x = NULL, y = "Variants found in cache") +
  journal_theme(9.5) +
  theme(legend.position = "none")

s1b <- ggplot(metrics, aes(x = strategy_label, y = time_saved_percent)) +
  geom_jitter(
    aes(colour = cohort_label, shape = cohort_label),
    width = 0.14, height = 0, size = 1.35, alpha = 0.7, stroke = 0.25
  ) +
  stat_summary(
    fun = median, geom = "point", shape = 95, size = 7.5, colour = ink
  ) +
  facet_wrap(~tool_label, nrow = 1) +
  scale_colour_manual(values = cohort_colours, name = "Cohort") +
  scale_shape_manual(values = cohort_shapes, name = "Cohort") +
  scale_y_continuous(
    breaks = seq(20, 100, 20), limits = c(15, 100),
    labels = function(x) paste0(x, "%")
  ) +
  labs(x = NULL, y = "Wall time saved per genome") +
  journal_theme(9.5)

supplementary1 <- s1a / s1b +
  plot_layout(heights = c(1, 1.05)) +
  plot_annotation(tag_levels = "A") &
  panel_tag_theme

save_plot(
  supplementary1,
  file.path(
    draft_dir, "supplement", "figures", "supplementary_figure1_external_wgs_v2"
  ),
  width = 7.2, height = 7.3, dpi = 600
)
save_plot(
  s1a + labs(tag = "A") + panel_tag_theme,
  file.path(
    draft_dir, "supplement", "figures", "panels",
    "supplementary_figure1_A_cache_coverage"
  ),
  width = 7.2, height = 3.5, dpi = 600
)
save_plot(
  s1b + labs(tag = "B") + panel_tag_theme,
  file.path(
    draft_dir, "supplement", "figures", "panels",
    "supplementary_figure1_B_time_saved"
  ),
  width = 7.2, height = 4.0, dpi = 600
)

write_tsv(
  metrics[setdiff(names(metrics), "statistics_mode")],
  file.path(draft_dir, "figures", "source_data", "external_wgs_metrics.tsv")
)
write_tsv(
  complexity[setdiff(names(complexity), "statistics_mode")],
  file.path(draft_dir, "figures", "source_data", "pipeline_complexity.tsv")
)

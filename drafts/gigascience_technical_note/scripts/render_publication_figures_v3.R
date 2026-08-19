arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) != 1) {
  stop("Usage: Rscript --vanilla scripts/render_publication_figures_v3.R <draft-directory>")
}

draft_dir <- normalizePath(arguments[[1]], mustWork = TRUE)
repo_dir <- normalizePath(file.path(draft_dir, "..", ".."), mustWork = TRUE)
source(file.path(repo_dir, "benchmarks", "figures", "R", "common.R"))

complexity_path <- file.path(
  repo_dir, "benchmarks", "figures", "output", "final",
  "2026-08-12T223200Z-pipeline-complexity-dual-cache-final-v2", "source",
  "dual_cache_pipeline_spectrum.tsv"
)
complexity <- read_tsv(complexity_path)
if (nrow(complexity) != 12 || !all(as_logical_semantic(complexity$semantic_pass))) {
  stop("Expected 12 semantically validated pipeline-complexity observations")
}

ink <- "#172A34"
blue <- "#356A82"
orange <- "#C96F2D"
green <- "#3E8153"
purple <- "#6D5790"
grey <- "#65747B"
light_blue <- "#E2EDF2"
light_orange <- "#F7E9DC"
light_green <- "#E1EEE4"
light_purple <- "#EBE6F2"
light_grey <- "#EEF1F2"
line_grey <- "#D7DEE1"

journal_theme <- function(base_size = 10) {
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

box_layer <- function(
  xmin, xmax, ymin, ymax, fill, label, border = ink, size = 2.75
) {
  list(
    annotate(
      "rect", xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
      fill = fill, color = border, linewidth = 0.55
    ),
    annotate(
      "text", x = (xmin + xmax) / 2, y = (ymin + ymax) / 2,
      label = label, family = "Lato", color = ink, size = size,
      lineheight = 0.94
    )
  )
}

arrow_layer <- function(
  x, xend, y, yend, colour = ink, linetype = "solid", linewidth = 0.62
) {
  annotate(
    "segment", x = x, xend = xend, y = y, yend = yend,
    colour = colour, linetype = linetype, linewidth = linewidth,
    arrow = grid::arrow(length = grid::unit(0.105, "inches"), type = "closed")
  )
}

# Figure 1 v3: connected one-time setup and per-sample paths.
figure1_v3_base <- ggplot() +
  annotate(
    "text", x = 0.05, y = 8.45, label = "A", hjust = 0,
    family = "Lato", fontface = "bold", size = 4.2, color = ink
  ) +
  annotate(
    "text", x = 0.05, y = 3.65, label = "B", hjust = 0,
    family = "Lato", fontface = "bold", size = 4.2, color = ink
  ) +
  annotate(
    "segment", x = 0.25, xend = 16.65, y = 3.85, yend = 3.85,
    colour = line_grey, linewidth = 0.5
  ) +
  box_layer(
    0.55, 2.65, 7.15, 8.05, light_blue,
    "Ready-to-Use\nRecipe Cache\nfrom Zenodo", blue
  ) +
  box_layer(
    0.55, 2.65, 5.65, 6.55, light_blue,
    "Blueprint\nfrom Zenodo", blue
  ) +
  box_layer(
    0.55, 2.65, 4.15, 5.05, light_orange,
    "Existing Cohort\nVariants", orange
  ) +
  arrow_layer(2.72, 3.25, 7.60, 7.60, blue) +
  arrow_layer(2.72, 3.25, 6.10, 6.10, blue) +
  arrow_layer(2.72, 3.25, 4.60, 4.60, orange) +
  box_layer(
    3.30, 5.45, 7.15, 8.05, light_grey,
    "Match Build, Tool,\nVersion and Recipe", grey, 2.25
  ) +
  box_layer(
    3.30, 5.45, 5.65, 6.55, light_purple,
    "Annotate Once with\nYour Fixed Recipe", purple, 2.3
  ) +
  box_layer(
    3.30, 5.45, 4.15, 5.05, light_grey,
    "Build Cohort\nBlueprint", grey
  ) +
  arrow_layer(5.50, 6.15, 7.60, 6.35, blue) +
  arrow_layer(5.50, 6.15, 6.10, 6.10, purple) +
  annotate(
    "curve", x = 4.38, xend = 4.38, y = 5.10, yend = 5.62,
    curvature = -0.25, colour = grey, linewidth = 0.62,
    arrow = grid::arrow(length = grid::unit(0.105, "inches"), type = "closed")
  ) +
  box_layer(
    6.20, 8.75, 5.45, 6.75, light_green,
    "Recipe-Specific\nCache Ready", green, 2.75
  ) +
  arrow_layer(7.48, 7.48, 5.38, 3.45, green, "dashed", 0.72) +
  box_layer(0.55, 2.30, 2.15, 3.15, light_orange, "New Sample", orange, 3.0) +
  arrow_layer(2.37, 2.90, 2.65, 2.65, orange) +
  box_layer(
    2.95, 5.20, 2.15, 3.15, light_grey,
    "Normalize and\nSplit Records", grey, 2.8
  ) +
  arrow_layer(5.27, 5.85, 2.65, 2.65, grey) +
  box_layer(5.90, 8.55, 2.15, 3.15, light_blue, "Cache Lookup", blue, 3.0) +
  arrow_layer(8.62, 9.15, 2.78, 3.30, green) +
  arrow_layer(8.62, 9.15, 2.52, 1.38, purple) +
  box_layer(
    9.20, 11.70, 2.85, 3.75, light_green,
    "Reuse Cached\nAnnotations", green, 2.85
  ) +
  box_layer(
    9.20, 11.70, 0.75, 1.65, light_purple,
    "Annotate Cache Misses\nwith the Same Recipe", purple, 2.15
  ) +
  arrow_layer(11.77, 12.25, 3.30, 2.70, green) +
  arrow_layer(11.77, 12.25, 1.20, 2.35, purple) +
  box_layer(
    12.30, 14.15, 2.05, 3.05, light_grey,
    "Merge, Sort\nand Index", grey, 2.8
  ) +
  arrow_layer(14.22, 14.72, 2.55, 2.55, grey) +
  box_layer(
    14.77, 16.75, 2.05, 3.05, light_green,
    "Complete\nAnnotated BCF", green, 2.35
  ) +
  theme_void(base_family = "Lato") +
  theme(plot.margin = margin(10, 8, 10, 8))

figure1_v3 <- figure1_v3_base +
  coord_cartesian(xlim = c(0, 16.9), ylim = c(0.45, 8.55), clip = "off")

save_plot(
  figure1_v3,
  file.path(draft_dir, "figures", "main", "figure1_workflow_v3"),
  width = 7.2, height = 4.9, dpi = 600
)
save_plot(
  figure1_v3,
  file.path(draft_dir, "figures", "main", "panels", "figure1_v3_complete_workflow"),
  width = 7.2, height = 4.9, dpi = 600
)
figure1a_v3 <- figure1_v3_base +
  coord_cartesian(xlim = c(0, 9.1), ylim = c(3.95, 8.55), clip = "off")
figure1b_v3 <- figure1_v3_base +
  coord_cartesian(xlim = c(0, 16.9), ylim = c(0.45, 3.75), clip = "off")
save_plot(
  figure1a_v3,
  file.path(draft_dir, "figures", "main", "panels", "figure1_v3_A_cache_paths"),
  width = 7.2, height = 3.7, dpi = 600
)
save_plot(
  figure1b_v3,
  file.path(draft_dir, "figures", "main", "panels", "figure1_v3_B_new_sample"),
  width = 7.2, height = 2.9, dpi = 600
)

# Figure 3 v3: make increasing absolute savings and bounded relative speedup explicit.
complexity$cache_strategy <- factor(
  complexity$cache_strategy,
  levels = c("AF ≥ 10% cache", "AF ≥ 1% cache")
)
levels(complexity$cache_strategy) <- c("gnomAD AF ≥ 10%", "gnomAD AF ≥ 1%")
complexity$pipeline_label <- factor(
  paste0(complexity$delay_ms, " ms"),
  levels = paste0(sort(unique(complexity$delay_ms), decreasing = TRUE), " ms")
)
complexity$hours_saved <- complexity$time_saved_seconds / 3600
complexity$hours_saved_label <- sprintf("%.1f h", complexity$hours_saved)

cache_colours <- c("gnomAD AF ≥ 10%" = orange, "gnomAD AF ≥ 1%" = green)

p3a_v3 <- ggplot(complexity, aes(y = pipeline_label, colour = cache_strategy)) +
  geom_segment(
    aes(x = direct_hours, xend = cached_hours, yend = pipeline_label),
    linewidth = 0.9,
    arrow = grid::arrow(length = grid::unit(0.075, "inches"), type = "closed")
  ) +
  geom_point(aes(x = direct_hours), shape = 21, fill = "white", size = 2.4, stroke = 0.7) +
  geom_point(aes(x = cached_hours), shape = 17, size = 2.7) +
  geom_label(
    aes(x = (direct_hours + cached_hours) / 2, label = hours_saved_label),
    fill = "white", linewidth = 0, size = 2.45, show.legend = FALSE
  ) +
  facet_wrap(~cache_strategy, nrow = 1) +
  scale_colour_manual(values = cache_colours, guide = "none") +
  scale_x_continuous(breaks = c(0, 5, 10, 15, 20), limits = c(0, 24)) +
  labs(
    tag = "A", x = "Observed runtime per genome (hours)",
    y = "Added work per consequence"
  ) +
  journal_theme(9.2)

ceilings <- do.call(
  rbind,
  lapply(levels(complexity$cache_strategy), function(strategy) {
    values <- complexity[complexity$cache_strategy == strategy, ]
    hit <- median(values$cache_hit_rate)
    data.frame(
      cache_strategy = factor(strategy, levels = levels(complexity$cache_strategy)),
      ceiling = 1 / (1 - hit),
      label = sprintf("%.1f%% hits: %.1f× ceiling", 100 * hit, 1 / (1 - hit))
    )
  })
)

p3b_v3 <- ggplot(
  complexity,
  aes(x = direct_hours, y = speedup, colour = cache_strategy, group = cache_strategy)
) +
  geom_hline(
    data = ceilings, aes(yintercept = ceiling, colour = cache_strategy),
    linetype = "dashed", linewidth = 0.65, show.legend = FALSE
  ) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.3) +
  geom_text(
    data = ceilings,
    aes(x = 23.2, y = ceiling, label = label, colour = cache_strategy),
    inherit.aes = FALSE, hjust = 1, vjust = -0.55, size = 2.55,
    show.legend = FALSE
  ) +
  scale_colour_manual(values = cache_colours) +
  scale_x_continuous(breaks = c(5, 10, 15, 20), limits = c(3.5, 24)) +
  scale_y_continuous(breaks = c(2, 4, 6, 8, 10), limits = c(1, 11.1)) +
  labs(
    tag = "B", x = "Direct pipeline runtime (hours)",
    y = "Speedup over direct annotation"
  ) +
  journal_theme(9.2) +
  theme(legend.position = "bottom")

figure3_v3 <- p3a_v3 / p3b_v3 +
  plot_layout(heights = c(1.08, 1), guides = "collect") &
  theme(legend.position = "bottom")

save_plot(
  figure3_v3,
  file.path(draft_dir, "figures", "main", "figure3_pipeline_complexity_v3"),
  width = 7.2, height = 7.2, dpi = 600
)
save_plot(
  p3a_v3,
  file.path(draft_dir, "figures", "main", "panels", "figure3_v3_A_runtime_returned"),
  width = 7.2, height = 3.8, dpi = 600
)
save_plot(
  p3b_v3,
  file.path(draft_dir, "figures", "main", "panels", "figure3_v3_B_speedup_ceiling"),
  width = 5.0, height = 3.8, dpi = 600
)

# Figure 2 v3: rendered once the independent matched assay extension is complete.
assay_path <- file.path(
  repo_dir, "benchmarks", "external_assay_v3", "source_data",
  "external_assay_v3_metrics.tsv"
)
if (file.exists(assay_path)) {
  assay <- read_tsv(assay_path)
  if (nrow(assay) != 72 || !all(as_logical_semantic(assay$semantic_pass))) {
    stop("Expected 72 semantically validated matched assay observations")
  }
  assay$assay_label <- factor(
    toupper(assay$assay), levels = c("PANEL", "WES", "WGS"),
    labels = c("Panel", "WES", "WGS")
  )
  assay$tool_label <- factor(
    ifelse(assay$tool == "vep", "VEP", "fastVEP"),
    levels = c("VEP", "fastVEP")
  )
  assay$time_saved_percent <- 100 * (1 - assay$relative_runtime)

  summary_keys <- unique(assay[c("tool_label", "assay_label")])
  assay_summary <- do.call(
    rbind,
    lapply(seq_len(nrow(summary_keys)), function(index) {
      key <- summary_keys[index, ]
      values <- assay[
        assay$tool_label == key$tool_label & assay$assay_label == key$assay_label,
      ]
      data.frame(
        tool_label = key$tool_label,
        assay_label = key$assay_label,
        direct_seconds = median(values$uncached_wall_seconds),
        cached_seconds = median(values$cached_wall_seconds),
        speedup = median(values$speedup),
        time_saved_percent = median(values$time_saved_percent)
      )
    })
  )
  assay_summary$direction <- factor(
    ifelse(
      assay_summary$cached_seconds < assay_summary$direct_seconds,
      "Time returned", "Additional overhead"
    ),
    levels = c("Time returned", "Additional overhead")
  )
  assay_summary$tile_label <- ifelse(
    assay_summary$time_saved_percent >= 0,
    sprintf(
      "%.2f×\n%.0f%% saved",
      assay_summary$speedup, assay_summary$time_saved_percent
    ),
    sprintf(
      "%.2f×\n%.0f%% longer",
      assay_summary$speedup, abs(assay_summary$time_saved_percent)
    )
  )

  direction_colours <- c(
    "Time returned" = green,
    "Additional overhead" = "#A5564B"
  )
  runtime_breaks <- c(5, 15, 60, 300, 1800, 7200)
  runtime_labels <- c("5 s", "15 s", "1 min", "5 min", "30 min", "2 h")

  p2a_v3 <- ggplot(
    assay_summary,
    aes(y = assay_label, colour = direction)
  ) +
    geom_segment(
      aes(
        x = direct_seconds, xend = cached_seconds,
        yend = assay_label
      ),
      linewidth = 1.0,
      arrow = grid::arrow(length = grid::unit(0.075, "inches"), type = "closed")
    ) +
    geom_point(
      aes(x = direct_seconds), shape = 21, fill = "white",
      size = 2.5, stroke = 0.7, show.legend = FALSE
    ) +
    geom_point(
      aes(x = cached_seconds), shape = 17, size = 2.8,
      show.legend = FALSE
    ) +
    facet_wrap(~tool_label, nrow = 1) +
    scale_colour_manual(values = direction_colours) +
    scale_x_log10(
      breaks = runtime_breaks, labels = runtime_labels,
      limits = c(4, 15000)
    ) +
    labs(
      tag = "A", x = "Median runtime per sample (log scale)",
      y = NULL, colour = NULL
    ) +
    journal_theme(9.2) +
    theme(legend.position = "bottom")

  p2b_v3 <- ggplot(
    assay_summary,
    aes(x = assay_label, y = tool_label, fill = time_saved_percent)
  ) +
    geom_tile(colour = "white", linewidth = 2.0) +
    geom_text(
      aes(label = tile_label), family = "Lato", color = ink,
      fontface = "bold", size = 3.0, lineheight = 0.95
    ) +
    scale_fill_gradient2(
      low = "#C98578", mid = "#F5F1EC", high = "#76A984",
      midpoint = 0, limits = c(-200, 100), oob = scales::squish,
      labels = function(x) paste0(x, "%"), name = "Median wall time saved"
    ) +
    scale_y_discrete(limits = c("fastVEP", "VEP")) +
    labs(tag = "B", x = NULL, y = NULL) +
    theme_minimal(base_size = 9.2, base_family = "Lato") +
    theme(
      text = element_text(color = ink),
      axis.text = element_text(color = ink),
      panel.grid = element_blank(),
      legend.position = "bottom",
      legend.title = element_text(face = "plain"),
      plot.tag = element_text(face = "bold", size = 12, color = ink),
      plot.tag.position = c(0.005, 0.995),
      plot.margin = margin(8, 8, 8, 8)
    )

  figure2_v3 <- p2a_v3 + p2b_v3 +
    plot_layout(widths = c(1.25, 1), guides = "collect") &
    theme(legend.position = "bottom")
  save_plot(
    figure2_v3,
    file.path(draft_dir, "figures", "main", "figure2_assay_annotator_v3"),
    width = 7.2, height = 4.5, dpi = 600
  )
  save_plot(
    p2a_v3,
    file.path(draft_dir, "figures", "main", "panels", "figure2_v3_A_runtime"),
    width = 4.8, height = 4.0, dpi = 600
  )
  save_plot(
    p2b_v3,
    file.path(draft_dir, "figures", "main", "panels", "figure2_v3_B_use_case"),
    width = 4.2, height = 4.0, dpi = 600
  )
  write_tsv(
    assay,
    file.path(draft_dir, "figures", "source_data", "external_assay_v3_metrics.tsv")
  )

  cohort_colours <- c(KPGP = blue, SGDP = orange)
  cohort_shapes <- c(KPGP = 16, SGDP = 17)
  assay$cohort_label <- factor(toupper(assay$cohort), levels = c("KPGP", "SGDP"))
  coverage <- assay[assay$tool == "vep", ]

  s3a <- ggplot(
    coverage,
    aes(x = assay_label, y = 100 * cache_hit_rate)
  ) +
    geom_jitter(
      aes(colour = cohort_label, shape = cohort_label),
      width = 0.10, height = 0, size = 1.6, alpha = 0.78, stroke = 0.25
    ) +
    stat_summary(fun = median, geom = "point", shape = 95, size = 8, colour = ink) +
    scale_colour_manual(values = cohort_colours, name = "Cohort") +
    scale_shape_manual(values = cohort_shapes, name = "Cohort") +
    scale_y_continuous(labels = function(x) paste0(x, "%")) +
    labs(tag = "A", x = NULL, y = "Variants found in AF ≥ 1% cache") +
    journal_theme(9.2) +
    guides(colour = "none", shape = "none") +
    theme(legend.position = "none")

  s3b <- ggplot(assay, aes(x = assay_label, y = speedup)) +
    geom_hline(yintercept = 1, colour = grey, linetype = "dashed", linewidth = 0.55) +
    geom_jitter(
      aes(colour = cohort_label, shape = cohort_label),
      width = 0.10, height = 0, size = 1.45, alpha = 0.75, stroke = 0.25
    ) +
    stat_summary(fun = median, geom = "point", shape = 95, size = 7.5, colour = ink) +
    facet_wrap(~tool_label, nrow = 1) +
    scale_colour_manual(values = cohort_colours, name = "Cohort") +
    scale_shape_manual(values = cohort_shapes, name = "Cohort") +
    scale_y_log10(
      breaks = c(0.25, 0.5, 1, 2, 5, 10),
      labels = c("0.25×", "0.5×", "1×", "2×", "5×", "10×")
    ) +
    labs(tag = "B", x = NULL, y = "Speedup over direct annotation") +
    journal_theme(9.2)

  s3b <- s3b +
    guides(
      colour = "none",
      shape = guide_legend(
        override.aes = list(colour = unname(cohort_colours))
      )
    )

  supplementary3_v3 <- s3a / s3b +
    plot_layout(heights = c(0.9, 1.1), guides = "collect") &
    theme(legend.position = "bottom")
  save_plot(
    supplementary3_v3,
    file.path(
      draft_dir, "supplement", "figures",
      "supplementary_figure3_matched_assays_v3"
    ),
    width = 7.2, height = 6.8, dpi = 600
  )
  save_plot(
    s3a,
    file.path(
      draft_dir, "supplement", "figures", "panels",
      "supplementary_figure3_v3_A_cache_coverage"
    ),
    width = 4.6, height = 3.5, dpi = 600
  )
  save_plot(
    s3b,
    file.path(
      draft_dir, "supplement", "figures", "panels",
      "supplementary_figure3_v3_B_speedup"
    ),
    width = 6.2, height = 3.8, dpi = 600
  )
} else {
  message("Figure 2 v3 deferred until the independent matched assay table exists")
}

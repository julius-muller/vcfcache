render_assay_figure <- function(input_dir, output_dir, snapshot) {
  primary <- read_tsv(file.path(input_dir, "primary_wgs_metrics.tsv"))
  assay <- read_tsv(file.path(input_dir, "assay_metrics.tsv"))
  rows <- rbind(primary, assay)

  if (!all(as_logical_semantic(rows$semantic_pass))) {
    stop("Assay figure input contains a failed semantic comparison")
  }

  assay_levels <- c("panel", "wes", "wgs")
  assay_labels <- c(
    panel = "1000G-derived panel",
    wes = "1000G-derived WES",
    wgs = "1000G WGS"
  )
  assay_colors <- c(
    panel = vcf_colors[["green"]],
    wes = vcf_colors[["orange"]],
    wgs = vcf_colors[["blue"]]
  )
  rows <- rows[rows$assay %in% assay_levels, ]
  rows$assay <- factor(rows$assay, levels = assay_levels)
  rows$assay_label <- factor(
    assay_labels[as.character(rows$assay)],
    levels = unname(assay_labels[assay_levels])
  )
  rows$miss_rate <- 1 - rows$cache_hit_rate

  summaries <- do.call(
    rbind,
    lapply(seq_along(assay_levels), function(index) {
      name <- assay_levels[[index]]
      group <- rows[rows$assay == name, ]
      interval <- bootstrap_median_interval(
        group$relative_runtime,
        seed = 20260803 + index
      )
      data.frame(
        assay = name,
        assay_label = assay_labels[[name]],
        samples = length(unique(group$sample)),
        median_input_records = median(group$input_records),
        median_hit_rate = median(group$cache_hit_rate),
        hit_rate_q1 = unname(quantile(group$cache_hit_rate, 0.25)),
        hit_rate_q3 = unname(quantile(group$cache_hit_rate, 0.75)),
        median_relative_runtime = median(group$relative_runtime),
        relative_runtime_q1 = unname(quantile(group$relative_runtime, 0.25)),
        relative_runtime_q3 = unname(quantile(group$relative_runtime, 0.75)),
        relative_runtime_bootstrap_low = interval[[1]],
        relative_runtime_bootstrap_high = interval[[2]],
        median_speedup = median(group$speedup),
        median_uncached_seconds = median(group$uncached_wall_seconds),
        median_cached_seconds = median(group$cached_wall_seconds),
        total_wall_hours_saved = sum(group$wall_seconds_saved) / 3600,
        semantic_passes = sum(as_logical_semantic(group$semantic_pass)),
        unexpected_annotation_mismatches = sum(group$annotation_mismatches),
        ignored_known_vep_mismatches = sum(group$ignored_annotation_mismatches),
        annotation_order_only = sum(group$annotation_order_only),
        stringsAsFactors = FALSE
      )
    })
  )
  write_tsv(summaries, file.path(output_dir, "source", "assay_summary.tsv"))

  p_runtime <- ggplot(rows, aes(x = assay_label, y = relative_runtime, color = assay)) +
    geom_hline(yintercept = 1, color = vcf_colors[["red"]], linewidth = 0.6, linetype = 2) +
    geom_boxplot(width = 0.55, outlier.shape = NA, color = vcf_colors[["ink"]], fill = "white") +
    geom_point(
      position = position_jitter(width = 0.16, height = 0, seed = 20260805),
      alpha = 0.58, size = 1.8
    ) +
    scale_color_manual(values = assay_colors, guide = "none") +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1), expand = expansion(mult = c(0.03, 0.08))) +
    labs(
      title = "Observed time remaining",
      subtitle = "Source-overlap calibration only; line = no runtime benefit",
      x = NULL,
      y = "Cached / uncached wall time"
    ) +
    theme_vcfcache()

  absolute <- rbind(
    data.frame(
      assay = summaries$assay,
      assay_label = summaries$assay_label,
      mode = "Without VCFcache",
      seconds = summaries$median_uncached_seconds
    ),
    data.frame(
      assay = summaries$assay,
      assay_label = summaries$assay_label,
      mode = "With VCFcache",
      seconds = summaries$median_cached_seconds
    )
  )
  absolute$assay_label <- factor(absolute$assay_label, levels = unname(assay_labels))
  segments <- data.frame(
    assay_label = factor(summaries$assay_label, levels = unname(assay_labels)),
    before = summaries$median_uncached_seconds / 60,
    after = summaries$median_cached_seconds / 60
  )
  p_absolute <- ggplot(absolute, aes(x = assay_label, y = seconds / 60)) +
    geom_segment(
      data = segments,
      aes(x = assay_label, xend = assay_label, y = before, yend = after),
      inherit.aes = FALSE, color = vcf_colors[["light_grey"]], linewidth = 3
    ) +
    geom_point(aes(fill = mode), shape = 21, size = 4, color = "white", stroke = 0.8) +
    scale_fill_manual(values = c("Without VCFcache" = vcf_colors[["grey"]], "With VCFcache" = vcf_colors[["green"]])) +
    scale_y_log10(labels = scales::label_number(accuracy = 0.1), breaks = c(0.25, 0.5, 1, 3, 10, 30, 100, 300)) +
    labs(
      title = "Median wait per sample",
      subtitle = "The small panel workload exposes fixed lookup overhead",
      x = NULL,
      y = "Wall time (minutes, log scale)",
      fill = NULL
    ) +
    theme_vcfcache() +
    theme(legend.position = "bottom")

  ideal <- data.frame(miss_rate = seq(0, 0.30, length.out = 200))
  ideal$relative_runtime <- ideal$miss_rate
  p_mechanism <- ggplot(rows, aes(x = miss_rate, y = relative_runtime, color = assay)) +
    geom_line(
      data = ideal,
      aes(x = miss_rate, y = relative_runtime),
      inherit.aes = FALSE,
      color = vcf_colors[["grey"]], linewidth = 0.8, linetype = 2
    ) +
    geom_point(alpha = 0.62, size = 2) +
    scale_color_manual(values = assay_colors, labels = assay_labels, name = NULL) +
    scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    labs(
      title = "Misses drive annotation work",
      subtitle = "Runtime mechanics conditional on observed hit rate",
      x = "Variants not found in cache",
      y = "Cached / uncached wall time"
    ) +
    theme_vcfcache()

  total_passes <- sum(summaries$semantic_passes)
  unexpected <- sum(summaries$unexpected_annotation_mismatches)
  ignored <- sum(summaries$ignored_known_vep_mismatches)
  trust <- ggplot() +
    annotate("rect", xmin = 0, xmax = 1, ymin = 0, ymax = 1, fill = "#EAF6EE", color = NA) +
    annotate(
      "text", x = 0.03, y = 0.68, hjust = 0, size = 6.2, fontface = "bold",
      color = vcf_colors[["green"]], label = sprintf("✓  %d/%d paired outputs pass", total_passes, total_passes)
    ) +
    annotate(
      "text", x = 0.03, y = 0.30, hjust = 0, size = 3.6,
      color = vcf_colors[["ink"]],
      label = sprintf(
        "Unexpected annotation mismatches: %d     |     Known upstream VEP HGNC_ID differences reported separately: %d",
        unexpected, ignored
      )
    ) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
    theme_void()

  combined <- (p_runtime | p_absolute | p_mechanism) / trust +
    plot_layout(heights = c(1, 0.23)) +
    plot_annotation(
      title = "Source-overlap calibration across assay sizes",
      subtitle = "1000 Genomes-derived inputs · bundled Zenodo gnomAD AF ≥ 1% cache · VEP 115.2 --everything",
      caption = paste(
        preliminary_caption(snapshot),
        "CALIBRATION ONLY: these identities occur in the gnomAD source universe; their hit rates and speedups are not real-world estimates.",
        "Use the external KPGP, SGDP and PGP cohorts for the primary WGS benefit estimate.",
        sep = "\n"
      ),
      tag_levels = "A",
      theme = theme(
        plot.title = element_text(size = 18, face = "bold", color = vcf_colors[["ink"]]),
        plot.subtitle = element_text(size = 11, color = vcf_colors[["grey"]]),
        plot.caption = element_text(size = 8, color = vcf_colors[["grey"]], hjust = 0),
        plot.tag = element_text(face = "bold", color = vcf_colors[["blue"]])
      )
    )

  panel_dir <- file.path(
    output_dir, "manuscript", "panels", "assay_source_overlap_calibration"
  )
  save_plot(p_runtime, file.path(panel_dir, "A_time_remaining"), 5.2, 4.3)
  save_plot(p_absolute, file.path(panel_dir, "B_median_wait"), 5.2, 4.3)
  save_plot(p_mechanism, file.path(panel_dir, "C_miss_fraction"), 5.2, 4.3)
  save_plot(trust, file.path(panel_dir, "D_output_correctness"), 10.5, 2.1)
  figure_name <- "assay_source_overlap_calibration"
  save_plot(combined, file.path(output_dir, "manuscript", figure_name), 15, 6.5)
  invisible(list(plot = combined, summary = summaries))
}

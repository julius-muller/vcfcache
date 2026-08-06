render_user_impact_figure <- function(input_dir, output_dir, snapshot) {
  scenarios <- data.frame(
    situation = c("Distant WGS", "Related/cohort WGS", "WES", "Very high reuse"),
    hit_rate = c(0.50, 0.80, 0.90, 0.95),
    stringsAsFactors = FALSE
  )
  scenarios$situation <- factor(scenarios$situation, levels = rev(scenarios$situation))
  scenarios$relative_runtime <- 1 - scenarios$hit_rate
  scenarios$speedup <- 1 / scenarios$relative_runtime

  scenario_bars <- rbind(
    data.frame(
      panel = "time_left",
      label = scenarios$situation,
      segment = "Time remaining",
      fraction = scenarios$relative_runtime,
      hit_rate = scenarios$hit_rate,
      sample_count = 1,
      baseline_seconds = 100,
      cached_seconds = scenarios$relative_runtime * 100,
      saved_seconds = scenarios$hit_rate * 100
    ),
    data.frame(
      panel = "time_left",
      label = scenarios$situation,
      segment = "Time returned",
      fraction = scenarios$hit_rate,
      hit_rate = scenarios$hit_rate,
      sample_count = 1,
      baseline_seconds = 100,
      cached_seconds = scenarios$relative_runtime * 100,
      saved_seconds = scenarios$hit_rate * 100
    )
  )
  scenario_bars$segment <- factor(
    scenario_bars$segment,
    levels = c("Time remaining", "Time returned")
  )

  pipeline <- data.frame(
    label = c("10 min", "1 hour", "10 hours"),
    baseline_seconds = c(10, 60, 600) * 60,
    hit_rate = 0.80,
    stringsAsFactors = FALSE
  )
  pipeline$cached_seconds <- pipeline$baseline_seconds * (1 - pipeline$hit_rate)
  pipeline$saved_seconds <- pipeline$baseline_seconds - pipeline$cached_seconds
  pipeline$panel <- "pipeline_cost"
  pipeline$sample_count <- 1

  sample_scale <- data.frame(
    sample_count = c(10, 100, 1000),
    baseline_seconds = c(10, 100, 1000) * 3600,
    hit_rate = 0.80
  )
  sample_scale$cached_seconds <- sample_scale$baseline_seconds * (1 - sample_scale$hit_rate)
  sample_scale$saved_seconds <- sample_scale$baseline_seconds - sample_scale$cached_seconds
  sample_scale$label <- paste(format(sample_scale$sample_count, big.mark = ",", scientific = FALSE), "samples")
  sample_scale$panel <- "sample_scale"

  source_rows <- rbind(
    scenario_bars,
    data.frame(
      panel = pipeline$panel,
      label = pipeline$label,
      segment = "",
      fraction = pipeline$cached_seconds / pipeline$baseline_seconds,
      hit_rate = pipeline$hit_rate,
      sample_count = pipeline$sample_count,
      baseline_seconds = pipeline$baseline_seconds,
      cached_seconds = pipeline$cached_seconds,
      saved_seconds = pipeline$saved_seconds
    ),
    data.frame(
      panel = sample_scale$panel,
      label = sample_scale$label,
      segment = "",
      fraction = sample_scale$cached_seconds / sample_scale$baseline_seconds,
      hit_rate = sample_scale$hit_rate,
      sample_count = sample_scale$sample_count,
      baseline_seconds = sample_scale$baseline_seconds,
      cached_seconds = sample_scale$cached_seconds,
      saved_seconds = sample_scale$saved_seconds
    )
  )
  write_tsv(source_rows, file.path(output_dir, "source", "user_impact_model_preview.tsv"))

  p_scenarios <- ggplot(scenario_bars, aes(x = fraction, y = label, fill = segment)) +
    geom_col(
      width = 0.62,
      color = "white",
      linewidth = 0.5,
      position = position_stack(reverse = TRUE)
    ) +
    geom_text(
      data = scenarios,
      aes(
        x = 1.03,
        y = situation,
        label = sprintf("%.0f%% hit  →  %.1f× faster", hit_rate * 100, speedup)
      ),
      inherit.aes = FALSE,
      hjust = 0, size = 3.8, fontface = "bold", color = vcf_colors[["ink"]]
    ) +
    scale_fill_manual(
      values = c("Time remaining" = vcf_colors[["blue"]], "Time returned" = vcf_colors[["light_green"]]),
      breaks = c("Time remaining", "Time returned"), name = NULL
    ) +
    scale_x_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1.42), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
    labs(
      title = "Find your expected cache reuse",
      subtitle = "Dark = annotation work left · green = waiting time returned",
      x = NULL,
      y = NULL
    ) +
    theme_vcfcache(base_size = 12) +
    theme(legend.position = "top", panel.grid = element_blank())

  pipeline$today <- paste("Today:", pipeline$label)
  pipeline$with_cache <- paste("With cache:", format_duration(pipeline$cached_seconds))
  pipeline$returned <- paste("Get back", format_duration(pipeline$saved_seconds))
  pipeline$x <- seq_len(nrow(pipeline))
  p_pipeline <- ggplot(pipeline, aes(x = x, y = 1)) +
    geom_tile(width = 0.88, height = 0.90, fill = "#EDF5F8", color = "white", linewidth = 1.5) +
    geom_text(aes(y = 1.22, label = today), size = 3.8, color = vcf_colors[["grey"]]) +
    geom_text(aes(y = 1.00, label = with_cache), size = 4.5, fontface = "bold", color = vcf_colors[["blue"]]) +
    geom_text(aes(y = 0.78, label = returned), size = 3.8, fontface = "bold", color = vcf_colors[["green"]]) +
    coord_cartesian(xlim = c(0.5, 3.5), ylim = c(0.5, 1.5), clip = "off") +
    labs(
      title = "Fast or slow pipeline?",
      subtitle = "At 80% reuse, the same fraction of expensive annotation work disappears"
    ) +
    theme_void(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 14, color = vcf_colors[["ink"]]),
      plot.subtitle = element_text(size = 10, color = vcf_colors[["grey"]], margin = margin(b = 8)),
      plot.margin = margin(12, 14, 12, 14)
    )

  sample_scale$x <- seq_len(nrow(sample_scale))
  sample_scale$returned <- paste("Get back", format_duration(sample_scale$saved_seconds))
  p_scale <- ggplot(sample_scale, aes(x = x, y = 1)) +
    geom_tile(width = 0.90, height = 0.88, fill = "white", color = vcf_colors[["light_grey"]], linewidth = 1.3) +
    geom_text(aes(y = 1.18, label = label), size = 4.3, fontface = "bold", color = vcf_colors[["ink"]]) +
    geom_text(aes(y = 0.88, label = returned), size = 3.8, fontface = "bold", color = vcf_colors[["green"]]) +
    coord_cartesian(xlim = c(0.45, 3.55), ylim = c(0.55, 1.45), clip = "off") +
    labs(
      title = "Ten samples or one thousand: savings add up",
      subtitle = "Example: 1-hour pipeline, 80% reuse, bundled cache with no build cost"
    ) +
    theme_void(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 14, color = vcf_colors[["ink"]]),
      plot.subtitle = element_text(size = 10, color = vcf_colors[["grey"]], margin = margin(b = 8)),
      plot.margin = margin(12, 14, 12, 14)
    )

  validated <- snapshot$sample_counts$primary_wgs +
    sum(unlist(snapshot$sample_counts$assay_by_type)) +
    snapshot$sample_counts$external_semantically_validated
  p_trust <- ggplot() +
    annotate("rect", xmin = 0, xmax = 1, ymin = 0, ymax = 1, fill = "#EAF6EE", color = NA) +
    annotate(
      "text", x = 0.06, y = 0.72, hjust = 0, size = 7, fontface = "bold",
      color = vcf_colors[["green"]], label = "✓ Same annotations"
    ) +
    annotate(
      "text", x = 0.06, y = 0.43, hjust = 0, size = 4.2,
      color = vcf_colors[["ink"]], label = sprintf("%d/%d benchmark samples pass output-equivalence checks", validated, validated)
    ) +
    annotate(
      "text", x = 0.06, y = 0.20, hjust = 0, size = 3.4,
      color = vcf_colors[["grey"]], label = "Any unexpected annotation change fails the task."
    ) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
    theme_void()

  combined <- (p_scenarios | p_pipeline) / (p_scale | p_trust) +
    plot_layout(widths = c(1.12, 0.88), heights = c(1.08, 0.72)) +
    plot_annotation(
      title = "How much annotation time do you get back?",
      subtitle = "New time ≈ lookup + (1 − cache hit rate) × your current annotation time",
      caption = paste(
        "MODELED PREVIEW: zero lookup overhead is shown until the controlled-runtime campaign supplies the measured value.",
        "Very short pipelines can be dominated by lookup overhead; source-overlap calibration runs demonstrate this break-even behavior.",
        preliminary_caption(snapshot),
        sep = "\n"
      ),
      theme = theme(
        plot.background = element_rect(fill = vcf_colors[["paper"]], color = NA),
        plot.title = element_text(size = 24, face = "bold", color = vcf_colors[["ink"]]),
        plot.subtitle = element_text(size = 13, color = vcf_colors[["grey"]]),
        plot.caption = element_text(size = 8, color = vcf_colors[["red"]], hjust = 0)
      )
    )

  panel_dir <- file.path(output_dir, "repository", "panels", "user_impact")
  save_plot(p_scenarios, file.path(panel_dir, "A_cache_reuse"), 7.4, 5.0)
  save_plot(p_pipeline, file.path(panel_dir, "B_pipeline_duration"), 7.4, 5.0)
  save_plot(p_scale, file.path(panel_dir, "C_sample_scale"), 7.4, 3.7)
  save_plot(p_trust, file.path(panel_dir, "D_semantic_equivalence"), 7.4, 3.7)
  save_plot(combined, file.path(output_dir, "repository", "user_impact_model_preview"), 15, 9.2)
  invisible(list(plot = combined, source = source_rows))
}

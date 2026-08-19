#!/usr/bin/env Rscript

arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) != 2) {
  stop(
    "Usage: Rscript --vanilla benchmarks/figures/R/render_pipeline_complexity_dual.R ",
    "<dual-cache-snapshot-directory> <output-directory>"
  )
}

full_arguments <- commandArgs(trailingOnly = FALSE)
file_argument <- full_arguments[grepl("^--file=", full_arguments)]
if (length(file_argument) != 1) {
  stop("Could not determine the renderer location")
}
script_dir <- dirname(normalizePath(sub("^--file=", "", file_argument)))
source(file.path(script_dir, "common.R"))

input_dir <- normalizePath(arguments[[1]], mustWork = TRUE)
output_dir <- arguments[[2]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
snapshot <- jsonlite::fromJSON(file.path(input_dir, "SNAPSHOT.json"))
if (!identical(snapshot$status, "PIPELINE_COMPLEXITY_DUAL_CACHE_FINAL")) {
  stop("Unexpected dual-cache snapshot status: ", snapshot$status)
}

af10 <- read_tsv(file.path(input_dir, "wgs_pipeline_spectrum_af10.tsv"))
af1 <- read_tsv(file.path(input_dir, "wgs_pipeline_spectrum_af1.tsv"))
annotators <- read_tsv(file.path(input_dir, "annotator_external_wgs_metrics.tsv"))
validate_rows <- function(rows) {
  nrow(rows) == 6 &&
    all(as_logical_semantic(rows$semantic_pass)) &&
    all(rows$statistics_mode == "light") &&
    length(unique(rows$sample)) == 1 && unique(rows$sample) == "KPGP-00319"
}
if (!validate_rows(af10) || !validate_rows(af1)) {
  stop("One dual-cache spectrum is incomplete or invalid")
}
if (!all(af10$delay_us == af1$delay_us) ||
    !isTRUE(all.equal(af10$uncached_wall_seconds, af1$uncached_wall_seconds))) {
  stop("Dual-cache spectra do not share identical direct baselines")
}

af10$cache_strategy <- "AF ≥ 10% cache"
af1$cache_strategy <- "AF ≥ 1% cache"
measurement_columns <- c(
  "sample", "cohort", "assembly", "pipeline", "delay_us", "input_records",
  "cache_hit_rate", "uncached_wall_seconds", "cached_wall_seconds",
  "relative_runtime", "speedup", "time_saved_seconds", "semantic_pass",
  "semantic_comparator", "statistics_mode", "cache_strategy"
)
combined <- rbind(af10[, measurement_columns], af1[, measurement_columns])
combined$cache_strategy <- factor(
  combined$cache_strategy,
  levels = c("AF ≥ 10% cache", "AF ≥ 1% cache")
)
combined$delay_ms <- combined$delay_us / 1000
delay_levels <- paste0(sort(unique(combined$delay_ms)), " ms")
combined$pipeline_label <- factor(
  paste0(combined$delay_ms, " ms"), levels = delay_levels
)
combined$hours_saved <- combined$time_saved_seconds / 3600
combined$direct_hours <- combined$uncached_wall_seconds / 3600
combined$cached_hours <- combined$cached_wall_seconds / 3600
combined$speedup_label <- sprintf("%.1f×", combined$speedup)
combined$remaining_label <- scales::percent(combined$relative_runtime, accuracy = 0.1)
cache_colors <- c(
  "AF ≥ 10% cache" = vcf_colors[["blue"]],
  "AF ≥ 1% cache" = vcf_colors[["green"]]
)

p_remaining <- ggplot(
  combined,
  aes(x = pipeline_label, y = relative_runtime, fill = cache_strategy)
) +
  geom_col(position = position_dodge(width = 0.76), width = 0.66, color = "white") +
  geom_text(
    aes(label = remaining_label),
    position = position_dodge(width = 0.76), vjust = -0.42,
    size = 3.0, fontface = "bold"
  ) +
  scale_fill_manual(values = cache_colors) +
  scale_y_continuous(
    limits = c(0, 0.255), breaks = seq(0, 0.25, 0.05), labels = scales::percent,
    expand = expansion(mult = c(0, 0.01))
  ) +
  labs(
    title = "A  More cache coverage removes more waiting",
    subtitle = "Measured share of direct wall time remaining",
    x = "Synthetic delay per transcript consequence", y = "Wall time remaining",
    fill = NULL
  ) +
  theme_vcfcache() +
  theme(legend.position = "top")

p_saved <- ggplot(
  combined,
  aes(x = pipeline_label, y = hours_saved, fill = cache_strategy)
) +
  geom_col(position = position_dodge(width = 0.76), width = 0.66, color = "white") +
  geom_text(
    aes(label = sprintf("%.1f h", hours_saved)),
    position = position_dodge(width = 0.76), vjust = -0.42,
    size = 3.0, fontface = "bold"
  ) +
  scale_fill_manual(values = cache_colors) +
  scale_y_continuous(
    limits = c(0, max(combined$hours_saved) * 1.14),
    expand = expansion(mult = c(0, 0.01))
  ) +
  labs(
    title = "B  Absolute savings grow with pipeline cost",
    subtitle = "Up to 20.5 hours returned for one real WGS",
    x = "Synthetic delay per transcript consequence", y = "Hours returned per genome",
    fill = NULL
  ) +
  theme_vcfcache() +
  theme(legend.position = "top")

ceilings <- aggregate(
  cache_hit_rate ~ cache_strategy, combined,
  function(value) 1 / (1 - unique(value))
)
p_speedup <- ggplot(
  combined,
  aes(
    x = direct_hours, y = speedup,
    color = cache_strategy, shape = cache_strategy
  )
) +
  geom_hline(
    data = ceilings,
    aes(yintercept = cache_hit_rate, color = cache_strategy),
    linetype = 2, linewidth = 0.75, show.legend = FALSE
  ) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 3.8, stroke = 1.0) +
  annotate(
    "text", x = max(combined$direct_hours),
    y = ceilings$cache_hit_rate[ceilings$cache_strategy == "AF ≥ 1% cache"] + 0.12,
    label = "90.26% hits → ~10.3× ceiling", hjust = 1,
    color = cache_colors[["AF ≥ 1% cache"]], size = 3.2
  ) +
  annotate(
    "text", x = max(combined$direct_hours),
    y = ceilings$cache_hit_rate[ceilings$cache_strategy == "AF ≥ 10% cache"] - 0.18,
    label = "80.23% hits → ~5.1× ceiling", hjust = 1,
    color = cache_colors[["AF ≥ 10% cache"]], size = 3.2
  ) +
  scale_color_manual(values = cache_colors) +
  scale_shape_manual(values = c("AF ≥ 10% cache" = 16, "AF ≥ 1% cache" = 17)) +
  scale_y_continuous(limits = c(4, 11), breaks = 4:11) +
  labs(
    title = "C  Speedup approaches the hit-rate ceiling",
    subtitle = "The larger bundled cache nearly doubles the slow-pipeline ceiling",
    x = "Direct annotation wall time (hours)", y = "Measured speedup",
    color = NULL, shape = NULL
  ) +
  theme_vcfcache() +
  theme(legend.position = "top")

matched <- annotators[
  annotators$sample == "KPGP-00319" & annotators$tool == "vep" &
    annotators$strategy %in% c("gnomad_af_0.1", "gnomad_af_0.01"),
]
lookup <- data.frame(
  cache_strategy = factor(
    c("AF ≥ 10% cache", "AF ≥ 1% cache"),
    levels = levels(combined$cache_strategy)
  ),
  overhead = c(
    matched$cached_wall_seconds[matched$strategy == "gnomad_af_0.1"] -
      (1 - matched$cache_hit_rate[matched$strategy == "gnomad_af_0.1"]) *
        matched$uncached_wall_seconds[matched$strategy == "gnomad_af_0.1"],
    matched$cached_wall_seconds[matched$strategy == "gnomad_af_0.01"] -
      (1 - matched$cache_hit_rate[matched$strategy == "gnomad_af_0.01"]) *
        matched$uncached_wall_seconds[matched$strategy == "gnomad_af_0.01"]
  )
)
model <- merge(combined, lookup, by = "cache_strategy")
model$predicted_relative_runtime <- (
  pmax(0, model$overhead) +
    (1 - model$cache_hit_rate) * model$uncached_wall_seconds
) / model$uncached_wall_seconds
model$residual_percentage_points <- 100 * (
  model$relative_runtime - model$predicted_relative_runtime
)
p_model <- ggplot(
  model,
  aes(
    x = predicted_relative_runtime, y = relative_runtime,
    color = cache_strategy, shape = cache_strategy
  )
) +
  geom_abline(slope = 1, intercept = 0, color = vcf_colors[["grey"]], linetype = 2) +
  geom_point(size = 3.8, stroke = 1.0) +
  scale_color_manual(values = cache_colors) +
  scale_shape_manual(values = c("AF ≥ 10% cache" = 16, "AF ≥ 1% cache" = 17)) +
  scale_x_continuous(labels = scales::percent, limits = c(0.08, 0.24)) +
  scale_y_continuous(labels = scales::percent, limits = c(0.08, 0.24)) +
  coord_equal() +
  labs(
    title = "D  Model agreement",
    subtitle = "Prediction versus measurement",
    x = "Predicted wall time", y = "Measured wall time",
    color = NULL, shape = NULL
  ) +
  theme_vcfcache() +
  theme(legend.position = "top")

supplement <- (p_remaining | p_saved) / (p_speedup | p_model) +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(
    title = "Bundled cache coverage determines how much pipeline time returns",
    subtitle = paste(
      "One 4.80-million-variant real WGS · matched direct baselines ·",
      "bundled Zenodo caches · light statistics · no technical repeats"
    ),
    caption = paste(
      "AF ≥ 10%: 80.23% hits; AF ≥ 1%: 90.26% hits. All 12 cached outputs passed semantic validation.",
      "Synthetic delays emit no annotation fields and affect only variants submitted to VEP."
    ),
    theme = theme(
      plot.title = element_text(face = "bold", size = 17, color = vcf_colors[["ink"]]),
      plot.subtitle = element_text(size = 11, color = vcf_colors[["grey"]]),
      plot.caption = element_text(size = 8.5, color = vcf_colors[["grey"]], hjust = 0)
    )
  )

selected_delays <- c(0, 4, 10)
vep_rows <- matched
vep_examples <- data.frame(
  example = "Standard VEP --everything",
  cache_strategy = factor(
    ifelse(vep_rows$strategy == "gnomad_af_0.1", "AF ≥ 10% cache", "AF ≥ 1% cache"),
    levels = levels(combined$cache_strategy)
  ),
  direct_hours = vep_rows$uncached_wall_seconds / 3600,
  cached_hours = vep_rows$cached_wall_seconds / 3600,
  stringsAsFactors = FALSE
)
pipeline_examples <- combined[combined$delay_ms %in% selected_delays[-1], ]
pipeline_examples$example <- ifelse(
  pipeline_examples$delay_ms == 4,
  "Heavier pipeline (4 ms)",
  "All-day pipeline (10 ms)"
)
headline <- rbind(
  vep_examples,
  pipeline_examples[, c("example", "cache_strategy", "direct_hours", "cached_hours")]
)
headline$example <- factor(
  headline$example,
  levels = c("Standard VEP --everything", "Heavier pipeline (4 ms)", "All-day pipeline (10 ms)")
)
headline$speedup <- headline$direct_hours / headline$cached_hours
headline$label <- paste0(
  sprintf("%.1f× · ", headline$speedup),
  format_duration(headline$cached_hours * 3600)
)
p_headline <- ggplot(
  headline,
  aes(x = example, y = cached_hours, fill = cache_strategy)
) +
  geom_col(position = position_dodge(width = 0.76), width = 0.65, color = "white") +
  geom_text(
    aes(label = label),
    position = position_dodge(width = 0.76), hjust = -0.1,
    fontface = "bold", size = 4.0
  ) +
  geom_point(
    data = headline[!duplicated(headline$example), ],
    aes(x = example, y = direct_hours), inherit.aes = FALSE,
    shape = 23, fill = "white", color = vcf_colors[["ink"]], size = 3.6, stroke = 1.2
  ) +
  coord_flip(clip = "off") +
  scale_fill_manual(values = cache_colors) +
  scale_y_continuous(
    limits = c(0, max(headline$direct_hours) * 1.22),
    expand = expansion(mult = c(0, 0.01)),
    labels = function(value) format_duration(value * 3600)
  ) +
  labs(
    title = "More cache coverage gives you more of the day back",
    subtitle = "Measured on the same 4.80-million-variant WGS; diamonds mark direct annotation time",
    x = NULL, y = "End-to-end wall time", fill = NULL,
    caption = paste(
      "Labels report measured speedup and remaining wall time.",
      "AF ≥ 10% cache: 80.23% hits; AF ≥ 1% cache: 90.26% hits.",
      "Both are bundled VCFcache downloads from Zenodo.",
      sep = "\n"
    )
  ) +
  theme_vcfcache(base_size = 14) +
  theme(
    legend.position = "top",
    panel.grid.major.y = element_blank(),
    axis.text.y = element_text(face = "bold", size = 12),
    plot.title = element_text(size = 21),
    plot.margin = margin(16, 85, 16, 16)
  )

incremental <- merge(
  af10[, c("delay_us", "time_saved_seconds")],
  af1[, c("delay_us", "time_saved_seconds")],
  by = "delay_us", suffixes = c("_af10", "_af1")
)
incremental$delay_ms <- incremental$delay_us / 1000
incremental$additional_hours <- (
  incremental$time_saved_seconds_af1 - incremental$time_saved_seconds_af10
) / 3600
p_incremental <- ggplot(
  incremental,
  aes(x = factor(delay_ms), y = additional_hours)
) +
  geom_col(fill = vcf_colors[["green"]], width = 0.68, color = "white") +
  geom_text(
    aes(label = sprintf("+%.1f h", additional_hours)),
    vjust = -0.45, fontface = "bold", size = 3.8
  ) +
  scale_y_continuous(
    limits = c(0, max(incremental$additional_hours) * 1.18),
    expand = expansion(mult = c(0, 0.01))
  ) +
  labs(
    title = "Additional hours returned by the AF ≥ 1% cache",
    subtitle = "Increment beyond the smaller AF ≥ 10% bundled cache",
    x = "Synthetic delay per transcript consequence", y = "Additional hours per genome"
  ) +
  theme_vcfcache()

write_tsv(combined, file.path(output_dir, "source", "dual_cache_pipeline_spectrum.tsv"))
write_tsv(model, file.path(output_dir, "source", "dual_cache_runtime_model.tsv"))
write_tsv(headline, file.path(output_dir, "source", "dual_cache_headline_examples.tsv"))
write_tsv(incremental, file.path(output_dir, "source", "af1_incremental_hours.tsv"))

panel_dir <- file.path(output_dir, "manuscript", "panels", "supplement_dual_cache_pipeline")
save_plot(p_remaining, file.path(panel_dir, "A_relative_runtime"), 6.2, 4.8)
save_plot(p_saved, file.path(panel_dir, "B_absolute_time_saved"), 6.2, 4.8)
save_plot(p_speedup, file.path(panel_dir, "C_speedup_ceiling"), 6.2, 4.8)
save_plot(p_model, file.path(panel_dir, "D_runtime_model"), 6.2, 4.8)
save_plot(
  supplement,
  file.path(output_dir, "manuscript", "supplement_dual_cache_pipeline"),
  12.4, 9.2
)
save_plot(
  p_headline,
  file.path(output_dir, "repository", "headliner_cache_coverage"),
  12, 6.6
)
save_plot(
  p_incremental,
  file.path(output_dir, "alternatives", "af1_additional_hours"),
  8.6, 5.4
)

writeLines(capture.output(sessionInfo()), file.path(output_dir, "R_SESSION_INFO_DUAL_CACHE.txt"))
writeLines(
  c(
    paste("snapshot", input_dir),
    paste("rendered_at", format(Sys.time(), tz = "UTC", usetz = TRUE)),
    paste("renderer", normalizePath(file.path(script_dir, "render_pipeline_complexity_dual.R")))
  ),
  file.path(output_dir, "RENDERED_FROM_DUAL_CACHE.txt")
)
cat("Rendered dual-cache pipeline figures under", output_dir, "\n")

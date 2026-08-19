#!/usr/bin/env Rscript

arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) != 2) {
  stop(
    "Usage: Rscript --vanilla benchmarks/figures/R/render_pipeline_complexity.R ",
    "<pipeline-snapshot-directory> <output-directory>"
  )
}

full_arguments <- commandArgs(trailingOnly = FALSE)
file_argument <- full_arguments[grepl("^--file=", full_arguments)]
if (length(file_argument) != 1) {
  stop("Could not determine the render_pipeline_complexity.R location")
}
script_dir <- dirname(normalizePath(sub("^--file=", "", file_argument)))
source(file.path(script_dir, "common.R"))

input_dir <- normalizePath(arguments[[1]], mustWork = TRUE)
output_dir <- arguments[[2]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
snapshot <- jsonlite::fromJSON(file.path(input_dir, "SNAPSHOT.json"))
if (!identical(snapshot$status, "PIPELINE_COMPLEXITY_WGS_FINAL")) {
  stop("Unexpected pipeline snapshot status: ", snapshot$status)
}

wgs <- read_tsv(file.path(input_dir, "wgs_pipeline_spectrum.tsv"))
annotators <- read_tsv(file.path(input_dir, "annotator_external_wgs_metrics.tsv"))
if (nrow(wgs) != 6 ||
    !all(as_logical_semantic(wgs$semantic_pass)) ||
    !all(wgs$statistics_mode == "light") ||
    length(unique(wgs$sample)) != 1 || unique(wgs$sample) != "KPGP-00319") {
  stop("WGS spectrum table is incomplete, invalid, or not light mode")
}

wgs <- wgs[order(wgs$delay_us), ]
wgs$delay_ms <- wgs$delay_us / 1000
wgs$pipeline_label <- factor(
  paste0(wgs$delay_ms, " ms/consequence"),
  levels = paste0(wgs$delay_ms, " ms/consequence")
)
wgs$fraction_saved <- 1 - wgs$relative_runtime
wgs$speedup_label <- sprintf("%.1f×", wgs$speedup)
wgs$saved_label <- format_duration(wgs$time_saved_seconds)
wgs$delay_group <- factor(wgs$pipeline_label, levels = levels(wgs$pipeline_label))
delay_colors <- setNames(
  grDevices::colorRampPalette(c(vcf_colors[["blue"]], vcf_colors[["purple"]]))(nrow(wgs)),
  levels(wgs$pipeline_label)
)

p_relative <- ggplot(
  wgs,
  aes(x = pipeline_label, y = relative_runtime, fill = pipeline_label)
) +
  geom_col(width = 0.68, color = "white") +
  geom_text(
    aes(label = speedup_label), vjust = -0.45,
    fontface = "bold", size = 3.4
  ) +
  scale_fill_manual(values = delay_colors, guide = "none") +
  scale_y_continuous(
    limits = c(0, 0.35), breaks = seq(0, 0.3, 0.1), labels = scales::percent,
    expand = expansion(mult = c(0, 0.01))
  ) +
  labs(
    title = "A  About 80% reuse removes most of every pipeline",
    subtitle = "Numbers above bars are measured speedups",
    x = "Synthetic plugin delay", y = "Wall time remaining"
  ) +
  theme_vcfcache() +
  theme(axis.text.x = element_text(angle = 18, hjust = 1))

p_saved <- ggplot(
  wgs,
  aes(x = pipeline_label, y = time_saved_seconds / 3600, fill = pipeline_label)
) +
  geom_col(width = 0.64, color = vcf_colors[["ink"]], linewidth = 0.35) +
  geom_text(
    aes(label = saved_label), hjust = -0.12,
    fontface = "bold", size = 3.8
  ) +
  coord_flip(clip = "off") +
  scale_fill_manual(values = delay_colors, guide = "none") +
  scale_y_continuous(
    limits = c(0, max(wgs$time_saved_seconds / 3600) * 1.24),
    expand = expansion(mult = c(0, 0.01))
  ) +
  labs(
    title = "B  Slower pipelines return hours, not minutes",
    subtitle = "The same real WGS and cache for every load",
    x = NULL, y = "Hours returned per genome"
  ) +
  theme_vcfcache() +
  theme(legend.position = "none", panel.grid.major.y = element_blank())

matched_anchor <- annotators[
  annotators$sample == "KPGP-00319" &
    annotators$strategy == "gnomad_af_0.1",
]
vep_anchor <- matched_anchor[matched_anchor$tool == "vep", ]
fastvep_anchor <- matched_anchor[matched_anchor$tool == "fastvep", ]
lookup_overhead <- max(
  0,
  vep_anchor$cached_wall_seconds -
    (1 - vep_anchor$cache_hit_rate) * vep_anchor$uncached_wall_seconds
)
model_rows <- wgs
model_rows$predicted_cached_seconds <- lookup_overhead +
  (1 - model_rows$cache_hit_rate) * model_rows$uncached_wall_seconds
model_rows$predicted_relative_runtime <-
  model_rows$predicted_cached_seconds / model_rows$uncached_wall_seconds
model_rows$model_residual_percentage_points <- 100 * (
  model_rows$relative_runtime - model_rows$predicted_relative_runtime
)
p_model <- ggplot(
  model_rows,
  aes(x = predicted_relative_runtime, y = relative_runtime, color = pipeline_label)
) +
  geom_abline(slope = 1, intercept = 0, color = vcf_colors[["grey"]], linetype = 2) +
  geom_point(size = 4.0, stroke = 1.0) +
  scale_color_manual(values = delay_colors, guide = "none") +
  scale_x_continuous(labels = scales::percent, limits = c(0, 0.28)) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 0.28)) +
  coord_equal() +
  labs(
    title = "C  Model agreement",
    subtitle = "Predicted versus measured relative wall time",
    x = "Predicted wall time", y = "Measured wall time",
    color = NULL
  ) +
  theme_vcfcache() +
  theme(legend.position = "bottom")

hit_rate_ceiling <- 1 / (1 - unique(wgs$cache_hit_rate))
p_speedup <- ggplot(
  wgs,
  aes(x = uncached_wall_seconds / 3600, y = speedup, color = pipeline_label)
) +
  geom_hline(
    yintercept = hit_rate_ceiling, color = vcf_colors[["green"]],
    linewidth = 0.75, linetype = 2
  ) +
  geom_line(aes(group = 1), color = vcf_colors[["grey"]], linewidth = 0.7) +
  geom_point(size = 4.0, stroke = 1.0) +
  annotate(
    "label",
    x = max(wgs$uncached_wall_seconds / 3600), y = hit_rate_ceiling,
    label = sprintf("80.23%% hit-rate limit ≈ %.1f×", hit_rate_ceiling),
    hjust = 1, vjust = -0.55, size = 3.3,
    color = vcf_colors[["green"]], fill = "white", linewidth = 0
  ) +
  scale_color_manual(values = delay_colors, guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0.04, 0.04))) +
  scale_y_continuous(limits = c(4.35, 5.25), breaks = seq(4.4, 5.2, 0.2)) +
  labs(
    title = "D  Pipeline-speed dependence",
    subtitle = "Speedup saturates; absolute time saved keeps growing",
    x = "Direct annotation wall time (hours)", y = "Measured speedup"
  ) +
  theme_vcfcache()

supplement <- (p_relative | p_saved) / (p_model | p_speedup) +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(
    title = "VCFcache converts annotation complexity into time saved",
    subtitle = paste(
      "One 4.80-million-variant real WGS · bundled Zenodo gnomAD AF ≥ 10% cache ·",
      "~80% hits · light statistics · no technical repeats"
    ),
    caption = paste(
      "All six cached outputs were semantically equivalent to their pipeline-specific uncached output.",
      "Synthetic delays run only for transcript consequences submitted to VEP and emit no annotation fields."
    ),
    theme = theme(
      plot.title = element_text(face = "bold", size = 17, color = vcf_colors[["ink"]]),
      plot.subtitle = element_text(size = 11, color = vcf_colors[["grey"]]),
      plot.caption = element_text(size = 8.5, color = vcf_colors[["grey"]], hjust = 0)
    )
  )

anchor <- rbind(
  data.frame(
    example = c("fastVEP 0.3.0", "VEP 115.2 --everything"),
    source = "Matched completed campaign",
    hit_rate = c(fastvep_anchor$cache_hit_rate, vep_anchor$cache_hit_rate),
    uncached_seconds = c(
      fastvep_anchor$uncached_wall_seconds,
      vep_anchor$uncached_wall_seconds
    ),
    cached_seconds = c(
      fastvep_anchor$cached_wall_seconds,
      vep_anchor$cached_wall_seconds
    ),
    stringsAsFactors = FALSE
  ),
  data.frame(
    example = paste0("VEP + ", wgs$delay_ms, " ms/consequence"),
    source = "Virtual no-output plugin load",
    hit_rate = wgs$cache_hit_rate,
    uncached_seconds = wgs$uncached_wall_seconds,
    cached_seconds = wgs$cached_wall_seconds,
    stringsAsFactors = FALSE
  )
)
anchor$saved_seconds <- anchor$uncached_seconds - anchor$cached_seconds
anchor$relative_runtime <- anchor$cached_seconds / anchor$uncached_seconds
anchor$speedup <- anchor$uncached_seconds / anchor$cached_seconds
anchor$remaining <- anchor$relative_runtime
anchor$returned <- 1 - anchor$relative_runtime
anchor <- anchor[order(anchor$uncached_seconds), ]
anchor$example <- factor(anchor$example, levels = anchor$example)
anchor$label <- paste0(sprintf("%.1f×", anchor$speedup), "  ·  ", format_duration(anchor$saved_seconds), " back")

headline_rows <- rbind(
  transform(anchor, segment = "Time still spent", fraction = remaining),
  transform(anchor, segment = "Time returned", fraction = returned)
)
headline_rows$segment <- factor(
  headline_rows$segment,
  levels = c("Time still spent", "Time returned")
)
p_headline <- ggplot(headline_rows, aes(x = example, y = fraction, fill = segment)) +
  geom_col(width = 0.68, color = "white", linewidth = 0.8) +
  geom_text(
    data = anchor,
    aes(x = example, y = 1.025, label = label),
    inherit.aes = FALSE, hjust = 0, fontface = "bold", size = 4.1,
    color = vcf_colors[["ink"]]
  ) +
  coord_flip(clip = "off") +
  scale_fill_manual(
    values = c("Time still spent" = vcf_colors[["ink"]], "Time returned" = vcf_colors[["green"]])
  ) +
  scale_y_continuous(
    limits = c(0, 1.43), breaks = seq(0, 1, 0.25), labels = scales::percent,
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    title = "From fastVEP to an all-day pipeline: cache hits keep returning time",
    subtitle = "The same real WGS and approximately 80% reuse across the complete measured spectrum",
    x = NULL, y = "Share of the original wait", fill = NULL,
    caption = paste(
      "Bars are normalized within each measured example; labels report measured speedup and absolute time returned.",
      "Every row uses KPGP-00319 (4.80M variants) and the AF ≥ 10% cache strategy.",
      "The synthetic plugin emits no fields; it changes pipeline cost without changing annotations.",
      sep = "\n"
    )
  ) +
  theme_vcfcache(base_size = 14) +
  theme(
    panel.grid.major.y = element_blank(),
    axis.text.y = element_text(face = "bold", size = 11),
    plot.title = element_text(size = 21),
    plot.subtitle = element_text(size = 11.5),
    plot.caption = element_text(size = 8.5),
    legend.position = "top",
    plot.margin = margin(16, 80, 16, 16)
  )

p_absolute <- ggplot(anchor, aes(y = example)) +
  geom_segment(
    aes(x = cached_seconds / 60, xend = uncached_seconds / 60, yend = example),
    linewidth = 3.0, color = vcf_colors[["light_green"]], lineend = "round"
  ) +
  geom_point(aes(x = uncached_seconds / 60), size = 3.8, color = vcf_colors[["ink"]]) +
  geom_point(aes(x = cached_seconds / 60), size = 3.8, color = vcf_colors[["green"]]) +
  scale_x_log10(labels = function(x) format_duration(x * 60)) +
  labs(
    title = "Absolute waiting time removed across pipeline costs",
    subtitle = "Dark point: direct annotation · green point: with VCFcache · segment: time returned",
    x = "End-to-end wall time (log scale)", y = NULL,
    caption = "All examples use the same real WGS and approximately 80% cache hits."
  ) +
  theme_vcfcache(base_size = 12) +
  theme(panel.grid.major.y = element_blank(), axis.text.y = element_text(face = "bold"))

write_tsv(wgs, file.path(output_dir, "source", "wgs_pipeline_complexity.tsv"))
write_tsv(model_rows, file.path(output_dir, "source", "wgs_pipeline_model.tsv"))
write_tsv(anchor, file.path(output_dir, "source", "pipeline_spectrum_examples.tsv"))

panel_dir <- file.path(output_dir, "manuscript", "panels", "supplement_pipeline_complexity")
save_plot(p_relative, file.path(panel_dir, "A_relative_runtime"), 6.2, 4.6)
save_plot(p_saved, file.path(panel_dir, "B_absolute_time_saved"), 6.2, 4.6)
save_plot(p_model, file.path(panel_dir, "C_runtime_model"), 6.2, 4.8)
save_plot(p_speedup, file.path(panel_dir, "D_speedup_vs_pipeline_runtime"), 6.2, 4.8)
save_plot(
  supplement,
  file.path(output_dir, "manuscript", "supplement_pipeline_complexity"),
  12.4, 9.0
)
save_plot(
  p_headline,
  file.path(output_dir, "repository", "headliner_pipeline_spectrum"),
  12, 7.4
)
save_plot(
  p_absolute,
  file.path(output_dir, "alternatives", "pipeline_time_returned"),
  9.2, 5.8
)

writeLines(capture.output(sessionInfo()), file.path(output_dir, "R_SESSION_INFO_PIPELINE.txt"))
writeLines(
  c(
    paste("snapshot", input_dir),
    paste("rendered_at", format(Sys.time(), tz = "UTC", usetz = TRUE)),
    paste("renderer", normalizePath(file.path(script_dir, "render_pipeline_complexity.R")))
  ),
  file.path(output_dir, "RENDERED_FROM_PIPELINE.txt")
)
cat("Rendered pipeline-complexity ggplot2 figures under", output_dir, "\n")

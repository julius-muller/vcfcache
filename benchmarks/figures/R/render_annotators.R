arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) != 2) {
  stop(
    "Usage: Rscript --vanilla benchmarks/figures/R/render_annotators.R ",
    "<annotator-snapshot-directory> <output-directory>"
  )
}

full_arguments <- commandArgs(trailingOnly = FALSE)
file_argument <- full_arguments[grepl("^--file=", full_arguments)]
if (length(file_argument) != 1) {
  stop("Could not determine the render_annotators.R location")
}
script_dir <- dirname(normalizePath(sub("^--file=", "", file_argument)))
source(file.path(script_dir, "common.R"))

input_dir <- normalizePath(arguments[[1]], mustWork = TRUE)
output_dir <- arguments[[2]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

snapshot <- jsonlite::fromJSON(
  file.path(input_dir, "SNAPSHOT.json"),
  simplifyVector = FALSE
)
allowed_statuses <- c(
  "FASTVEP_FINAL_VEP_COMPARISON_PROVISIONAL",
  "ANNOTATOR_COMPARISON_FINAL",
  "ANNOTATOR_LIGHT_MATCHED_FINAL"
)
if (!snapshot$status %in% allowed_statuses) {
  stop("Unexpected annotator snapshot status: ", snapshot$status)
}
matched_light_comparison <- identical(
  snapshot$status,
  "ANNOTATOR_LIGHT_MATCHED_FINAL"
)
final_comparison <- snapshot$status %in% c(
  "ANNOTATOR_COMPARISON_FINAL",
  "ANNOTATOR_LIGHT_MATCHED_FINAL"
)

fastvep <- read_tsv(file.path(input_dir, "fastvep_external_wgs_metrics.tsv"))
combined <- read_tsv(file.path(input_dir, "annotator_external_wgs_metrics.tsv"))
calibration <- if (final_comparison) {
  read_tsv(
    file.path(
      input_dir,
      if (matched_light_comparison) {
        "vep_statistics_full_cohort_comparison.tsv"
      } else {
        "vep_statistics_calibration.tsv"
      }
    )
  )
} else {
  data.frame()
}
if (nrow(fastvep) != 156 || nrow(combined) != 312) {
  stop("Expected 156 fastVEP rows and 312 paired annotator rows")
}
if (length(unique(fastvep$sample)) != 52 ||
    !all(as_logical_semantic(fastvep$semantic_pass)) ||
    !all(grepl("^fastvep_complete_record_and_header_", fastvep$semantic_comparator))) {
  stop("The fastVEP source table is incomplete or not strictly validated")
}
expected_calibration_rows <- if (matched_light_comparison) 156 else 18
expected_calibration_samples <- if (matched_light_comparison) 52 else 6
if (final_comparison &&
    (nrow(calibration) != expected_calibration_rows ||
      length(unique(calibration$sample)) != expected_calibration_samples ||
      !all(as_logical_semantic(calibration$semantic_pass)))) {
  stop("The VEP statistics calibration is incomplete or invalid")
}

strategy_levels <- c("gnomad_af_0.1", "gnomad_af_0.01", "cohort_3_genomes")
strategy_labels <- c(
  gnomad_af_0.1 = "gnomAD AF ≥ 10%",
  gnomad_af_0.01 = "gnomAD AF ≥ 1%",
  cohort_3_genomes = "Three-genome cohort cache"
)
strategy_short <- c(
  gnomad_af_0.1 = "gnomAD\nAF ≥ 10%",
  gnomad_af_0.01 = "gnomAD\nAF ≥ 1%",
  cohort_3_genomes = "Cohort cache\n3 genomes"
)
strategy_colors <- c(
  gnomad_af_0.1 = vcf_colors[["orange"]],
  gnomad_af_0.01 = vcf_colors[["blue"]],
  cohort_3_genomes = vcf_colors[["green"]]
)
cohort_labels <- c(
  kpgp = "KPGP · GRCh38",
  sgdp = "SGDP · GRCh38",
  pgp = "PGP · GRCh37"
)
cohort_shapes <- c(kpgp = 21, sgdp = 22, pgp = 24)
tool_labels <- c(
  vep = if (matched_light_comparison) {
    "VEP 115.2\nstatistics: light"
  } else {
    "VEP 115.2\nstatistics: full"
  },
  fastvep = "fastVEP 0.3.0\nstatistics: light"
)
tool_colors <- c(vep = vcf_colors[["purple"]], fastvep = vcf_colors[["blue"]])

prepare_rows <- function(rows) {
  rows$strategy <- factor(rows$strategy, levels = strategy_levels)
  rows$strategy_label <- factor(
    strategy_short[as.character(rows$strategy)],
    levels = unname(strategy_short[strategy_levels])
  )
  rows$cohort_label <- factor(
    cohort_labels[rows$cohort],
    levels = unname(cohort_labels[c("kpgp", "sgdp", "pgp")])
  )
  rows$tool_label <- factor(
    tool_labels[rows$tool],
    levels = unname(tool_labels[c("vep", "fastvep")])
  )
  rows
}

fastvep$tool <- "fastvep"
fastvep <- prepare_rows(fastvep)
combined <- prepare_rows(combined)

paired_hit_rates <- reshape(
  combined[c("sample", "cohort", "strategy", "tool", "cache_hit_rate")],
  idvar = c("sample", "cohort", "strategy"),
  timevar = "tool",
  direction = "wide"
)
if (any(abs(paired_hit_rates$cache_hit_rate.vep - paired_hit_rates$cache_hit_rate.fastvep) > 1e-12)) {
  stop("Paired VEP and fastVEP cache hit rates differ")
}

group_summary <- function(rows, include_tool = FALSE) {
  grouping <- if (include_tool) c("tool", "cohort", "strategy") else c("cohort", "strategy")
  keys <- unique(rows[grouping])
  do.call(
    rbind,
    lapply(seq_len(nrow(keys)), function(index) {
      keep <- rep(TRUE, nrow(rows))
      for (column in grouping) {
        keep <- keep & as.character(rows[[column]]) == as.character(keys[[column]][index])
      }
      group <- rows[keep, ]
      result <- as.list(keys[index, , drop = FALSE])
      result$samples <- length(unique(group$sample))
      result$median_hit_rate <- median(group$cache_hit_rate)
      result$median_relative_runtime <- median(group$relative_runtime)
      result$median_speedup <- median(group$speedup)
      result$median_uncached_seconds <- median(group$uncached_wall_seconds)
      result$median_cached_seconds <- median(group$cached_wall_seconds)
      result$median_seconds_saved <- median(group$uncached_wall_seconds - group$cached_wall_seconds)
      as.data.frame(result, stringsAsFactors = FALSE)
    })
  )
}

fastvep_summary <- group_summary(fastvep)
annotator_summary <- group_summary(combined, include_tool = TRUE)
write_tsv(fastvep_summary, file.path(output_dir, "source", "fastvep_real_world_summary.tsv"))
write_tsv(annotator_summary, file.path(output_dir, "source", "vep_fastvep_paired_summary.tsv"))

bootstrap_keys <- unique(combined[c("tool", "strategy")])
bootstrap_summary <- do.call(
  rbind,
  lapply(seq_len(nrow(bootstrap_keys)), function(index) {
    key <- bootstrap_keys[index, ]
    group <- combined[
      combined$tool == key$tool & as.character(combined$strategy) == as.character(key$strategy),
    ]
    interval <- stratified_bootstrap_median_interval(
      group$speedup,
      group$cohort,
      seed = 20260807 + index,
      replicates = 10000
    )
    data.frame(
      tool = key$tool,
      strategy = as.character(key$strategy),
      samples = length(unique(group$sample)),
      median_speedup = median(group$speedup),
      speedup_low = interval[[1]],
      speedup_high = interval[[2]],
      stringsAsFactors = FALSE
    )
  })
)
bootstrap_summary$strategy_label <- factor(
  strategy_short[bootstrap_summary$strategy],
  levels = unname(strategy_short[strategy_levels])
)
bootstrap_summary$tool_label <- factor(
  tool_labels[bootstrap_summary$tool],
  levels = unname(tool_labels[c("vep", "fastvep")])
)
write_tsv(
  bootstrap_summary,
  file.path(output_dir, "source", "vep_fastvep_stratified_bootstrap.tsv")
)

calibration_summary <- data.frame()
if (final_comparison) {
  calibration$strategy <- factor(calibration$strategy, levels = strategy_levels)
  calibration$strategy_label <- factor(
    strategy_short[as.character(calibration$strategy)],
    levels = unname(strategy_short[strategy_levels])
  )
  calibration$cohort_label <- factor(
    cohort_labels[calibration$cohort],
    levels = unname(cohort_labels[c("kpgp", "sgdp", "pgp")])
  )
  calibration_keys <- unique(calibration["strategy"])
  calibration_summary <- do.call(
    rbind,
    lapply(seq_len(nrow(calibration_keys)), function(index) {
      strategy <- as.character(calibration_keys$strategy[index])
      group <- calibration[as.character(calibration$strategy) == strategy, ]
      interval <- stratified_bootstrap_median_interval(
        group$speedup_ratio_light_over_full,
        group$cohort,
        seed = 20260820 + index,
        replicates = 10000
      )
      data.frame(
        strategy = strategy,
        samples = length(unique(group$sample)),
        median_full_speedup = median(group$full_speedup),
        median_light_speedup = median(group$light_speedup),
        median_speedup_ratio_light_over_full = median(group$speedup_ratio_light_over_full),
        ratio_low = interval[[1]],
        ratio_high = interval[[2]],
        stringsAsFactors = FALSE
      )
    })
  )
  calibration_summary$strategy_label <- factor(
    strategy_short[calibration_summary$strategy],
    levels = unname(strategy_short[strategy_levels])
  )
  write_tsv(
    calibration_summary,
    file.path(output_dir, "source", "vep_statistics_light_calibration_summary.tsv")
  )
}

# Final standalone fastVEP figure -------------------------------------------------
p_fast_speedup <- ggplot(
  fastvep,
  aes(x = strategy_label, y = speedup, fill = strategy)
) +
  geom_hline(yintercept = 1, color = vcf_colors[["red"]], linewidth = 0.6, linetype = 2) +
  geom_boxplot(
    aes(group = strategy_label), width = 0.58, outlier.shape = NA,
    color = vcf_colors[["ink"]], fill = "white"
  ) +
  geom_point(
    aes(shape = cohort),
    position = position_jitter(width = 0.15, height = 0, seed = 20260807),
    color = vcf_colors[["ink"]], size = 2.1, stroke = 0.4, alpha = 0.70
  ) +
  scale_fill_manual(values = strategy_colors, guide = "none") +
  scale_shape_manual(values = cohort_shapes, labels = toupper(names(cohort_shapes)), name = NULL) +
  scale_y_continuous(labels = scales::label_number(suffix = "×", accuracy = 0.1)) +
  labs(
    title = "VCFcache accelerates fastVEP across all cache strategies",
    subtitle = "Each point is one independent evaluation genome; box lines show medians",
    x = NULL,
    y = "Speedup over direct fastVEP"
  ) +
  guides(shape = guide_legend(override.aes = list(fill = "white", alpha = 1, size = 3))) +
  theme_vcfcache() +
  theme(axis.text.x = element_text(size = 8.5))

direct_fastvep <- fastvep[!duplicated(fastvep[c("sample", "cohort")]), ]
direct_fastvep$configuration <- "Direct fastVEP"
direct_fastvep$wall_minutes <- direct_fastvep$uncached_wall_seconds / 60
cached_fastvep <- fastvep
cached_fastvep$configuration <- strategy_labels[as.character(cached_fastvep$strategy)]
cached_fastvep$wall_minutes <- cached_fastvep$cached_wall_seconds / 60
fast_time <- rbind(
  direct_fastvep[c("sample", "cohort", "cohort_label", "configuration", "wall_minutes")],
  cached_fastvep[c("sample", "cohort", "cohort_label", "configuration", "wall_minutes")]
)
configuration_levels <- c("Direct fastVEP", unname(strategy_labels[strategy_levels]))
fast_time$configuration <- factor(fast_time$configuration, levels = configuration_levels)
configuration_colors <- c(
  "Direct fastVEP" = vcf_colors[["grey"]],
  setNames(unname(strategy_colors[strategy_levels]), unname(strategy_labels[strategy_levels]))
)
p_fast_time <- ggplot(fast_time, aes(x = configuration, y = wall_minutes, fill = configuration)) +
  geom_boxplot(width = 0.62, outlier.shape = NA, color = vcf_colors[["ink"]]) +
  geom_point(
    position = position_jitter(width = 0.13, height = 0, seed = 20260808),
    shape = 21, color = vcf_colors[["ink"]], size = 1.55, stroke = 0.3, alpha = 0.55
  ) +
  facet_wrap(~cohort_label, nrow = 1) +
  scale_fill_manual(values = configuration_colors, guide = "none") +
  scale_y_continuous(labels = scales::label_number(suffix = " min", accuracy = 1)) +
  labs(
    title = "Minutes per whole genome fall even for a fast annotator",
    subtitle = "Identical input genomes and fastVEP configuration",
    x = NULL,
    y = "Wall time per genome"
  ) +
  theme_vcfcache() +
  theme(axis.text.x = element_text(angle = 28, hjust = 1, size = 7.5))

ideal <- data.frame(cache_hit_rate = seq(0.70, 0.98, length.out = 250))
ideal$relative_runtime <- 1 - ideal$cache_hit_rate
p_fast_model <- ggplot(
  fastvep,
  aes(x = cache_hit_rate, y = relative_runtime, fill = strategy, shape = cohort)
) +
  geom_line(
    data = ideal,
    aes(x = cache_hit_rate, y = relative_runtime),
    inherit.aes = FALSE,
    color = vcf_colors[["grey"]], linewidth = 0.8, linetype = 2
  ) +
  geom_point(color = vcf_colors[["ink"]], stroke = 0.45, size = 2.5, alpha = 0.70) +
  scale_fill_manual(values = strategy_colors, labels = strategy_labels, name = NULL) +
  scale_shape_manual(values = cohort_shapes, labels = toupper(names(cohort_shapes)), name = NULL) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, NA)) +
  labs(
    title = "Fixed cache-processing cost limits the theoretical maximum",
    subtitle = "Dashed line is the zero-overhead miss-fraction rule",
    x = "Observed cache hit rate",
    y = "Cached / direct fastVEP wall time"
  ) +
  guides(
    fill = guide_legend(
      order = 1, nrow = 1,
      override.aes = list(shape = 21, color = vcf_colors[["ink"]], size = 3.5, alpha = 1)
    ),
    shape = guide_legend(
      order = 2, nrow = 1,
      override.aes = list(fill = "white", color = vcf_colors[["ink"]], size = 3.5, alpha = 1)
    )
  ) +
  theme_vcfcache() +
  theme(legend.text = element_text(size = 8))

fast_combined <- (p_fast_speedup | p_fast_model) / p_fast_time +
  plot_layout(heights = c(1, 0.92)) +
  plot_annotation(
    title = "VCFcache remains useful with fastVEP",
    subtitle = "52 real-world WGS evaluation genomes · one run per genome and condition",
    caption = paste(
      "fastVEP 0.3.0 · VCFcache --statistics light · KPGP, SGDP and PGP.",
      "All 156 cached outputs passed strict complete-record and relevant-header comparison against direct fastVEP.",
      "Bundled-cache rows use locally built fastVEP caches derived from the corresponding VCFcache Zenodo blueprints; cohort caches use three separate build genomes.",
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

fast_panel_dir <- file.path(output_dir, "manuscript", "panels", "fastvep_real_world_wgs")
save_plot(p_fast_speedup, file.path(fast_panel_dir, "A_speedup"), 7.4, 5.2)
save_plot(p_fast_model, file.path(fast_panel_dir, "B_hit_rate_and_overhead"), 7.4, 5.2)
save_plot(p_fast_time, file.path(fast_panel_dir, "C_wall_time"), 14.8, 5.0)
save_plot(fast_combined, file.path(output_dir, "manuscript", "fastvep_real_world_wgs"), 15, 10.5)

# Provisional paired VEP/fastVEP figure ------------------------------------------
cohort_points <- annotator_summary
cohort_points$strategy_label <- factor(
  strategy_short[cohort_points$strategy],
  levels = unname(strategy_short[strategy_levels])
)
cohort_points$tool_label <- factor(
  tool_labels[cohort_points$tool],
  levels = unname(tool_labels[c("vep", "fastvep")])
)
p_compare_speedup <- ggplot(
  bootstrap_summary,
  aes(x = strategy_label, y = median_speedup, color = tool, group = tool)
) +
  geom_hline(yintercept = 1, color = vcf_colors[["red"]], linewidth = 0.6, linetype = 2) +
  geom_errorbar(
    aes(ymin = speedup_low, ymax = speedup_high),
    position = position_dodge(width = 0.44), width = 0.13, linewidth = 0.75
  ) +
  geom_point(position = position_dodge(width = 0.44), size = 3.3) +
  geom_point(
    data = cohort_points,
    aes(x = strategy_label, y = median_speedup, color = tool, shape = cohort),
    position = position_dodge(width = 0.44),
    inherit.aes = FALSE, size = 2.1, stroke = 0.5, alpha = 0.58
  ) +
  scale_color_manual(values = tool_colors, labels = tool_labels, name = NULL) +
  scale_shape_manual(values = cohort_shapes, labels = toupper(names(cohort_shapes)), name = NULL) +
  scale_y_continuous(labels = scales::label_number(suffix = "×", accuracy = 0.1)) +
  labs(
    title = "Caching helps both annotators",
    subtitle = "Large points: stratified median with 95% bootstrap interval; small symbols: cohort medians",
    x = NULL,
    y = "Speedup over direct annotation"
  ) +
  theme_vcfcache() +
  theme(axis.text.x = element_text(size = 8.5), legend.text = element_text(size = 8))

p_compare_model <- ggplot(
  combined,
  aes(x = cache_hit_rate, y = relative_runtime, fill = strategy, shape = cohort)
) +
  geom_line(
    data = ideal,
    aes(x = cache_hit_rate, y = relative_runtime),
    inherit.aes = FALSE,
    color = vcf_colors[["grey"]], linewidth = 0.8, linetype = 2
  ) +
  geom_point(color = vcf_colors[["ink"]], stroke = 0.4, size = 2.0, alpha = 0.60) +
  facet_wrap(~tool_label, nrow = 1) +
  scale_fill_manual(values = strategy_colors, labels = strategy_labels, name = NULL) +
  scale_shape_manual(values = cohort_shapes, labels = toupper(names(cohort_shapes)), name = NULL) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 10), limits = c(0, NA)) +
  labs(
    title = "The same hits remove more work from the slower pipeline",
    subtitle = "Dashed line: zero-overhead miss-fraction rule; axes are identical",
    x = "Observed cache hit rate",
    y = "Cached / direct wall time"
  ) +
  guides(
    fill = guide_legend(
      order = 1, nrow = 1,
      override.aes = list(shape = 21, color = vcf_colors[["ink"]], size = 3.2, alpha = 1)
    ),
    shape = guide_legend(
      order = 2, nrow = 1,
      override.aes = list(fill = "white", color = vcf_colors[["ink"]], size = 3.2, alpha = 1)
    )
  ) +
  theme_vcfcache() +
  theme(legend.text = element_text(size = 7.8))

time_pair <- annotator_summary
time_pair$direct_minutes <- time_pair$median_uncached_seconds / 60
time_pair$cached_minutes <- time_pair$median_cached_seconds / 60
time_pair$strategy_label <- factor(
  strategy_short[time_pair$strategy],
  levels = unname(strategy_short[strategy_levels])
)
time_pair$cohort_label <- factor(
  cohort_labels[time_pair$cohort],
  levels = unname(cohort_labels[c("kpgp", "sgdp", "pgp")])
)
time_pair$tool_label <- factor(
  tool_labels[time_pair$tool],
  levels = unname(tool_labels[c("vep", "fastvep")])
)
p_compare_time <- ggplot(time_pair, aes(y = strategy_label, color = tool)) +
  geom_segment(
    aes(x = cached_minutes, xend = direct_minutes, yend = strategy_label),
    linewidth = 1.0, alpha = 0.65,
    position = position_dodge(width = 0.46)
  ) +
  geom_point(
    aes(x = direct_minutes), shape = 21, fill = "white", size = 2.8,
    stroke = 0.9, position = position_dodge(width = 0.46)
  ) +
  geom_point(
    aes(x = cached_minutes), shape = 19, size = 2.8,
    position = position_dodge(width = 0.46)
  ) +
  facet_wrap(~cohort_label, nrow = 1) +
  scale_color_manual(values = tool_colors, labels = tool_labels, name = NULL) +
  scale_x_log10(
    breaks = c(1, 2, 5, 10, 30, 60, 120, 240),
    labels = scales::label_number(suffix = " min")
  ) +
  labs(
    title = "Absolute wall time spans minutes to hours",
    subtitle = "Open point: direct; filled point: cached; line: time returned per genome",
    x = "Median wall time per genome (log scale)",
    y = NULL
  ) +
  theme_vcfcache() +
  theme(axis.text.y = element_text(size = 8), legend.text = element_text(size = 8))

p_calibration <- NULL
if (final_comparison) {
  calibration$relative_speedup_change <- calibration$speedup_ratio_light_over_full - 1
  calibration_summary$relative_speedup_change <-
    calibration_summary$median_speedup_ratio_light_over_full - 1
  calibration_summary$change_low <- calibration_summary$ratio_low - 1
  calibration_summary$change_high <- calibration_summary$ratio_high - 1
  calibration_summary$change_label <- scales::percent(
    calibration_summary$relative_speedup_change,
    accuracy = 0.1,
    prefix = "+"
  )
  calibration_summary$label_y <- vapply(
    calibration_summary$strategy,
    function(strategy_name) {
      max(
        calibration$relative_speedup_change[
          as.character(calibration$strategy) == as.character(strategy_name)
        ]
      ) + 0.025
    },
    numeric(1)
  )
  p_calibration <- ggplot(
    calibration,
    aes(
      x = strategy_label,
      y = relative_speedup_change,
      fill = strategy,
      shape = cohort
    )
  ) +
    geom_hline(yintercept = 0, color = vcf_colors[["grey"]], linewidth = 0.7) +
    geom_point(
      position = position_jitter(width = 0.10, height = 0, seed = 20260809),
      color = vcf_colors[["ink"]], size = 2.8, stroke = 0.45, alpha = 0.72
    ) +
    geom_errorbar(
      data = calibration_summary,
      aes(
        x = strategy_label,
        ymin = change_low,
        ymax = change_high
      ),
      inherit.aes = FALSE,
      width = 0.12,
      color = vcf_colors[["ink"]],
      linewidth = 0.75
    ) +
    geom_point(
      data = calibration_summary,
      aes(x = strategy_label, y = relative_speedup_change),
      inherit.aes = FALSE,
      color = vcf_colors[["ink"]], fill = "white", shape = 21, size = 3.5, stroke = 0.9
    ) +
    geom_text(
      data = calibration_summary,
      aes(x = strategy_label, y = label_y, label = change_label),
      inherit.aes = FALSE,
      color = vcf_colors[["ink"]], fontface = "bold", size = 3.3
    ) +
    scale_fill_manual(values = strategy_colors, guide = "none") +
    scale_shape_manual(values = cohort_shapes, labels = toupper(names(cohort_shapes)), name = NULL) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1), expand = expansion(mult = c(0.08, 0.20))) +
    labs(
      title = if (matched_light_comparison) {
        "The full-cohort rerun confirms that the legacy rescan was conservative"
      } else {
        "Light statistics make the historical VEP estimate slightly conservative"
      },
      subtitle = if (matched_light_comparison) {
        "All 52 genomes paired by sample and strategy; points compare light/full speedup, labels show medians"
      } else {
        "Six paired genomes (two per cohort); points compare light/full speedup, labels show medians"
      },
      x = NULL,
      y = "Change in measured VEP speedup"
    ) +
    guides(shape = guide_legend(override.aes = list(fill = "white", alpha = 1, size = 3))) +
    theme_vcfcache() +
    theme(axis.text.x = element_text(size = 8.5))
}

comparison_body <- if (final_comparison) {
  (p_compare_speedup | p_compare_model) / p_compare_time / p_calibration +
    plot_layout(heights = c(1, 0.86, 0.72))
} else {
  (p_compare_speedup | p_compare_model) / p_compare_time +
    plot_layout(heights = c(1, 0.92))
}
comparison_combined <- comparison_body +
  plot_annotation(
    title = "VCFcache reduces runtime for both VEP and fastVEP",
    subtitle = if (matched_light_comparison) {
      "Exactly paired genomes, variants, cache blueprints and light-statistics timing mode"
    } else if (final_comparison) {
      "Exactly paired genomes, variants and cache blueprints · completed statistics sensitivity calibration"
    } else {
      "Exactly paired genomes, variants and cache blueprints · comparison draft"
    },
    caption = if (matched_light_comparison) {
      paste(
        "VEP 115.2 and fastVEP 0.3.0 were both timed with VCFcache --statistics light across all 52 genomes.",
        "Panel D pairs the new VEP light measurements with the historical full-rescan values; no observation was rescaled.",
        "Each genome/tool/strategy combination was run once; stratified bootstrap intervals treat genomes as inferential units.",
        sep = "\n"
      )
    } else if (final_comparison) {
      paste(
        "VEP 115.2 values are the measured --statistics full results; fastVEP 0.3.0 used --statistics light.",
        "The six-genome paired calibration is shown directly and is not used to rescale the 52 original VEP observations.",
        "Each genome/tool/strategy combination was run once; stratified bootstrap intervals treat genomes as inferential units.",
        sep = "\n"
      )
    } else {
      paste(
        "PROVISIONAL CROSS-ANNOTATOR VIEW: VEP 115.2 was timed with --statistics full; fastVEP 0.3.0 with --statistics light.",
        "The standalone fastVEP results are final. Replace the VEP timings after the six-genome light-statistics calibration is complete.",
        "Each genome/tool/strategy combination was run once; genomes, stratified by cohort, are the inferential units.",
        sep = "\n"
      )
    },
    tag_levels = "A",
    theme = theme(
      plot.title = element_text(size = 18, face = "bold", color = vcf_colors[["ink"]]),
      plot.subtitle = element_text(size = 11, color = vcf_colors[["grey"]]),
      plot.caption = element_text(
        size = 8,
        color = if (final_comparison) vcf_colors[["grey"]] else vcf_colors[["red"]],
        hjust = 0
      ),
      plot.tag = element_text(face = "bold", color = vcf_colors[["blue"]])
    )
  )

comparison_name <- if (final_comparison) {
  if (matched_light_comparison) {
    "vep_fastvep_light_matched_final"
  } else {
    "vep_fastvep_impact_final"
  }
} else {
  "vep_fastvep_impact_provisional"
}
comparison_panel_dir <- file.path(output_dir, "manuscript", "panels", comparison_name)
save_plot(p_compare_speedup, file.path(comparison_panel_dir, "A_speedup"), 7.4, 5.2)
save_plot(p_compare_model, file.path(comparison_panel_dir, "B_hit_rate_and_runtime"), 7.4, 5.2)
save_plot(p_compare_time, file.path(comparison_panel_dir, "C_absolute_wall_time"), 14.8, 5.0)
if (final_comparison) {
  save_plot(p_calibration, file.path(comparison_panel_dir, "D_statistics_calibration"), 14.8, 4.5)
}
save_plot(
  comparison_combined,
  file.path(output_dir, "manuscript", comparison_name),
  15, if (final_comparison) 14.5 else 10.5
)

# Repository headline alternatives ------------------------------------------------
headline_strategy <- "gnomad_af_0.01"
headline <- do.call(
  rbind,
  lapply(c("vep", "fastvep"), function(tool_name) {
    group <- combined[
      combined$tool == tool_name & as.character(combined$strategy) == headline_strategy,
    ]
    data.frame(
      tool = tool_name,
      tool_display = ifelse(tool_name == "vep", "VEP", "fastVEP"),
      relative_runtime = median(group$relative_runtime),
      speedup = median(group$speedup),
      samples = length(unique(group$sample)),
      stringsAsFactors = FALSE
    )
  })
)
headline$tool_display <- factor(headline$tool_display, levels = c("fastVEP", "VEP"))
headline$remaining_label <- scales::percent(headline$relative_runtime, accuracy = 1)
headline$saved_label <- scales::percent(1 - headline$relative_runtime, accuracy = 1)
headline$speedup_label <- sprintf("%.1f× faster", headline$speedup)
write_tsv(headline, file.path(output_dir, "source", "repository_headline_before_after.tsv"))

p_headline_before_after <- ggplot(headline, aes(x = tool_display)) +
  geom_col(aes(y = 1), width = 0.56, fill = vcf_colors[["light_grey"]], color = vcf_colors[["ink"]]) +
  geom_col(
    aes(y = relative_runtime, fill = tool),
    width = 0.56, color = vcf_colors[["ink"]], linewidth = 0.65
  ) +
  geom_text(
    aes(y = relative_runtime / 2, label = paste0(remaining_label, "\nleft")),
    color = "white", fontface = "bold", size = 5.0, lineheight = 0.9
  ) +
  geom_text(
    aes(y = relative_runtime + (1 - relative_runtime) / 2, label = paste0(saved_label, "\nskipped")),
    color = vcf_colors[["ink"]], fontface = "bold", size = 4.3, lineheight = 0.9
  ) +
  geom_text(
    aes(y = 1.08, label = speedup_label, color = tool),
    fontface = "bold", size = 6.0
  ) +
  coord_flip(clip = "off") +
  scale_fill_manual(values = tool_colors, guide = "none") +
  scale_color_manual(values = tool_colors, guide = "none") +
  scale_y_continuous(limits = c(0, 1.16), breaks = c(0, 0.25, 0.5, 0.75, 1), labels = scales::percent) +
  labs(
    title = "Skip the annotation work you have already done",
    subtitle = "Real-world WGS · bundled gnomAD AF ≥ 1% cache · median across 52 genomes",
    x = NULL,
    y = "Relative wall time",
    caption = "✓ Same annotations after semantic validation   ·   One run per genome and condition"
  ) +
  theme_vcfcache(base_size = 15) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(color = "white", linewidth = 1.1),
    axis.text.y = element_text(face = "bold", size = 16),
    axis.text.x = element_blank(),
    axis.ticks = element_blank(),
    axis.title.x = element_blank(),
    plot.title = element_text(size = 24),
    plot.subtitle = element_text(size = 13),
    plot.caption = element_text(size = 11, face = "bold", color = vcf_colors[["green"]]),
    plot.margin = margin(18, 35, 18, 18)
  )

empirical_overhead <- do.call(
  rbind,
  lapply(c("vep", "fastvep"), function(tool_name) {
    group <- combined[combined$tool == tool_name, ]
    data.frame(
      tool = tool_name,
      overhead_fraction = median(group$relative_runtime - (1 - group$cache_hit_rate)),
      stringsAsFactors = FALSE
    )
  })
)
scenario_definitions <- data.frame(
  scenario = c("Distant WGS\n50% hits", "Related WGS\n80% hits", "WES\n90% hits", "Very high reuse\n95% hits"),
  hit_rate = c(0.50, 0.80, 0.90, 0.95),
  stringsAsFactors = FALSE
)
scenario_rows <- merge(empirical_overhead, scenario_definitions, all = TRUE)
scenario_rows$relative_runtime <- pmin(
  1,
  scenario_rows$overhead_fraction + 1 - scenario_rows$hit_rate
)
scenario_rows$speedup <- 1 / scenario_rows$relative_runtime
scenario_rows$tool_display <- ifelse(
  scenario_rows$tool == "vep",
  "Slower pipeline\n(VEP-like overhead)",
  "Very fast pipeline\n(fastVEP-like overhead)"
)
scenario_rows$scenario <- factor(
  scenario_rows$scenario,
  levels = scenario_definitions$scenario
)
scenario_rows$tool_display <- factor(
  scenario_rows$tool_display,
  levels = c("Very fast pipeline\n(fastVEP-like overhead)", "Slower pipeline\n(VEP-like overhead)")
)
scenario_rows$speedup_label <- sprintf("%.1f×", scenario_rows$speedup)
scenario_rows$remaining_label <- paste0(
  scales::percent(scenario_rows$relative_runtime, accuracy = 1),
  " time left"
)
write_tsv(scenario_rows, file.path(output_dir, "source", "repository_headline_scenarios.tsv"))

p_headline_scenarios <- ggplot(
  scenario_rows,
  aes(x = scenario, y = tool_display, fill = speedup)
) +
  geom_tile(color = "white", linewidth = 5) +
  geom_text(aes(label = speedup_label), color = "white", fontface = "bold", size = 8) +
  geom_text(
    aes(label = remaining_label),
    color = "white", size = 3.5, vjust = 2.4
  ) +
  scale_fill_gradientn(
    colors = c(vcf_colors[["grey"]], vcf_colors[["blue"]], vcf_colors[["green"]]),
    guide = "none"
  ) +
  labs(
    title = "More reuse means less waiting—whatever your pipeline",
    subtitle = "Choose the column closest to your samples; bigger numbers are better",
    x = NULL,
    y = NULL,
    caption = paste(
      "Illustrative runtime model using median overhead measured in 52 real WGS genomes.",
      "Actual speedup depends on your cache hit rate and annotation pipeline. Cache build cost excluded.",
      sep = "\n"
    )
  ) +
  theme_void(base_size = 15) +
  theme(
    plot.background = element_rect(fill = vcf_colors[["paper"]], color = NA),
    plot.title = element_text(size = 24, face = "bold", color = vcf_colors[["ink"]], margin = margin(b = 8)),
    plot.subtitle = element_text(size = 13, color = vcf_colors[["grey"]], margin = margin(b = 16)),
    plot.caption = element_text(size = 9, color = vcf_colors[["grey"]], hjust = 0, margin = margin(t = 14)),
    axis.text.x = element_text(size = 12, face = "bold", color = vcf_colors[["ink"]], margin = margin(t = 8)),
    axis.text.y = element_text(size = 12, face = "bold", color = vcf_colors[["ink"]], margin = margin(r = 8)),
    plot.margin = margin(18, 18, 18, 18)
  )

workflow_nodes <- data.frame(
  x = c(1.4, 4.5, 4.5, 8.0, 11.5),
  y = c(5.0, 5.0, 8.0, 5.0, 5.0),
  label = c(
    "New sample\nVCF",
    "VCFcache\nlookup",
    "Bundled gnomAD\nor your cohort",
    "VEP / fastVEP\nannotate misses",
    "Complete\nannotated VCF  ✓"
  ),
  kind = c("input", "cache", "source", "annotator", "output"),
  stringsAsFactors = FALSE
)
workflow_colors <- c(
  input = "#FFF4E8",
  cache = "#D9EAF2",
  source = "#E4F4E8",
  annotator = "#EFE8F7",
  output = "#DDF3E5"
)
workflow_arrows <- data.frame(
  x = c(2.25, 5.15, 5.35, 8.90),
  y = c(5.0, 7.25, 5.0, 5.0),
  xend = c(3.55, 4.75, 7.05, 10.55),
  yend = c(5.0, 5.8, 5.0, 5.0),
  curvature = c(-0.04, 0.18, 0.04, -0.04)
)
workflow_claims <- data.frame(
  x = c(3.0, 7.1, 10.5),
  y = c(1.55, 1.55, 1.55),
  label = c("50–95%\nrealistic reuse", "1.8–8.6×\nmedian speedup", "Same annotations\nafter validation"),
  color = c(vcf_colors[["orange"]], vcf_colors[["blue"]], vcf_colors[["green"]]),
  stringsAsFactors = FALSE
)
write_tsv(workflow_nodes, file.path(output_dir, "source", "repository_headline_workflow_nodes.tsv"))

p_headline_workflow <- ggplot() +
  geom_curve(
    data = workflow_arrows,
    aes(x = x, y = y, xend = xend, yend = yend),
    curvature = 0.08,
    linewidth = 1.2,
    color = vcf_colors[["grey"]],
    arrow = grid::arrow(length = grid::unit(0.18, "inches"), type = "closed")
  ) +
  geom_label(
    data = workflow_nodes,
    aes(x = x, y = y, label = label, fill = kind),
    color = vcf_colors[["ink"]], fontface = "bold", size = 4.6,
    linewidth = 0.8, label.padding = grid::unit(0.38, "lines"),
    label.r = grid::unit(0.18, "lines"), lineheight = 0.95
  ) +
  geom_label(
    data = workflow_claims,
    aes(x = x, y = y, label = label, color = color),
    fill = "white", fontface = "bold", size = 4.6,
    linewidth = 0.8, label.padding = grid::unit(0.32, "lines"),
    label.r = grid::unit(0.18, "lines"), lineheight = 0.95,
    show.legend = FALSE
  ) +
  scale_fill_manual(values = workflow_colors, guide = "none") +
  scale_color_identity() +
  coord_cartesian(xlim = c(0.3, 12.8), ylim = c(0.4, 9.0), clip = "off") +
  labs(
    title = "Annotate what is new. Reuse what is known.",
    subtitle = "VCFcache sits in front of your annotator and sends only cache misses through the expensive path",
    caption = "Observed range across real-world WGS strategies and annotators · exact impact depends on cache coverage and pipeline cost"
  ) +
  theme_void(base_size = 15) +
  theme(
    plot.background = element_rect(fill = vcf_colors[["paper"]], color = NA),
    plot.title = element_text(size = 25, face = "bold", color = vcf_colors[["ink"]], margin = margin(b = 8)),
    plot.subtitle = element_text(size = 13, color = vcf_colors[["grey"]], margin = margin(b = 10)),
    plot.caption = element_text(size = 9, color = vcf_colors[["grey"]], hjust = 0),
    plot.margin = margin(20, 20, 16, 20)
  )

save_plot(
  p_headline_before_after,
  file.path(output_dir, "repository", "headliner_before_after"),
  12, 6.75
)
save_plot(
  p_headline_scenarios,
  file.path(output_dir, "repository", "headliner_scenario_map"),
  12, 6.75
)
save_plot(
  p_headline_workflow,
  file.path(output_dir, "repository", "headliner_workflow_doodle"),
  12, 6.75
)

session <- capture.output(sessionInfo())
writeLines(session, file.path(output_dir, "R_SESSION_INFO_ANNOTATORS.txt"))
writeLines(
  c(
    paste("snapshot", input_dir),
    paste("rendered_at", format(Sys.time(), tz = "UTC", usetz = TRUE)),
    paste("renderer", normalizePath(file.path(script_dir, "render_annotators.R"), mustWork = FALSE))
  ),
  file.path(output_dir, "RENDERED_FROM_ANNOTATORS.txt")
)

cat("Rendered fastVEP and paired annotator ggplot2 figures under", output_dir, "\n")

render_alternative_figures <- function(input_dir, output_dir, snapshot) {
  generic_metrics <- file.path(input_dir, "external_wgs_metrics.tsv")
  legacy_metrics <- file.path(input_dir, "external_wgs_preliminary_metrics.tsv")
  external <- read_tsv(if (file.exists(generic_metrics)) generic_metrics else legacy_metrics)
  if (!all(external$validation_status == "semantically_validated")) {
    warning("Alternative external-WGS figures include non-final rows")
  }

  strategy_levels <- c("gnomad_af_0.1", "gnomad_af_0.01", "cohort_3_genomes")
  strategy_labels <- c(
    gnomad_af_0.1 = "gnomAD ≥10%",
    gnomad_af_0.01 = "gnomAD ≥1%",
    cohort_3_genomes = "Cohort: 3 genomes"
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
  external$strategy_label <- factor(
    strategy_labels[external$strategy],
    levels = unname(strategy_labels[strategy_levels])
  )
  external$cohort_label <- factor(
    cohort_labels[external$cohort],
    levels = unname(cohort_labels[c("kpgp", "sgdp", "pgp")])
  )
  external$hours_saved <- (
    external$uncached_wall_seconds - external$cached_wall_seconds
  ) / 3600

  keys <- unique(external[c("cohort", "strategy")])
  medians <- do.call(
    rbind,
    lapply(seq_len(nrow(keys)), function(index) {
      key <- keys[index, ]
      group <- external[
        external$cohort == key$cohort &
          external$strategy == key$strategy,
      ]
      data.frame(
        cohort = key$cohort,
        cohort_label = cohort_labels[[key$cohort]],
        strategy = key$strategy,
        strategy_label = strategy_labels[[key$strategy]],
        samples = length(unique(group$sample)),
        median_uncached_hours = median(group$uncached_wall_seconds) / 3600,
        median_cached_hours = median(group$cached_wall_seconds) / 3600,
        median_hours_saved = median(group$hours_saved),
        median_hit_rate = median(group$cache_hit_rate),
        median_speedup = median(group$speedup),
        stringsAsFactors = FALSE
      )
    })
  )
  medians$strategy_label <- factor(
    medians$strategy_label,
    levels = unname(strategy_labels[strategy_levels])
  )
  medians$cohort_label <- factor(
    medians$cohort_label,
    levels = unname(cohort_labels[c("kpgp", "sgdp", "pgp")])
  )
  write_tsv(
    medians,
    file.path(output_dir, "source", "alternative_external_wgs_medians.tsv")
  )

  p_before_after <- ggplot(medians, aes(y = strategy_label)) +
    geom_segment(
      aes(
        x = median_cached_hours,
        xend = median_uncached_hours,
        yend = strategy_label
      ),
      linewidth = 2.2,
      color = vcf_colors[["light_grey"]],
      lineend = "round"
    ) +
    geom_point(
      aes(x = median_uncached_hours),
      size = 3.8,
      shape = 21,
      fill = "white",
      color = vcf_colors[["grey"]],
      stroke = 1
    ) +
    geom_point(
      aes(x = median_cached_hours, fill = strategy),
      size = 4.2,
      shape = 21,
      color = "white",
      stroke = 0.8
    ) +
    geom_text(
      aes(
        x = median_cached_hours,
        label = sprintf("%.1f×", median_speedup)
      ),
      nudge_y = 0.27,
      size = 3.2,
      fontface = "bold",
      color = vcf_colors[["ink"]]
    ) +
    facet_wrap(~cohort_label, nrow = 1) +
    scale_fill_manual(values = strategy_colors, guide = "none") +
    scale_x_log10(
      labels = scales::label_number(accuracy = 0.1, suffix = " h")
    ) +
    labs(
      title = "A workday becomes a much shorter wait",
      subtitle = "Median wall time per real-world genome · white = direct VEP · colored = VCFcache",
      x = "Median wall time per genome (log scale)",
      y = NULL
    ) +
    theme_vcfcache() +
    theme(axis.text.y = element_text(size = 8.5))

  p_hit_rate <- ggplot(
    external,
    aes(x = strategy_label, y = cache_hit_rate, color = strategy)
  ) +
    geom_boxplot(
      width = 0.58,
      outlier.shape = NA,
      color = vcf_colors[["ink"]],
      fill = "white"
    ) +
    geom_point(
      position = position_jitter(width = 0.13, height = 0, seed = 20260806),
      alpha = 0.62,
      size = 1.55
    ) +
    facet_wrap(~cohort_label, nrow = 1) +
    scale_color_manual(values = strategy_colors, guide = "none") +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    labs(
      title = "Cache reuse across external genomes",
      subtitle = "Held out from cohort-cache construction; no documented gnomAD project overlap",
      x = NULL,
      y = "Variants found in cache"
    ) +
    theme_vcfcache() +
    theme(axis.text.x = element_text(angle = 18, hjust = 1, size = 8))

  p_saved <- ggplot(
    external,
    aes(x = strategy_label, y = hours_saved, color = strategy)
  ) +
    geom_hline(yintercept = 0, color = vcf_colors[["red"]], linetype = 2) +
    geom_boxplot(
      width = 0.58,
      outlier.shape = NA,
      color = vcf_colors[["ink"]],
      fill = "white"
    ) +
    geom_point(
      position = position_jitter(width = 0.13, height = 0, seed = 20260807),
      alpha = 0.62,
      size = 1.55
    ) +
    facet_wrap(~cohort_label, nrow = 1) +
    scale_color_manual(values = strategy_colors, guide = "none") +
    labs(
      title = "Hours returned for every genome",
      subtitle = "Positive values are measured wall time saved",
      x = NULL,
      y = "Hours saved per genome"
    ) +
    theme_vcfcache() +
    theme(axis.text.x = element_text(angle = 18, hjust = 1, size = 8))

  decision_figure <- p_before_after / (p_hit_rate | p_saved) +
    plot_layout(heights = c(0.95, 1.05)) +
    plot_annotation(
      title = "Alternative view: real-world WGS benefit at a glance",
      subtitle = "One run per sample and condition · bundled Zenodo caches versus cohort-derived reuse",
      caption = paste(
        preliminary_caption(snapshot),
        "Cached outputs were compared with direct VEP from the same sample and recipe.",
        "No documented project overlap with gnomAD; individual absence from undisclosed contributors cannot be proven.",
        sep = "\n"
      ),
      tag_levels = "A",
      theme = theme(
        plot.title = element_text(size = 18, face = "bold", color = vcf_colors[["ink"]]),
        plot.subtitle = element_text(size = 11, color = vcf_colors[["grey"]]),
        plot.caption = element_text(size = 8, color = vcf_colors[["grey"]], hjust = 0)
      )
    )

  decision_panels <- file.path(
    output_dir, "alternatives", "panels", "external_wgs_decision"
  )
  save_plot(p_before_after, file.path(decision_panels, "A_before_after"), 14.8, 4.7)
  save_plot(p_hit_rate, file.path(decision_panels, "B_hit_rate_distribution"), 7.4, 5.0)
  save_plot(p_saved, file.path(decision_panels, "C_hours_saved"), 7.4, 5.0)
  save_plot(
    decision_figure,
    file.path(output_dir, "alternatives", "external_wgs_decision"),
    15,
    10
  )

  primary <- read_tsv(file.path(input_dir, "primary_wgs_metrics.tsv"))
  assay <- read_tsv(file.path(input_dir, "assay_metrics.tsv"))
  paired <- rbind(primary, assay)
  assay_labels <- c(
    panel = "1000G-derived panel",
    wes = "1000G-derived WES",
    wgs = "1000G WGS"
  )
  paired$evidence_group <- assay_labels[paired$assay]
  paired$source_type <- "Source-overlap calibration"
  external$evidence_group <- paste0(toupper(external$cohort), " WGS")
  external$source_type <- "External WGS"
  diagnostics <- rbind(
    paired[c(
      "sample", "evidence_group", "source_type", "input_records",
      "cache_hit_rate", "cached_wall_seconds", "uncached_wall_seconds",
      "relative_runtime", "speedup"
    )],
    external[c(
      "sample", "evidence_group", "source_type", "input_records",
      "cache_hit_rate", "cached_wall_seconds", "uncached_wall_seconds",
      "relative_runtime", "speedup"
    )]
  )
  diagnostics$miss_fraction <- 1 - diagnostics$cache_hit_rate
  diagnostics$inferred_overhead_seconds <- diagnostics$cached_wall_seconds -
    diagnostics$miss_fraction * diagnostics$uncached_wall_seconds
  diagnostics$ideal_speedup <- 1 / diagnostics$miss_fraction
  diagnostics$speedup_efficiency <- diagnostics$speedup / diagnostics$ideal_speedup
  write_tsv(
    diagnostics,
    file.path(output_dir, "source", "alternative_runtime_diagnostics.tsv")
  )

  group_levels <- c(
    "1000G-derived panel", "1000G-derived WES", "1000G WGS",
    "KPGP WGS", "SGDP WGS", "PGP WGS"
  )
  diagnostic_colors <- c(
    "1000G-derived panel" = vcf_colors[["green"]],
    "1000G-derived WES" = vcf_colors[["orange"]],
    "1000G WGS" = vcf_colors[["blue"]],
    "KPGP WGS" = vcf_colors[["purple"]],
    "SGDP WGS" = "#00838F",
    "PGP WGS" = vcf_colors[["red"]]
  )
  diagnostics$evidence_group <- factor(
    diagnostics$evidence_group,
    levels = group_levels
  )
  ideal <- data.frame(miss_fraction = seq(0, 0.55, length.out = 300))
  ideal$relative_runtime <- ideal$miss_fraction

  p_model <- ggplot(
    diagnostics,
    aes(
      x = miss_fraction,
      y = relative_runtime,
      color = evidence_group,
      shape = source_type
    )
  ) +
    geom_line(
      data = ideal,
      aes(x = miss_fraction, y = relative_runtime),
      inherit.aes = FALSE,
      color = vcf_colors[["grey"]],
      linewidth = 0.9,
      linetype = 2
    ) +
    geom_point(alpha = 0.62, size = 2) +
    scale_color_manual(values = diagnostic_colors, name = NULL, drop = FALSE) +
    scale_shape_manual(
      values = c("Source-overlap calibration" = 16, "External WGS" = 17),
      name = NULL
    ) +
    scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    labs(
      title = "The miss fraction predicts work remaining",
      subtitle = "Dashed line is the zero-overhead runtime model",
      x = "Variants requiring annotation",
      y = "Cached / direct wall time"
    ) +
    theme_vcfcache()

  p_overhead <- ggplot(
    diagnostics,
    aes(
      x = uncached_wall_seconds / 3600,
      y = inferred_overhead_seconds / 60,
      color = evidence_group,
      shape = source_type
    )
  ) +
    geom_hline(yintercept = 0, color = vcf_colors[["grey"]], linetype = 2) +
    geom_point(alpha = 0.62, size = 2) +
    scale_color_manual(values = diagnostic_colors, guide = "none", drop = FALSE) +
    scale_shape_manual(
      values = c("Source-overlap calibration" = 16, "External WGS" = 17),
      guide = "none"
    ) +
    scale_x_log10(
      labels = scales::label_number(accuracy = 0.1, suffix = " h")
    ) +
    labs(
      title = "Lookup overhead matters most for short pipelines",
      subtitle = "Inferred overhead = measured cached time − miss fraction × direct time",
      x = "Direct annotation time (log scale)",
      y = "Inferred overhead (minutes)"
    ) +
    theme_vcfcache()

  p_efficiency <- ggplot(
    diagnostics,
    aes(x = evidence_group, y = speedup_efficiency, color = evidence_group)
  ) +
    geom_hline(yintercept = 1, color = vcf_colors[["grey"]], linetype = 2) +
    geom_boxplot(
      width = 0.58,
      outlier.shape = NA,
      color = vcf_colors[["ink"]],
      fill = "white"
    ) +
    geom_point(
      position = position_jitter(width = 0.13, height = 0, seed = 20260808),
      alpha = 0.5,
      size = 1.5
    ) +
    scale_color_manual(values = diagnostic_colors, guide = "none", drop = FALSE) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    labs(
      title = "How much of the ideal speedup is realized?",
      subtitle = "100% means zero lookup and preprocessing overhead",
      x = NULL,
      y = "Observed / zero-overhead speedup"
    ) +
    theme_vcfcache() +
    theme(axis.text.x = element_text(angle = 22, hjust = 1))

  diagnostic_figure <- (p_model | p_overhead) / p_efficiency +
    plot_layout(heights = c(1, 0.85)) +
    plot_annotation(
      title = "Alternative view: testing the runtime-scaling model",
      subtitle = "Panels, exomes and genomes place the same cache-hit model on very different absolute time scales",
      caption = paste(
        preliminary_caption(snapshot),
        "1000 Genomes-derived points calibrate runtime conditional on hit rate; they do not estimate public-cache coverage.",
        "KPGP, SGDP and PGP provide the real-world WGS distributions and have no documented project overlap with gnomAD.",
        "Negative inferred overhead values can arise from ordinary runtime variation and condition ordering.",
        sep = "\n"
      ),
      tag_levels = "A",
      theme = theme(
        plot.title = element_text(size = 18, face = "bold", color = vcf_colors[["ink"]]),
        plot.subtitle = element_text(size = 11, color = vcf_colors[["grey"]]),
        plot.caption = element_text(size = 8, color = vcf_colors[["grey"]], hjust = 0)
      )
    )

  diagnostic_panels <- file.path(
    output_dir, "alternatives", "panels", "runtime_diagnostics"
  )
  save_plot(p_model, file.path(diagnostic_panels, "A_miss_fraction_model"), 7.4, 5.0)
  save_plot(p_overhead, file.path(diagnostic_panels, "B_inferred_overhead"), 7.4, 5.0)
  save_plot(p_efficiency, file.path(diagnostic_panels, "C_speedup_efficiency"), 14.8, 4.7)
  save_plot(
    diagnostic_figure,
    file.path(output_dir, "alternatives", "runtime_model_diagnostics"),
    15,
    9.8
  )

  invisible(list(decision = decision_figure, diagnostics = diagnostic_figure))
}

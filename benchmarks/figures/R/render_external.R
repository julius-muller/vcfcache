render_external_figure <- function(input_dir, output_dir, snapshot) {
  rows <- read_tsv(file.path(input_dir, "external_wgs_preliminary_metrics.tsv"))
  strategies_manifest <- jsonlite::fromJSON(
    file.path(input_dir, "external_strategies.json"),
    simplifyVector = FALSE
  )
  if (nrow(rows) == 0) {
    stop("No completed external WGS rows are available")
  }

  strategy_levels <- c("gnomad_af_0.1", "gnomad_af_0.01", "cohort_3_genomes")
  strategy_labels <- c(
    gnomad_af_0.1 = "gnomAD\n≥ 10%",
    gnomad_af_0.01 = "gnomAD\n≥ 1%",
    cohort_3_genomes = "Cohort\n3 genomes"
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

  rows$strategy <- factor(rows$strategy, levels = strategy_levels)
  rows$strategy_label <- factor(
    strategy_labels[as.character(rows$strategy)],
    levels = unname(strategy_labels[strategy_levels])
  )
  rows$cohort_label <- factor(
    cohort_labels[rows$cohort],
    levels = c(cohort_labels[["kpgp"]], cohort_labels[["sgdp"]], cohort_labels[["pgp"]])
  )
  rows$validation_class <- ifelse(
    rows$validation_status == "semantically_validated",
    "Formally validated",
    "Comparator post-processing pending"
  )

  keys <- unique(rows[c("cohort", "assembly", "strategy")])
  summaries <- do.call(
    rbind,
    lapply(seq_len(nrow(keys)), function(index) {
      key <- keys[index, ]
      group <- rows[
        rows$cohort == key$cohort &
          rows$assembly == key$assembly &
          rows$strategy == key$strategy,
      ]
      data.frame(
        cohort = key$cohort,
        assembly = key$assembly,
        strategy = as.character(key$strategy),
        samples = length(unique(group$sample)),
        formally_validated_samples = length(unique(group$sample[group$validation_status == "semantically_validated"])),
        median_hit_rate = median(group$cache_hit_rate),
        hit_rate_q1 = unname(quantile(group$cache_hit_rate, 0.25)),
        hit_rate_q3 = unname(quantile(group$cache_hit_rate, 0.75)),
        median_relative_runtime = median(group$relative_runtime),
        relative_runtime_q1 = unname(quantile(group$relative_runtime, 0.25)),
        relative_runtime_q3 = unname(quantile(group$relative_runtime, 0.75)),
        median_speedup = median(group$speedup),
        median_uncached_seconds = median(group$uncached_wall_seconds),
        median_cached_seconds = median(group$cached_wall_seconds),
        stringsAsFactors = FALSE
      )
    })
  )
  summaries <- summaries[order(summaries$cohort, match(summaries$strategy, strategy_levels)), ]
  write_tsv(summaries, file.path(output_dir, "source", "external_wgs_summary_preliminary.tsv"))

  p_strategy <- ggplot(rows, aes(x = strategy_label, y = relative_runtime, color = strategy)) +
    geom_boxplot(width = 0.60, outlier.shape = NA, color = vcf_colors[["ink"]], fill = "white") +
    geom_point(
      position = position_jitter(width = 0.15, height = 0, seed = 20260805),
      aes(alpha = validation_class), size = 1.7
    ) +
    facet_wrap(~cohort_label, nrow = 1) +
    scale_color_manual(values = strategy_colors, guide = "none") +
    scale_alpha_manual(
      values = c("Formally validated" = 0.85, "Comparator post-processing pending" = 0.38),
      name = NULL
    ) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, NA)) +
    labs(
      title = "Runtime remaining by cache strategy",
      subtitle = "Evaluation genomes were excluded from cohort-cache construction",
      x = NULL,
      y = "Cached / uncached wall time"
    ) +
    theme_vcfcache() +
    theme(axis.text.x = element_text(size = 8), legend.position = "bottom")

  ideal <- data.frame(
    cache_hit_rate = seq(0.45, 0.99, length.out = 300)
  )
  ideal$relative_runtime <- 1 - ideal$cache_hit_rate
  p_model <- ggplot(
    rows,
    aes(
      x = cache_hit_rate,
      y = relative_runtime,
      fill = strategy,
      shape = cohort
    )
  ) +
    geom_line(
      data = ideal,
      aes(x = cache_hit_rate, y = relative_runtime),
      inherit.aes = FALSE,
      color = vcf_colors[["grey"]], linewidth = 0.8, linetype = 2
    ) +
    geom_point(color = vcf_colors[["ink"]], stroke = 0.45, size = 2.5, alpha = 0.72) +
    scale_fill_manual(values = strategy_colors, labels = gsub("\n", " ", strategy_labels), name = NULL) +
    scale_shape_manual(values = cohort_shapes, labels = toupper(names(cohort_shapes)), name = NULL) +
    scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    labs(
      title = "Observed behavior follows the miss fraction",
      subtitle = "Dashed line: zero-overhead rule; vertical excess is lookup and preprocessing",
      x = "Observed cache hit rate",
      y = "Cached / uncached wall time"
    ) +
    guides(
      fill = guide_legend(
        order = 1,
        nrow = 1,
        override.aes = list(shape = 21, size = 3.5, color = vcf_colors[["ink"]], alpha = 1)
      ),
      shape = guide_legend(
        order = 2,
        nrow = 1,
        override.aes = list(size = 3.5, fill = vcf_colors[["light_grey"]], color = vcf_colors[["ink"]], alpha = 1)
      )
    ) +
    theme_vcfcache() +
    theme(
      legend.text = element_text(size = 8),
      legend.key.width = grid::unit(0.8, "lines"),
      legend.spacing.x = grid::unit(0.3, "lines")
    )

  build_costs <- vapply(
    strategies_manifest$cohort_strategies,
    function(value) as.numeric(value$build_wall_seconds),
    numeric(1)
  )
  custom <- summaries[summaries$strategy == "cohort_3_genomes", ]
  cohort_sizes <- unique(round(exp(seq(log(1), log(1000), length.out = 180))))
  economics <- do.call(
    rbind,
    lapply(seq_len(nrow(custom)), function(index) {
      row <- custom[index, ]
      build <- build_costs[[row$cohort]]
      effective <- row$median_cached_seconds + build / cohort_sizes
      data.frame(
        cohort = row$cohort,
        cohort_label = cohort_labels[[row$cohort]],
        cohort_size = cohort_sizes,
        cache_build_seconds = build,
        effective_seconds_per_sample = effective,
        effective_speedup = row$median_uncached_seconds / effective,
        stringsAsFactors = FALSE
      )
    })
  )
  break_even <- do.call(
    rbind,
    lapply(seq_len(nrow(custom)), function(index) {
      row <- custom[index, ]
      saving <- row$median_uncached_seconds - row$median_cached_seconds
      data.frame(
        cohort = row$cohort,
        assembly = row$assembly,
        cache_build_seconds = build_costs[[row$cohort]],
        median_seconds_saved_per_sample = saving,
        break_even_samples = ifelse(saving > 0, ceiling(build_costs[[row$cohort]] / saving), Inf),
        stringsAsFactors = FALSE
      )
    })
  )
  write_tsv(economics, file.path(output_dir, "source", "external_cohort_economics_preliminary.tsv"))
  write_tsv(break_even, file.path(output_dir, "source", "external_cohort_break_even_preliminary.tsv"))

  economics$cohort_label <- factor(economics$cohort_label, levels = cohort_labels[c("kpgp", "sgdp", "pgp")])
  p_economics <- ggplot(economics, aes(x = cohort_size, y = effective_speedup, color = cohort)) +
    geom_hline(yintercept = 1, color = vcf_colors[["red"]], linewidth = 0.6, linetype = 2) +
    geom_line(linewidth = 1.05) +
    facet_wrap(~cohort_label, nrow = 1) +
    scale_color_manual(values = c(kpgp = vcf_colors[["blue"]], sgdp = vcf_colors[["green"]], pgp = vcf_colors[["red"]]), guide = "none") +
    scale_x_log10(breaks = c(1, 2, 5, 10, 100, 1000), labels = scales::label_number()) +
    scale_y_continuous(labels = scales::label_number(suffix = "×", accuracy = 0.1)) +
    labs(
      title = "A custom cache repays its one-time build cost",
      subtitle = "Effective speedup includes measured three-genome cache construction",
      x = "Samples sharing the custom cache (log scale)",
      y = "Effective speedup per sample"
    ) +
    theme_vcfcache()

  combined <- (p_strategy | p_model) / p_economics +
    plot_layout(heights = c(1.05, 0.85)) +
    plot_annotation(
      title = "Independent real-world WGS: bundled versus cohort-derived caches",
      subtitle = "KPGP, SGDP and PGP evaluation genomes · excluded from cohort-cache construction · one run per sample",
      caption = paste(
        preliminary_caption(snapshot),
        "Greyed points are runtime-complete but await corrected semantic comparator post-processing; do not cite this draft.",
        "Assembly-specific cohorts remain separated. Bundled caches are the vcfcache Zenodo releases.",
        sep = "\n"
      ),
      tag_levels = "A",
      theme = theme(
        plot.title = element_text(size = 18, face = "bold", color = vcf_colors[["ink"]]),
        plot.subtitle = element_text(size = 11, color = vcf_colors[["grey"]]),
        plot.caption = element_text(size = 8, color = vcf_colors[["red"]], hjust = 0),
        plot.tag = element_text(face = "bold", color = vcf_colors[["blue"]])
      )
    )

  panel_dir <- file.path(output_dir, "manuscript", "panels", "external_wgs")
  save_plot(p_strategy, file.path(panel_dir, "A_cache_strategy"), 7.4, 5.0)
  save_plot(p_model, file.path(panel_dir, "B_hit_rate_model"), 7.4, 5.0)
  save_plot(p_economics, file.path(panel_dir, "C_build_amortization"), 14.8, 4.7)
  save_plot(combined, file.path(output_dir, "manuscript", "external_wgs_preliminary"), 15, 10)
  invisible(list(plot = combined, summary = summaries, economics = economics))
}

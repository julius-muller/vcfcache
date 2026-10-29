# Shared helpers for the final, self-contained R/ggplot2 analysis.

#### INPUT AND PRESENTATION ####

read_table <- function(path) {
  tibble::as_tibble(read.delim(path, check.names = FALSE, stringsAsFactors = FALSE,
                             na.strings = c("", "NA")))
}
write_table <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}
truth <- function(x) tolower(as.character(x)) %in% c("true", "1")
recipe_names <- c("VEP everything", "fastVEP basic", "fastVEP clinical",
                  "SnpEff core", "SnpEff clinical", "vep-rs core", "vep-rs clinical",
                  "VEP simulated delay")
cohort_colours <- c(KPGP = "#306B8D", SGDP = "#AB6036", PGP = "#706095")
strategy_colours <- c("gnomAD AF ≥ 10%" = "#3E718D", "gnomAD AF ≥ 1%" = "#2D826B",
                      "Three-genome cohort cache" = "#B57537")

enrich <- function(d) {
  d$validated <- truth(d$validated)
  d$is_alias <- truth(d$is_alias)
  d$direct_reused_from_failed_attempt <- truth(d$direct_reused_from_failed_attempt)
  d$carried_forward <- truth(d$carried_forward)
  d$recipe <- ifelse(d$tool == "vep", "VEP everything",
    ifelse(d$tool == "fastvep", ifelse(d$level == "clinical", "fastVEP clinical", "fastVEP basic"),
      ifelse(d$tool == "snpeff", paste("SnpEff", d$level), paste("vep-rs", d$level))))
  d$recipe[d$family == "synthetic_pipeline_spectrum"] <- "VEP simulated delay"
  d$recipe <- factor(d$recipe, levels = recipe_names)
  d$assay <- factor(toupper(d$assay), levels = c("PANEL", "WES", "WGS"),
                    labels = c("Panel", "WES", "WGS"))
  d$cohort <- factor(toupper(d$cohort), levels = c("KPGP", "SGDP", "PGP"))
  d$strategy_label <- ifelse(grepl("0.01$|1%", d$strategy) & !grepl("10%", d$strategy),
      "gnomAD AF ≥ 1%", ifelse(grepl("0.1$|10%", d$strategy),
      "gnomAD AF ≥ 10%", "Three-genome cohort cache"))
  d$strategy_label <- factor(d$strategy_label, levels = names(strategy_colours))
  d$cache_to_direct <- d$cached_seconds / d$direct_seconds
  d$saved_minutes <- d$saved_seconds / 60
  d$saved_hours <- d$saved_seconds / 3600
  # This is a model residual, NOT measured lookup overhead. Misses need not
  # have the mean annotation cost, and annotator startup can remain per call.
  d$equal_cost_model_residual_seconds <- d$cached_seconds -
    (1 - d$cache_hit_rate) * d$direct_seconds
  d$required_cache_side_reduction_for_parity <- pmax(0,
    (d$cached_seconds - d$direct_seconds) / d$cache_side_residual_seconds)
  d$cached_peak_gib <- d$cached_peak_bytes / 1024^3
  d$direct_peak_gib <- d$direct_peak_bytes / 1024^3
  d$delay_us <- NA_real_
  synthetic <- d$family == "synthetic_pipeline_spectrum"
  d$delay_us[synthetic] <- as.numeric(sub(".*delay_([0-9]+)us.*", "\\1", d$level[synthetic]))
  d
}

journal_theme <- function(size = 9) {
  ggplot2::theme_classic(base_size = size, base_family = "sans") +
    ggplot2::theme(axis.line = ggplot2::element_line(linewidth = .35, colour = "#313B42"),
      axis.ticks = ggplot2::element_line(linewidth = .3),
      axis.text = ggplot2::element_text(colour = "#27343C"),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold", size = size),
      legend.position = "bottom", legend.title = ggplot2::element_text(size = size),
      legend.key.width = grid::unit(12, "pt"),
      plot.title = ggplot2::element_blank(), plot.subtitle = ggplot2::element_blank(),
      plot.caption = ggplot2::element_blank(),
      plot.tag = ggplot2::element_text(face = "bold", size = size + 3),
      plot.margin = ggplot2::margin(7, 10, 5, 6))
}

save_figure <- function(plot, prefix, width, height) {
  dir.create(dirname(prefix), recursive = TRUE, showWarnings = FALSE)
  # Re-rendering must never silently overwrite a previously printed figure.
  if (any(file.exists(paste0(prefix, c(".pdf", ".svg", ".png"))))) {
    stop("Figure exists; select a new figure directory: ", prefix)
  }
  grDevices::cairo_pdf(paste0(prefix, ".pdf"), width = width, height = height,
                      family = "sans", bg = "white")
  print(plot); grDevices::dev.off()
  grDevices::svg(paste0(prefix, ".svg"), width = width, height = height,
                 family = "sans", bg = "white")
  print(plot); grDevices::dev.off()
  grDevices::png(paste0(prefix, ".png"), width = width, height = height,
                 units = "in", res = 600, type = "cairo", bg = "white")
  print(plot); grDevices::dev.off()
  invisible(prefix)
}

#### SAMPLE-LEVEL INFERENCE ####

bootstrap_interval <- function(values, cohort, replicates = 4000L) {
  if (length(values) < 3L) return(c(NA_real_, NA_real_))
  # Maintain the deliberately heterogeneous cohort allocation. One independent
  # sample per condition is the resampling unit, not technical repeats.
  strata <- split(seq_along(values), cohort, drop = TRUE)
  medians <- replicate(replicates, {
    index <- unlist(lapply(strata, function(i) i[sample.int(length(i), length(i), replace = TRUE)]),
                    use.names = FALSE)
    median(values[index])
  })
  as.numeric(quantile(medians, c(.025, .975), names = FALSE))
}

split_groups <- function(d, by) {
  keys <- do.call(interaction, c(unname(as.list(d[by])), list(drop = TRUE, lex.order = TRUE)))
  split(seq_len(nrow(d)), keys, drop = TRUE)
}

summarise_conditions <- function(d, by) {
  set.seed(20261029)
  groups <- split_groups(d, by)
  result <- lapply(groups, function(index) {
    x <- d[index, , drop = FALSE]
    if (anyDuplicated(x$sample)) stop("Repeated sample within summary: ", paste(by, collapse = "/"))
    ci <- bootstrap_interval(x$speedup, x$cohort)
    cbind(as.data.frame(x[1, by, drop = FALSE]), data.frame(
      n = nrow(x), median_speedup = median(x$speedup), speedup_ci_low = ci[1],
      speedup_ci_high = ci[2], speedup_p10 = unname(quantile(x$speedup, .1)),
      speedup_p90 = unname(quantile(x$speedup, .9)),
      min_speedup = min(x$speedup), max_speedup = max(x$speedup),
      n_faster = sum(x$speedup > 1), fraction_faster = mean(x$speedup > 1),
      median_direct_seconds = median(x$direct_seconds), median_cached_seconds = median(x$cached_seconds),
      median_saved_seconds = median(x$saved_seconds), median_saved_fraction = median(x$saved_fraction),
      median_hit_rate = median(x$cache_hit_rate), min_hit_rate = min(x$cache_hit_rate),
      max_hit_rate = max(x$cache_hit_rate),
      median_annotation_fraction = median(x$annotation_fraction),
      median_cache_side_seconds = median(x$cache_side_residual_seconds),
      median_zero_cache_side_gain = median(x$hypothetical_zero_non_annotation_gain, na.rm = TRUE),
      median_required_side_reduction_for_parity = median(x$required_cache_side_reduction_for_parity),
      total_direct_seconds = sum(x$direct_seconds), total_cached_seconds = sum(x$cached_seconds),
      total_saved_seconds = sum(x$saved_seconds),
      aggregate_throughput_gain = sum(x$direct_seconds) / sum(x$cached_seconds),
      median_historical_to_new_cached_gain = median(x$historical_to_new_speedup)
    ))
  })
  tibble::as_tibble(do.call(rbind, result))
}

phase_data <- function(d) {
  mapping <- c(annotate_misses_seconds = "Annotate misses", route_seconds = "Lookup / routing",
    preprocess_seconds = "Preprocessing", assemble_seconds = "Assembly", engine_other_seconds = "Other engine",
    cli_other_seconds = "CLI / other")
  result <- do.call(rbind, lapply(names(mapping), function(name) {
    x <- d
    x$phase <- mapping[[name]]
    x$phase_seconds <- x[[name]]
    x$fraction <- x$phase_seconds / x$cached_seconds
    x
  }))
  result$phase <- factor(result$phase, levels = unname(mapping))
  tibble::as_tibble(result)
}

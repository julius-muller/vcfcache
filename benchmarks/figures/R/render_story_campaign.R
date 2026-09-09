#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5) {
  stop(
    "Usage: render_story_campaign.R <native_source_directory> ",
    "<pipeline_complexity.tsv> <cohort_economics.tsv> ",
    "<cohort_break_even.tsv> <output_directory>"
  )
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(tibble)
})

full_arguments <- commandArgs(trailingOnly = FALSE)
file_argument <- full_arguments[grepl("^--file=", full_arguments)]
if (length(file_argument) != 1) {
  stop("Could not determine the render_story_campaign.R location")
}
script_dir <- dirname(normalizePath(sub("^--file=", "", file_argument)))
source(file.path(script_dir, "load_benchmark_results.R"))

source_dir <- normalizePath(args[[1]], mustWork = TRUE)
complexity_path <- normalizePath(args[[2]], mustWork = TRUE)
economics_path <- normalizePath(args[[3]], mustWork = TRUE)
break_even_path <- normalizePath(args[[4]], mustWork = TRUE)
output_dir <- args[[5]]

for (subdir in c(
  "source", "manuscript", "manuscript/panels", "supplement",
  "supplement/panels", "repository", "alternatives"
)) {
  dir.create(file.path(output_dir, subdir), recursive = TRUE, showWarnings = FALSE)
}

write_tsv <- function(value, path) {
  write.table(value, path, sep = "\t", row.names = FALSE, quote = FALSE, na = "")
}

save_plot <- function(plot, prefix, width, height, dpi = 600) {
  svglite::svglite(
    paste0(prefix, ".svg"), width = width, height = height,
    pointsize = 10, bg = "white"
  )
  print(plot)
  grDevices::dev.off()

  grDevices::cairo_pdf(
    paste0(prefix, ".pdf"), width = width, height = height,
    pointsize = 10, family = "sans", bg = "white"
  )
  print(plot)
  grDevices::dev.off()

  grDevices::png(
    paste0(prefix, ".png"), width = width, height = height,
    units = "in", res = dpi, type = "cairo", bg = "white"
  )
  print(plot)
  grDevices::dev.off()
}

ink <- "#172832"
blue <- "#2F7195"
orange <- "#D17B32"
green <- "#3C8B62"
purple <- "#7863A8"
red <- "#B55454"
grey <- "#6D7B81"
light_grey <- "#E8EDEF"
grid_grey <- "#DCE3E5"
paper <- "#FCFDFC"
tool_colours <- c("VEP" = blue, "fastVEP" = orange, "SnpEff" = purple)
assay_levels <- c("Panel", "WES", "WGS")
tool_levels <- c("VEP", "fastVEP", "SnpEff")

journal_theme <- function(base_size = 9.5) {
  theme_classic(base_size = base_size, base_family = "sans") +
    theme(
      text = element_text(colour = ink),
      axis.text = element_text(colour = ink),
      axis.title = element_text(colour = ink),
      axis.line = element_line(colour = ink, linewidth = 0.35),
      axis.ticks = element_line(colour = ink, linewidth = 0.35),
      panel.grid.major.y = element_line(colour = grid_grey, linewidth = 0.3),
      panel.grid.minor = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", colour = ink),
      legend.position = "bottom",
      legend.title = element_blank(),
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      plot.caption = element_blank(),
      plot.tag = element_text(face = "bold", size = 12, colour = ink),
      plot.tag.position = c(0.005, 0.995),
      plot.margin = margin(8, 8, 8, 8)
    )
}

format_duration <- function(seconds) {
  vapply(seconds, function(value) {
    if (value < 60) {
      sprintf("%.0f s", value)
    } else if (value < 3600) {
      sprintf("%.0f min", value / 60)
    } else {
      sprintf("%.1f h", value / 3600)
    }
  }, character(1))
}

summarize_groups <- function(data, fields) {
  keys <- unique(data[fields])
  result <- lapply(seq_len(nrow(keys)), function(index) {
    keep <- rep(TRUE, nrow(data))
    for (field in fields) {
      keep <- keep & data[[field]] == keys[[field]][index]
    }
    values <- data[keep, ]
    cbind(
      keys[index, , drop = FALSE],
      data.frame(
        n = nrow(values),
        direct_seconds = median(values$direct_seconds),
        cached_seconds = median(values$cached_seconds),
        speedup = median(values$speedup),
        relative_runtime = median(values$relative_runtime),
        time_saved_seconds = median(values$time_saved_seconds),
        hit_rate = median(values$hit_rate)
      )
    )
  })
  do.call(rbind, result)
}

three_path <- file.path(source_dir, "three_annotator_assay_pairs.tsv")
snpeff_path <- file.path(source_dir, "snpeff_load_pairs.tsv")
campaign_path <- file.path(source_dir, "campaign_metrics.tsv")
for (path in c(three_path, snpeff_path, campaign_path)) {
  if (!file.exists(path)) {
    stop("Missing source table: ", path)
  }
}

three <- load_benchmark_table(three_path)
three$tool_label <- factor(three$tool_label, levels = tool_levels)
three$assay_label <- factor(three$assay_label, levels = assay_levels)
three$cohort_label <- factor(three$cohort_label, levels = c("KPGP", "SGDP"))
three$remaining_percent <- 100 * three$relative_runtime
three$returned_percent <- 100 * (1 - three$relative_runtime)
three$direct_per_day <- 86400 / three$direct_seconds
three$cached_per_day <- 86400 / three$cached_seconds
three$ideal_remaining_percent <- 100 * (1 - three$hit_rate)

if (nrow(three) != 108 || !all(three$remaining_percent > 0)) {
  stop("Expected 108 positive matched native-annotator observations")
}

three_summary <- summarize_groups(three, c("tool_label", "assay_label"))
three_summary$tool_label <- factor(three_summary$tool_label, levels = tool_levels)
three_summary$assay_label <- factor(three_summary$assay_label, levels = assay_levels)
three_summary$remaining_percent <- 100 * three_summary$relative_runtime
three_summary$returned_percent <- 100 * (1 - three_summary$relative_runtime)
three_summary$direct_per_day <- 86400 / three_summary$direct_seconds
three_summary$cached_per_day <- 86400 / three_summary$cached_seconds
three_summary$throughput_gain <- three_summary$cached_per_day / three_summary$direct_per_day

complexity <- load_benchmark_table(complexity_path)
complexity <- complexity[tolower(complexity$semantic_pass) == "true", ]
if (nrow(complexity) != 12) {
  stop("Expected 12 validated pipeline-complexity observations")
}
complexity$cache_label <- factor(
  complexity$cache_strategy,
  levels = c("AF ≥ 10% cache", "AF ≥ 1% cache"),
  labels = c("gnomAD AF ≥ 10%", "gnomAD AF ≥ 1%")
)
complexity$direct_hours <- complexity$uncached_wall_seconds / 3600
complexity$cached_hours <- complexity$cached_wall_seconds / 3600
complexity$returned_hours <- complexity$time_saved_seconds / 3600
complexity$maximum_returnable_hours <-
  complexity$cache_hit_rate * complexity$direct_hours
complexity$return_efficiency <-
  100 * complexity$returned_hours / complexity$maximum_returnable_hours

snpeff <- load_benchmark_table(snpeff_path)
snpeff$recipe_label <- factor(snpeff$recipe_label, levels = c("Core", "Clinical"))
snpeff$assay_label <- factor(snpeff$assay_label, levels = assay_levels)

economics <- load_benchmark_table(economics_path)
economics$cohort_label <- factor(
  economics$cohort,
  levels = c("kpgp", "sgdp", "pgp"),
  labels = c("KPGP · GRCh38", "SGDP · GRCh38", "PGP · GRCh37")
)
economics$direct_seconds <-
  economics$effective_seconds_per_sample * economics$effective_speedup
economics$cached_seconds <-
  economics$effective_seconds_per_sample -
  economics$cache_build_seconds / economics$cohort_size

break_even <- load_benchmark_table(break_even_path)
break_even$cohort_label <- factor(
  break_even$cohort,
  levels = c("kpgp", "sgdp", "pgp"),
  labels = c("KPGP · GRCh38", "SGDP · GRCh38", "PGP · GRCh37")
)

economic_parameters <- do.call(rbind, lapply(split(economics, economics$cohort), function(values) {
  data.frame(
    cohort = values$cohort[[1]],
    cohort_label = values$cohort_label[[1]],
    direct_seconds = median(values$direct_seconds),
    cached_seconds = median(values$cached_seconds),
    build_seconds = median(values$cache_build_seconds)
  )
}))
economic_parameters$cohort_label <- factor(
  economic_parameters$cohort_label,
  levels = levels(economics$cohort_label)
)

economic_sizes <- c(10, 100, 1000)
economic_bars <- do.call(rbind, lapply(seq_len(nrow(economic_parameters)), function(index) {
  values <- economic_parameters[index, ]
  data.frame(
    cohort_label = values$cohort_label,
    cohort_size = economic_sizes,
    days_returned = (
      economic_sizes * values$direct_seconds -
        (values$build_seconds + economic_sizes * values$cached_seconds)
    ) / 86400
  )
}))
economic_bars$cohort_size_label <- factor(
  paste0("n = ", economic_bars$cohort_size),
  levels = paste0("n = ", economic_sizes)
)

campaign <- load_benchmark_table(campaign_path)
pair_fields <- c("sample", "cohort", "tool", "level", "assembly", "assay")
campaign_direct <- campaign[campaign$condition == "direct", ]
campaign_cached <- campaign[campaign$condition == "cached", ]
campaign_pairs <- merge(
  campaign_direct, campaign_cached, by = pair_fields,
  suffixes = c("_direct", "_cached")
)
campaign_pairs$audit_outcome <- ifelse(
  campaign_pairs$concordant_cached == "PASS",
  "Eligible and concordant", "Excluded: annotator order effect"
)
audit_counts <- as.data.frame(table(campaign_pairs$audit_outcome), stringsAsFactors = FALSE)
names(audit_counts) <- c("outcome", "n")
audit_counts$outcome <- factor(
  audit_counts$outcome,
  levels = c("Eligible and concordant", "Excluded: annotator order effect")
)
audit_counts$xmin <- c(0, 176)
audit_counts$xmax <- c(176, 180)

write_tsv(three, file.path(output_dir, "source", "native_assay_pairs.tsv"))
write_tsv(three_summary, file.path(output_dir, "source", "native_assay_summary.tsv"))
write_tsv(complexity, file.path(output_dir, "source", "pipeline_complexity.tsv"))
write_tsv(snpeff, file.path(output_dir, "source", "snpeff_load_pairs.tsv"))
write_tsv(economics, file.path(output_dir, "source", "cohort_economics.tsv"))
write_tsv(economic_bars, file.path(output_dir, "source", "cohort_time_returned.tsv"))
write_tsv(audit_counts, file.path(output_dir, "source", "semantic_audit.tsv"))

# Figure 1A: application-level workflow, drawn as a vector-native ggplot.
workflow <- ggplot() +
  annotate("rect", xmin = 0.15, xmax = 1.55, ymin = 2.2, ymax = 3.3,
    fill = "#E7F0F5", colour = blue, linewidth = 0.7) +
  annotate("text", x = 0.85, y = 2.75, label = "New\nsample VCF", colour = ink,
    fontface = "bold", size = 2.75) +
  annotate("rect", xmin = 2.05, xmax = 3.75, ymin = 2.2, ymax = 3.3,
    fill = "#EFF4F5", colour = grey, linewidth = 0.7) +
  annotate("text", x = 2.9, y = 2.75, label = "Match\nvariants to cache", colour = ink,
    fontface = "bold", size = 2.75) +
  annotate("rect", xmin = 2.0, xmax = 3.8, ymin = 4.1, ymax = 5.0,
    fill = "#EEEAF6", colour = purple, linewidth = 0.65) +
  annotate("text", x = 2.9, y = 4.55,
    label = "Zenodo recipe cache\nor blueprint", colour = ink, size = 2.75) +
  annotate("rect", xmin = 4.35, xmax = 6.15, ymin = 3.15, ymax = 4.15,
    fill = "#E5F2EA", colour = green, linewidth = 0.7) +
  annotate("text", x = 5.25, y = 3.65,
    label = "Reuse cached\nannotations", colour = ink, fontface = "bold", size = 3.0) +
  annotate("rect", xmin = 4.35, xmax = 6.15, ymin = 1.35, ymax = 2.35,
    fill = "#F8EBDF", colour = orange, linewidth = 0.7) +
  annotate("text", x = 5.25, y = 1.85,
    label = "Annotate only\nuncached variants", colour = ink, fontface = "bold", size = 2.85) +
  annotate("rect", xmin = 6.75, xmax = 8.0, ymin = 2.2, ymax = 3.3,
    fill = "#EFF4F5", colour = grey, linewidth = 0.7) +
  annotate("text", x = 7.375, y = 2.75, label = "Merge +\nsort", colour = ink,
    fontface = "bold", size = 3.0) +
  annotate("rect", xmin = 8.55, xmax = 10.5, ymin = 2.2, ymax = 3.3,
    fill = "#E7F0F5", colour = blue, linewidth = 0.7) +
  annotate("text", x = 9.525, y = 2.75, label = "Complete\nannotated VCF", colour = ink,
    fontface = "bold", size = 2.9) +
  annotate("segment", x = 1.55, xend = 2.0, y = 2.75, yend = 2.75,
    colour = ink, linewidth = 0.65,
    arrow = grid::arrow(type = "closed", length = grid::unit(0.07, "in"))) +
  annotate("segment", x = 2.9, xend = 2.9, y = 4.1, yend = 3.35,
    colour = purple, linewidth = 0.65,
    arrow = grid::arrow(type = "closed", length = grid::unit(0.07, "in"))) +
  annotate("segment", x = 3.75, xend = 4.3, y = 2.95, yend = 3.55,
    colour = green, linewidth = 0.75,
    arrow = grid::arrow(type = "closed", length = grid::unit(0.07, "in"))) +
  annotate("segment", x = 3.75, xend = 4.3, y = 2.55, yend = 1.95,
    colour = orange, linewidth = 0.75,
    arrow = grid::arrow(type = "closed", length = grid::unit(0.07, "in"))) +
  annotate("text", x = 4.0, y = 3.43, label = "hit", colour = green,
    fontface = "bold", size = 2.8) +
  annotate("text", x = 4.0, y = 2.07, label = "miss", colour = orange,
    fontface = "bold", size = 2.8) +
  annotate("segment", x = 6.15, xend = 6.7, y = 3.65, yend = 2.95,
    colour = green, linewidth = 0.75,
    arrow = grid::arrow(type = "closed", length = grid::unit(0.07, "in"))) +
  annotate("segment", x = 6.15, xend = 6.7, y = 1.85, yend = 2.55,
    colour = orange, linewidth = 0.75,
    arrow = grid::arrow(type = "closed", length = grid::unit(0.07, "in"))) +
  annotate("segment", x = 8.0, xend = 8.5, y = 2.75, yend = 2.75,
    colour = ink, linewidth = 0.65,
    arrow = grid::arrow(type = "closed", length = grid::unit(0.07, "in"))) +
  coord_cartesian(xlim = c(0, 10.7), ylim = c(0.8, 5.3), clip = "off") +
  labs(tag = "A") +
  theme_void(base_family = "sans") +
  theme(plot.tag = element_text(face = "bold", size = 12, colour = ink),
    plot.tag.position = c(0.005, 0.995), plot.margin = margin(12, 12, 8, 12))

drivers <- data.frame(
  x = c(1, 3, 5),
  label = c("Cache coverage", "Annotation cost", "Cohort size"),
  detail = c("More hits\nless new work", "Slower pipelines\nmore time returned", "One-time build cost\nis amortized"),
  colour = c(green, blue, purple)
)
driver_plot <- ggplot(drivers, aes(x = x, y = 1)) +
  geom_tile(aes(fill = colour), width = 1.55, height = 1.15,
    colour = "white", linewidth = 1.2, alpha = 0.14) +
  geom_point(aes(colour = colour), y = 1.3, size = 5.8) +
  geom_text(aes(label = label), y = 1.05, fontface = "bold", size = 3.5,
    colour = ink) +
  geom_text(aes(label = detail), y = 0.72, size = 3.0, lineheight = 0.95,
    colour = ink) +
  scale_fill_identity() +
  scale_colour_identity() +
  coord_cartesian(xlim = c(0.1, 5.9), ylim = c(0.4, 1.55), clip = "off") +
  labs(tag = "B") +
  theme_void(base_family = "sans") +
  theme(plot.tag = element_text(face = "bold", size = 12, colour = ink),
    plot.tag.position = c(0.005, 0.995), plot.margin = margin(4, 12, 12, 12))

figure1 <- workflow / driver_plot + plot_layout(heights = c(1.85, 1))
save_plot(figure1, file.path(output_dir, "manuscript", "figure_workflow_and_value"), 7.2, 5.4)
save_plot(workflow, file.path(output_dir, "manuscript", "panels", "figure_workflow_A"), 7.2, 3.4)
save_plot(driver_plot, file.path(output_dir, "manuscript", "panels", "figure_workflow_B"), 7.2, 2.0)

# Figure 2: fraction of time remaining and serial throughput.
p2a <- ggplot(three, aes(x = tool_label, y = remaining_percent, fill = tool_label)) +
  geom_hline(yintercept = 100, colour = grey, linetype = "dashed", linewidth = 0.55) +
  geom_violin(width = 0.72, trim = FALSE, alpha = 0.18, colour = NA) +
  geom_boxplot(width = 0.24, outlier.shape = NA, fill = "white", linewidth = 0.45) +
  geom_point(aes(shape = cohort_label), position = position_jitter(width = 0.1),
    size = 1.5, alpha = 0.72, colour = ink) +
  facet_grid(. ~ assay_label) +
  scale_fill_manual(values = tool_colours, guide = "none") +
  scale_shape_manual(values = c("KPGP" = 16, "SGDP" = 17)) +
  scale_y_log10(
    breaks = c(5, 10, 25, 50, 100, 200),
    labels = function(x) paste0(x, "%")
  ) +
  coord_cartesian(ylim = c(3, 350)) +
  labs(tag = "A", x = NULL, y = "Runtime remaining") +
  journal_theme(9.4) +
  theme(axis.text.x = element_text(angle = 24, hjust = 1))

throughput_long <- rbind(
  data.frame(three_summary, condition = "Direct", genomes_per_day = three_summary$direct_per_day),
  data.frame(three_summary, condition = "With VCFcache", genomes_per_day = three_summary$cached_per_day)
)
throughput_long$condition <- factor(throughput_long$condition,
  levels = c("Direct", "With VCFcache"))

p2b <- ggplot(three_summary, aes(y = tool_label, colour = tool_label)) +
  geom_segment(aes(x = direct_per_day, xend = cached_per_day, yend = tool_label),
    linewidth = 0.85,
    arrow = grid::arrow(type = "closed", length = grid::unit(0.07, "in"))) +
  geom_point(aes(x = direct_per_day), shape = 21, fill = "white", size = 2.7, stroke = 0.75) +
  geom_point(aes(x = cached_per_day), shape = 17, size = 2.9) +
  geom_text(aes(x = pmin(pmax(direct_per_day, cached_per_day) * 1.45, 7500),
    label = sprintf("%.1f×", throughput_gain)),
    hjust = 0, colour = ink, size = 2.5, show.legend = FALSE) +
  facet_grid(. ~ assay_label) +
  scale_colour_manual(values = tool_colours, guide = "none") +
  scale_x_log10(breaks = c(10, 100, 1000), labels = comma,
    limits = c(7, 12000)) +
  labs(tag = "B", x = "Serial genomes per 24 h", y = NULL) +
  journal_theme(9.4) +
  theme(panel.grid.major.y = element_blank())

figure2 <- p2a / p2b + plot_layout(heights = c(1.2, 1))
save_plot(figure2, file.path(output_dir, "manuscript", "figure_empirical_time_and_throughput"), 7.2, 6.3)
save_plot(p2a, file.path(output_dir, "manuscript", "panels", "figure_empirical_A_remaining"), 7.2, 3.4)
save_plot(p2b, file.path(output_dir, "manuscript", "panels", "figure_empirical_B_throughput"), 7.2, 3.0)

# Figure 3: empirical runtime rule and controlled pipeline-cost response.
ideal_lines <- do.call(rbind, lapply(assay_levels, function(assay) {
  data.frame(assay_label = assay, hit_rate = seq(0.68, 0.985, length.out = 120))
}))
ideal_lines$ideal_remaining_percent <- 100 * (1 - ideal_lines$hit_rate)
ideal_lines$assay_label <- factor(ideal_lines$assay_label, levels = assay_levels)

p3a <- ggplot(three, aes(x = 100 * hit_rate, y = remaining_percent, colour = tool_label)) +
  geom_line(data = ideal_lines,
    aes(x = 100 * hit_rate, y = ideal_remaining_percent, group = assay_label),
    inherit.aes = FALSE, colour = grey, linetype = "dashed", linewidth = 0.55) +
  geom_point(aes(size = direct_seconds), alpha = 0.76) +
  facet_grid(. ~ assay_label) +
  scale_colour_manual(values = tool_colours) +
  scale_size_continuous(range = c(1.2, 5.0), trans = "log10",
    breaks = c(10, 600, 7200), labels = format_duration,
    name = "Direct runtime") +
  scale_y_log10(breaks = c(5, 10, 25, 50, 100, 200),
    labels = function(x) paste0(x, "%")) +
  scale_x_continuous(breaks = c(70, 80, 90, 95),
    labels = function(x) paste0(x, "%")) +
  coord_cartesian(xlim = c(68, 99), ylim = c(2, 350)) +
  labs(tag = "A", x = "Cache hit rate", y = "Runtime remaining") +
  journal_theme(9.2) +
  guides(colour = guide_legend(order = 1), size = guide_legend(order = 2))

cache_colours <- c("gnomAD AF ≥ 10%" = orange, "gnomAD AF ≥ 1%" = green)
ideal_complexity <- do.call(rbind, lapply(levels(complexity$cache_label), function(label) {
  hit_rate <- median(complexity$cache_hit_rate[complexity$cache_label == label])
  data.frame(cache_label = label, direct_hours = seq(3.5, 24, length.out = 100), hit_rate = hit_rate)
}))
ideal_complexity$cache_label <- factor(ideal_complexity$cache_label,
  levels = levels(complexity$cache_label))
ideal_complexity$returned_hours <- ideal_complexity$hit_rate * ideal_complexity$direct_hours

p3b <- ggplot(complexity,
  aes(x = direct_hours, y = returned_hours, colour = cache_label, group = cache_label)) +
  geom_line(data = ideal_complexity, linetype = "dashed", linewidth = 0.6, alpha = 0.75) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.7, alpha = 0.9) +
  scale_colour_manual(values = cache_colours) +
  scale_x_continuous(breaks = c(5, 10, 15, 20), limits = c(3.5, 24)) +
  scale_y_continuous(breaks = c(5, 10, 15, 20), limits = c(2.5, 22)) +
  labs(tag = "B", x = "Direct pipeline runtime per genome (hours)",
    y = "Hours returned per genome") +
  journal_theme(9.2)

figure3 <- p3a / p3b + plot_layout(heights = c(1.12, 1))
save_plot(figure3, file.path(output_dir, "manuscript", "figure_runtime_rule"), 7.2, 6.6)
save_plot(p3a, file.path(output_dir, "manuscript", "panels", "figure_runtime_rule_A_empirical"), 7.2, 3.6)
save_plot(p3b, file.path(output_dir, "manuscript", "panels", "figure_runtime_rule_B_complexity"), 6.1, 3.3)

# Figure 4: cache-build amortization and cumulative time returned.
cohort_colours <- c(
  "KPGP · GRCh38" = blue,
  "SGDP · GRCh38" = green,
  "PGP · GRCh37" = purple
)
p4a <- ggplot(economics,
  aes(x = cohort_size, y = effective_speedup, colour = cohort_label)) +
  geom_hline(yintercept = 1, colour = grey, linetype = "dashed", linewidth = 0.55) +
  geom_line(linewidth = 0.85) +
  geom_point(data = break_even,
    aes(x = break_even_samples, y = 1, colour = cohort_label),
    inherit.aes = FALSE, shape = 21, fill = "white", size = 2.8, stroke = 0.8) +
  annotate("text", x = 2.15, y = 1.17, label = "all cohorts: break-even at 2",
    hjust = 0, size = 2.5, colour = ink) +
  scale_colour_manual(values = cohort_colours) +
  scale_x_log10(breaks = c(1, 2, 5, 10, 100, 1000), labels = comma) +
  scale_y_continuous(breaks = c(1, 2, 3, 4, 5), limits = c(0.45, 5.5),
    labels = function(x) paste0(x, "×")) +
  labs(tag = "A", x = "Samples processed after a three-genome cache build",
    y = "Effective speedup including build") +
  journal_theme(9.4) +
  theme(plot.tag.position = c(-0.035, 0.995),
    plot.margin = margin(16, 8, 8, 18))

p4b <- ggplot(economic_bars,
  aes(x = cohort_size_label, y = days_returned, fill = cohort_label)) +
  geom_col(position = position_dodge(width = 0.78), width = 0.68,
    colour = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.1f d", days_returned)),
    position = position_dodge(width = 0.78), vjust = -0.45, size = 2.45,
    colour = ink) +
  scale_fill_manual(values = cohort_colours) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(tag = "B", x = "Cohort size", y = "Cumulative compute days returned") +
  journal_theme(9.4) +
  guides(fill = "none") +
  theme(panel.grid.major.x = element_blank(),
    plot.tag.position = c(-0.035, 0.995), plot.margin = margin(16, 8, 8, 18))

figure4 <- p4a / p4b + plot_layout(heights = c(1.05, 1))
save_plot(figure4, file.path(output_dir, "manuscript", "figure_cohort_economics"), 7.2, 6.3)
save_plot(p4a, file.path(output_dir, "manuscript", "panels", "figure_economics_A_amortization"), 7.2, 3.2)
save_plot(p4b, file.path(output_dir, "manuscript", "panels", "figure_economics_B_days_returned"), 7.2, 3.2)

# Supplement: annotation load and semantic audit.
s1a <- ggplot(snpeff,
  aes(x = recipe_label, y = speedup, group = interaction(sample, assay_label))) +
  geom_hline(yintercept = 1, colour = grey, linetype = "dashed", linewidth = 0.55) +
  geom_line(colour = "#BDC6C9", linewidth = 0.42, alpha = 0.78) +
  geom_point(aes(colour = recipe_label), size = 2.0, alpha = 0.86) +
  facet_grid(. ~ assay_label) +
  scale_colour_manual(values = c("Core" = grey, "Clinical" = purple), guide = "none") +
  scale_y_log10(breaks = c(0.5, 1, 2, 5, 10, 20),
    labels = function(x) paste0(x, "×"), limits = c(0.7, 20)) +
  labs(tag = "A", x = NULL, y = "Speedup over direct SnpEff") +
  journal_theme(9.4)

s1b <- ggplot(audit_counts) +
  geom_rect(aes(xmin = xmin, xmax = xmax, ymin = 0.25, ymax = 0.85, fill = outcome),
    colour = "white", linewidth = 1.2) +
  annotate("text", x = 88, y = 0.55, label = "176 / 176 eligible pairs\nsemantically concordant",
    colour = "white", fontface = "bold", size = 4.1, lineheight = 0.95) +
  annotate("text", x = 158, y = 1.03,
    label = "4 excluded before inference:\nSnpEff input-order effect",
    colour = ink, hjust = 0.5, size = 2.9, lineheight = 0.95) +
  annotate("segment", x = 158, xend = 178, y = 0.98, yend = 0.87,
    colour = grey, linewidth = 0.45) +
  annotate("text", x = 90, y = 0.05, label = "0 cache-attributable annotation mismatches",
    colour = ink, fontface = "bold", size = 3.5) +
  scale_fill_manual(values = c(
    "Eligible and concordant" = green,
    "Excluded: annotator order effect" = "#D6DDE0"
  ), guide = "none") +
  coord_cartesian(xlim = c(0, 180), ylim = c(-0.05, 1.2), clip = "off") +
  labs(tag = "B") +
  theme_void(base_family = "sans") +
  theme(plot.tag = element_text(face = "bold", size = 12, colour = ink),
    plot.tag.position = c(0.005, 0.995), plot.margin = margin(16, 16, 12, 16))

supplement1 <- s1a / s1b + plot_layout(heights = c(1.5, 0.85))
save_plot(supplement1, file.path(output_dir, "supplement", "figure_annotation_load_and_fidelity"), 7.2, 6.0)
save_plot(s1a, file.path(output_dir, "supplement", "panels", "figure_load_A_snpeff"), 7.2, 3.4)
save_plot(s1b, file.path(output_dir, "supplement", "panels", "figure_fidelity_B_audit"), 7.2, 2.3)

# Repository headliner 1: one-day serial throughput for WGS.
wgs <- three_summary[three_summary$assay_label == "WGS", ]
wgs$configuration <- factor(
  c(VEP = "VEP\n--everything", fastVEP = "fastVEP\nclinical", SnpEff = "SnpEff\nclinical")[
    as.character(wgs$tool_label)
  ],
  levels = c("VEP\n--everything", "fastVEP\nclinical", "SnpEff\nclinical")
)
wgs_throughput <- rbind(
  data.frame(wgs, condition = "Direct", genomes_per_day = wgs$direct_per_day),
  data.frame(wgs, condition = "With VCFcache", genomes_per_day = wgs$cached_per_day)
)
wgs_throughput$condition <- factor(wgs_throughput$condition,
  levels = c("Direct", "With VCFcache"))

one_day <- ggplot(wgs_throughput,
  aes(x = condition, y = genomes_per_day, fill = tool_label)) +
  geom_col(width = 0.62, colour = ink, linewidth = 0.45) +
  geom_text(aes(label = sprintf("%.0f", genomes_per_day)), vjust = -0.45,
    colour = ink, fontface = "bold", size = 4.0) +
  facet_grid(. ~ configuration, scales = "free_y") +
  scale_fill_manual(values = tool_colours, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.16))) +
  labs(x = NULL, y = "Serial WGS annotations per 24 h") +
  theme_minimal(base_size = 12, base_family = "sans") +
  theme(
    text = element_text(colour = ink), panel.grid = element_blank(),
    axis.text.y = element_blank(), axis.ticks = element_blank(),
    axis.text.x = element_text(face = "bold", colour = ink),
    axis.title.y = element_text(face = "bold"),
    strip.text = element_text(face = "bold", size = 13, colour = ink),
    plot.margin = margin(16, 16, 10, 16)
  )
save_plot(one_day, file.path(output_dir, "repository", "headliner_one_day_of_compute"), 9.5, 4.8)

# Repository headliner 2: an intuitive direct-to-cached WGS clock.
wgs$direct_minutes <- wgs$direct_seconds / 60
wgs$cached_minutes <- wgs$cached_seconds / 60
wgs$returned_label <- paste0(format_duration(wgs$time_saved_seconds), " returned")
clock <- ggplot(wgs, aes(y = tool_label, colour = tool_label)) +
  geom_segment(aes(x = direct_minutes, xend = cached_minutes, yend = tool_label),
    linewidth = 8, alpha = 0.16, lineend = "round") +
  geom_point(aes(x = direct_minutes), shape = 21, fill = "white", size = 6, stroke = 1.3) +
  geom_point(aes(x = cached_minutes), size = 5.5) +
  geom_label(aes(x = sqrt(direct_minutes * cached_minutes),
    label = sprintf("%.1f× faster\n%s", speedup, returned_label)),
    fill = "white", colour = ink, linewidth = 0.5, size = 3.4,
    lineheight = 0.95, show.legend = FALSE) +
  scale_colour_manual(values = tool_colours, guide = "none") +
  scale_x_log10(breaks = c(2, 5, 15, 30, 60, 180),
    labels = function(x) paste0(x, " min"), limits = c(1.5, 220)) +
  labs(x = "Median end-to-end WGS runtime", y = NULL) +
  theme_minimal(base_size = 12, base_family = "sans") +
  theme(
    text = element_text(colour = ink),
    axis.text.y = element_text(face = "bold", size = 12, colour = ink),
    axis.title.x = element_text(face = "bold"),
    panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
    plot.margin = margin(18, 22, 18, 18)
  )
save_plot(clock, file.path(output_dir, "repository", "headliner_wgs_clock"), 9.0, 4.4)

# Alternative: concise cohort-value card suitable for presentations.
economic_cards <- economic_bars[economic_bars$cohort_size == 100, ]
value_card <- ggplot(economic_cards,
  aes(x = cohort_label, y = days_returned, fill = cohort_label)) +
  geom_col(width = 0.62, colour = ink, linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.1f days", days_returned)), vjust = -0.5,
    fontface = "bold", size = 4.1, colour = ink) +
  scale_fill_manual(values = cohort_colours, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(x = NULL, y = "Compute time returned across 100 WGS samples") +
  theme_minimal(base_size = 12, base_family = "sans") +
  theme(
    text = element_text(colour = ink), panel.grid = element_blank(),
    axis.text.x = element_text(face = "bold", colour = ink),
    axis.text.y = element_blank(), axis.ticks = element_blank(),
    axis.title.y = element_text(face = "bold"),
    plot.margin = margin(16, 16, 12, 16)
  )
save_plot(value_card, file.path(output_dir, "alternatives", "cohort_100_time_returned"), 8.4, 4.5)

index <- c(
  "# Story-driven final figure set",
  "",
  "This is an additive candidate set. It does not replace or overwrite earlier snapshots.",
  "All plots are vector-native R/ggplot2 outputs and are also exported as PDF and 600 dpi PNG.",
  "",
  "## Manuscript figures",
  "",
  "**Figure 1: VCFcache reuses annotation work without changing the completed record set.** (A) Each input VCF is matched against a recipe-specific cache. Existing annotations are reused, only cache misses are sent to the configured annotator, and both streams are merged and sorted into a complete output VCF. Users may download a ready-to-use recipe cache or a blueprint from Zenodo. (B) Benefit is governed by cache coverage, the cost of the underlying annotation recipe and the number of samples over which any local cache build is amortized.",
  "",
  "**Figure 2: VCFcache returns most WGS annotation time across three native annotator configurations.** Twelve independent GRCh38 genomes (six KPGP and six SGDP) were tested as ACMG SF v3.3 Panel, Twist Human Core Exome WES and WGS inputs with the bundled gnomAD AF ≥1% strategy. VEP used `--everything`; fastVEP and SnpEff used their clinical recipes. Configurations are compared only against their own direct baseline. (A) Distribution of end-to-end runtime remaining after caching; the dashed line marks equal runtime. (B) Median serial throughput calculated as 86,400 seconds divided by median per-sample wall time; open circles denote direct annotation and arrowheads denote VCFcache. Cohort shapes identify independent biological samples. KPGP, Korean Personal Genome Project; SGDP, Simons Genome Diversity Project; VEP, Variant Effect Predictor; WES, whole-exome sequencing; WGS, whole-genome sequencing.",
  "",
  "**Figure 3: Cache coverage defines the theoretical opportunity, while annotation cost determines how much time is returned.** (A) Observed runtime remaining for the 108 matched GRCh38 sample-tool-assay pairs. Dashed curves show the zero-overhead expectation, `1 − f`, for hit rate `f`; point area represents direct annotation time. Fixed lookup, preprocessing and output costs dominate very short panel jobs but become small relative to WGS annotation. (B) One 4.80-million-variant KPGP genome was processed with six controlled VEP loads and bundled gnomAD AF ≥10% and AF ≥1% caches. Solid lines show observed hours returned as pipeline cost increases; dashed lines show the cache-coverage opportunity, `f ×` direct runtime. All 12 outputs passed semantic comparison.",
  "",
  "**Figure 4: A cohort-derived cache repays its build cost rapidly and returns cumulative compute time.** Three held-out real-world cohorts were evaluated with caches built from three disjoint genomes. (A) Effective speedup includes the complete one-time cache-build duration divided across the indicated cohort size; circles mark empirical break-even. (B) Compute days returned after build cost for cohorts of 10, 100 and 1,000 samples, using the measured cohort-specific median direct and cached WGS runtimes. KPGP, Korean Personal Genome Project; PGP, Personal Genome Project; SGDP, Simons Genome Diversity Project.",
  "",
  "## Supplementary figure",
  "",
  "**Supplementary Figure 1: Annotation load increases reuse benefit without a cache-attributable change in output.** (A) Paired SnpEff speedups for core and clinical recipes on the same GRCh37 PGP samples. Lines connect identical samples; the dashed line marks equal runtime. The documented warning-only sample was excluded before inference for WES and WGS. (B) Of 180 direct/cached cells in the final clinical campaign, 176 were eligible and semantically concordant. Four SnpEff cells were excluded because the direct annotator itself changed output with input order; these are not counted as evidence for or against VCFcache. No eligible pair contained a cache-attributable mismatch.",
  "",
  "## Repository figures",
  "",
  "- `headliner_one_day_of_compute`: intuitive serial WGS throughput under each native benchmark recipe.",
  "- `headliner_wgs_clock`: direct-to-cached median WGS wall time with speedup and time returned.",
  "- `cohort_100_time_returned`: presentation-friendly cumulative compute-time view for 100 WGS samples.",
  "",
  "## Interpretation boundary",
  "",
  "The three annotators intentionally retain their native recipes. These panels estimate within-tool VCFcache benefit and must not be interpreted as rankings of annotation completeness or cross-tool semantic equivalence. No 1000 Genomes sample appears in Figures 2 or 3A. Throughput assumes serial continuous processing under the benchmark allocation and excludes cache download/build time; Figure 4 explicitly includes local cohort-cache build time.",
  ""
)
writeLines(index, file.path(output_dir, "FIGURE_INDEX.md"))

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg[[1]]), mustWork = TRUE)
provenance_paths <- c(
  three_path, snpeff_path, campaign_path, complexity_path,
  economics_path, break_even_path, script_path
)
provenance <- data.frame(
  role = c(
    "native assay pairs", "SnpEff paired loads", "final clinical campaign",
    "controlled VEP pipeline campaign", "cohort economics",
    "cohort break-even", "R/ggplot2 renderer"
  ),
  path = provenance_paths,
  md5 = unname(tools::md5sum(provenance_paths)),
  stringsAsFactors = FALSE
)
write_tsv(provenance, file.path(output_dir, "SOURCE_PROVENANCE.tsv"))
writeLines(capture.output(sessionInfo()), file.path(output_dir, "R_SESSION_INFO.txt"))

output_files <- sort(list.files(output_dir, recursive = TRUE, full.names = TRUE))
output_files <- output_files[!file.info(output_files)$isdir]
output_manifest <- data.frame(
  path = sub(paste0("^", normalizePath(output_dir), "/"), "", normalizePath(output_files)),
  size_bytes = file.info(output_files)$size,
  md5 = unname(tools::md5sum(output_files)),
  stringsAsFactors = FALSE
)
write_tsv(output_manifest, file.path(output_dir, "OUTPUT_MANIFEST.tsv"))

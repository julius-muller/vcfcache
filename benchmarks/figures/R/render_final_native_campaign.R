#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop(
    "Usage: render_final_native_campaign.R ",
    "<campaign_metrics.tsv> <external_assay_metrics.tsv> ",
    "<pipeline_complexity.tsv> <output_directory>"
  )
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(scales)
})

campaign_path <- normalizePath(args[[1]], mustWork = TRUE)
assay_path <- normalizePath(args[[2]], mustWork = TRUE)
complexity_path <- normalizePath(args[[3]], mustWork = TRUE)
output_dir <- args[[4]]

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
for (subdir in c(
  "source", "manuscript", "manuscript/panels", "repository",
  "alternatives", "alternatives/panels"
)) {
  dir.create(file.path(output_dir, subdir), recursive = TRUE, showWarnings = FALSE)
}

read_tsv <- function(path) {
  read.delim(
    path, sep = "\t", header = TRUE, stringsAsFactors = FALSE,
    check.names = FALSE, na.strings = c("", "NA")
  )
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

ink <- "#1B2A34"
blue <- "#2B6F92"
orange <- "#D07A32"
green <- "#3D8B62"
purple <- "#7660A7"
red <- "#B64F4F"
grey <- "#69777D"
light_grey <- "#E9EEF0"
grid_grey <- "#DDE3E5"
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

format_seconds <- function(value) {
  vapply(value, function(seconds) {
    if (seconds < 60) {
      sprintf("%.0f s", seconds)
    } else if (seconds < 3600) {
      sprintf("%.0f min", seconds / 60)
    } else {
      sprintf("%.1f h", seconds / 3600)
    }
  }, character(1))
}

campaign <- read_tsv(campaign_path)
pair_keys <- c("sample", "cohort", "tool", "level", "assembly", "assay")
direct <- campaign[campaign$condition == "direct", ]
cached <- campaign[campaign$condition == "cached", ]
pairs <- merge(direct, cached, by = pair_keys, suffixes = c("_direct", "_cached"))
pairs <- pairs[
  pairs$concordant_cached == "PASS" &
    pairs$input_records_direct == pairs$output_records_direct &
    pairs$input_records_cached == pairs$output_records_cached,
]
pairs$direct_seconds <- pairs$wall_seconds_direct
pairs$cached_seconds <- pairs$wall_seconds_cached
pairs$speedup <- pairs$direct_seconds / pairs$cached_seconds
pairs$relative_runtime <- pairs$cached_seconds / pairs$direct_seconds
pairs$time_saved_seconds <- pairs$direct_seconds - pairs$cached_seconds
pairs$hit_rate <- 1 - pairs$missing_variants_cached / pairs$input_records_cached

external <- read_tsv(assay_path)
external <- external[
  external$tool == "vep" & tolower(external$semantic_pass) == "true",
]
vep <- data.frame(
  sample = external$sample,
  cohort = external$cohort,
  tool = "vep",
  level = "everything",
  assembly = external$assembly,
  assay = external$assay,
  direct_seconds = external$uncached_wall_seconds,
  cached_seconds = external$cached_wall_seconds,
  speedup = external$speedup,
  relative_runtime = external$relative_runtime,
  time_saved_seconds = external$wall_seconds_saved,
  hit_rate = external$cache_hit_rate,
  stringsAsFactors = FALSE
)

native <- pairs[
  pairs$assembly == "GRCh38" & pairs$level == "clinical" &
    pairs$tool %in% c("fastvep", "snpeff"),
  names(vep)
]
three <- rbind(vep, native)
three$tool_label <- factor(
  c(vep = "VEP", fastvep = "fastVEP", snpeff = "SnpEff")[three$tool],
  levels = tool_levels
)
three$assay_label <- factor(
  c(panel = "Panel", wes = "WES", wgs = "WGS")[three$assay],
  levels = assay_levels
)
three$cohort_label <- factor(
  toupper(three$cohort), levels = c("KPGP", "SGDP")
)
three$time_saved_percent <- 100 * (1 - three$relative_runtime)

if (nrow(three) != 108) {
  stop("Expected 108 matched three-annotator observations, found ", nrow(three))
}
sample_sets <- split(as.character(three$sample), three$tool_label)
if (!all(vapply(sample_sets[-1], function(x) {
  identical(sort(unique(x)), sort(unique(sample_sets[[1]])))
}, logical(1)))) {
  stop("The three annotator arms do not use identical GRCh38 sample sets")
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
        time_saved_percent = median(values$time_saved_percent),
        hit_rate = median(values$hit_rate)
      )
    )
  })
  do.call(rbind, result)
}

three_summary <- summarize_groups(three, c("tool_label", "assay_label"))
three_summary$tool_label <- factor(three_summary$tool_label, levels = tool_levels)
three_summary$assay_label <- factor(three_summary$assay_label, levels = assay_levels)
three_summary$gain <- ifelse(
  three_summary$speedup >= 1, "Time returned", "Workflow overhead"
)
three_summary$gain <- factor(
  three_summary$gain, levels = c("Time returned", "Workflow overhead")
)
three_summary$tile_label <- ifelse(
  three_summary$time_saved_percent >= 0,
  sprintf("%.2f×\n%.0f%% saved", three_summary$speedup, three_summary$time_saved_percent),
  sprintf("%.2f×\n%.0f%% longer", three_summary$speedup, abs(three_summary$time_saved_percent))
)

write_tsv(three, file.path(output_dir, "source", "three_annotator_assay_pairs.tsv"))
write_tsv(
  three_summary,
  file.path(output_dir, "source", "three_annotator_assay_summary.tsv")
)
file.copy(campaign_path, file.path(output_dir, "source", "campaign_metrics.tsv"))
file.copy(assay_path, file.path(output_dir, "source", "external_assay_v3_metrics.tsv"))
file.copy(
  complexity_path,
  file.path(output_dir, "source", "dual_cache_pipeline_spectrum.tsv")
)

# Figure 1: matched assay-scale results for three native annotator recipes.
p1a <- ggplot(three_summary, aes(y = tool_label, colour = tool_label)) +
  geom_segment(
    aes(x = direct_seconds, xend = cached_seconds, yend = tool_label),
    linewidth = 0.85,
    arrow = grid::arrow(length = grid::unit(0.07, "inches"), type = "closed")
  ) +
  geom_point(aes(x = direct_seconds), shape = 21, fill = "white", size = 2.6, stroke = 0.75) +
  geom_point(aes(x = cached_seconds), shape = 17, size = 2.8) +
  geom_label(
    aes(x = sqrt(direct_seconds * cached_seconds), label = sprintf("%.2f×", speedup)),
    fill = "white", colour = ink, linewidth = 0, size = 2.35,
    show.legend = FALSE
  ) +
  facet_grid(. ~ assay_label) +
  scale_colour_manual(values = tool_colours, guide = "none") +
  scale_x_log10(
    breaks = c(5, 60, 3600), labels = format_seconds,
    limits = c(3.5, 15000)
  ) +
  labs(tag = "A", x = "Median end-to-end runtime", y = NULL) +
  journal_theme(9.3) +
  theme(panel.grid.major.y = element_blank())

p1b <- ggplot(three, aes(x = tool_label, y = speedup, colour = tool_label)) +
  geom_hline(yintercept = 1, colour = grey, linetype = "dashed", linewidth = 0.55) +
  geom_boxplot(width = 0.52, outlier.shape = NA, linewidth = 0.55, alpha = 0.12) +
  geom_point(
    aes(shape = cohort_label), position = position_jitter(width = 0.11, height = 0),
    size = 1.75, alpha = 0.82
  ) +
  facet_grid(. ~ assay_label) +
  scale_colour_manual(values = tool_colours, guide = "none") +
  scale_shape_manual(values = c("KPGP" = 16, "SGDP" = 17)) +
  scale_y_log10(
    breaks = c(0.25, 0.5, 1, 2, 5, 10, 20),
    labels = function(x) paste0(x, "×"), limits = c(0.2, 22)
  ) +
  labs(tag = "B", x = NULL, y = "Speedup over direct annotation") +
  journal_theme(9.3) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))

p1c <- ggplot(three_summary, aes(x = assay_label, y = tool_label, fill = time_saved_percent)) +
  geom_tile(colour = "white", linewidth = 1.2) +
  geom_text(aes(label = tile_label), size = 2.65, lineheight = 0.92, colour = ink) +
  scale_fill_gradient2(
    low = "#D7A19A", mid = paper, high = "#86C19B", midpoint = 0,
    limits = c(-250, 95), oob = squish, guide = "none"
  ) +
  labs(tag = "C", x = NULL, y = NULL) +
  journal_theme(9.3) +
  theme(
    axis.line = element_blank(), axis.ticks = element_blank(),
    panel.grid = element_blank()
  )

figure1 <- p1a / (p1b | p1c) +
  plot_layout(heights = c(0.85, 1.25), widths = c(1.7, 1))

save_plot(
  figure1,
  file.path(output_dir, "manuscript", "figure_assay_three_annotators"),
  7.2, 6.8
)
save_plot(p1a, file.path(output_dir, "manuscript", "panels", "figure_assay_A_runtime"), 7.2, 2.5)
save_plot(p1b, file.path(output_dir, "manuscript", "panels", "figure_assay_B_samples"), 7.2, 3.4)
save_plot(p1c, file.path(output_dir, "manuscript", "panels", "figure_assay_C_summary"), 4.2, 3.4)

# Figure 2: clean within-SnpEff annotation-load contrast on GRCh37.
snpeff <- pairs[
  pairs$tool == "snpeff" & pairs$assembly == "GRCh37" &
    pairs$level %in% c("core", "clinical"),
]
snpeff$recipe_label <- factor(
  c(core = "Core", clinical = "Clinical")[snpeff$level],
  levels = c("Core", "Clinical")
)
snpeff$assay_label <- factor(
  c(panel = "Panel", wes = "WES", wgs = "WGS")[snpeff$assay],
  levels = assay_levels
)
snpeff$time_saved_percent <- 100 * (1 - snpeff$relative_runtime)
snpeff_summary <- summarize_groups(snpeff, c("recipe_label", "assay_label"))
snpeff_summary$recipe_label <- factor(
  snpeff_summary$recipe_label, levels = c("Core", "Clinical")
)
snpeff_summary$assay_label <- factor(
  snpeff_summary$assay_label, levels = assay_levels
)

write_tsv(snpeff, file.path(output_dir, "source", "snpeff_load_pairs.tsv"))
write_tsv(snpeff_summary, file.path(output_dir, "source", "snpeff_load_summary.tsv"))

recipe_colours <- c("Core" = grey, "Clinical" = purple)
p2a <- ggplot(
  snpeff,
  aes(x = recipe_label, y = speedup, group = interaction(sample, assay_label))
) +
  geom_hline(yintercept = 1, colour = grey, linetype = "dashed", linewidth = 0.55) +
  geom_line(colour = "#BBC3C6", linewidth = 0.42, alpha = 0.8) +
  geom_point(aes(colour = recipe_label), size = 2.1, alpha = 0.88) +
  facet_grid(. ~ assay_label) +
  scale_colour_manual(values = recipe_colours, guide = "none") +
  scale_y_log10(
    breaks = c(0.5, 1, 2, 5, 10, 20),
    labels = function(x) paste0(x, "×"), limits = c(0.7, 20)
  ) +
  labs(tag = "A", x = NULL, y = "Speedup over direct SnpEff") +
  journal_theme(9.4)

p2b <- ggplot(snpeff_summary, aes(y = recipe_label, colour = recipe_label)) +
  geom_segment(
    aes(x = direct_seconds, xend = cached_seconds, yend = recipe_label),
    linewidth = 0.9,
    arrow = grid::arrow(length = grid::unit(0.075, "inches"), type = "closed")
  ) +
  geom_point(aes(x = direct_seconds), shape = 21, fill = "white", size = 2.6, stroke = 0.75) +
  geom_point(aes(x = cached_seconds), shape = 17, size = 2.8) +
  geom_label(
    aes(x = sqrt(direct_seconds * cached_seconds), label = sprintf("%.2f×", speedup)),
    fill = "white", colour = ink, linewidth = 0, size = 2.35,
    show.legend = FALSE
  ) +
  facet_grid(. ~ assay_label) +
  scale_colour_manual(values = recipe_colours, guide = "none") +
  scale_x_log10(
    breaks = c(60, 600, 3600), labels = format_seconds,
    limits = c(20, 8000)
  ) +
  labs(tag = "B", x = "Median end-to-end runtime", y = NULL) +
  journal_theme(9.4) +
  theme(panel.grid.major.y = element_blank())

figure2 <- p2a / p2b + plot_layout(heights = c(1.15, 1))
save_plot(
  figure2,
  file.path(output_dir, "manuscript", "figure_snpeff_annotation_load"),
  7.2, 6.2
)
save_plot(p2a, file.path(output_dir, "manuscript", "panels", "figure_snpeff_A_samples"), 7.2, 3.1)
save_plot(p2b, file.path(output_dir, "manuscript", "panels", "figure_snpeff_B_runtime"), 7.2, 2.8)

# Figure 3: controlled VEP pipeline load at two biological cache depths.
complexity <- read_tsv(complexity_path)
complexity <- complexity[tolower(complexity$semantic_pass) == "true", ]
if (nrow(complexity) != 12) {
  stop("Expected 12 validated pipeline-complexity rows")
}
complexity$cache_label <- factor(
  complexity$cache_strategy,
  levels = c("AF ≥ 10% cache", "AF ≥ 1% cache"),
  labels = c("gnomAD AF ≥ 10%", "gnomAD AF ≥ 1%")
)
complexity$direct_hours <- complexity$uncached_wall_seconds / 3600
complexity$cached_hours <- complexity$cached_wall_seconds / 3600
complexity$saved_hours <- complexity$time_saved_seconds / 3600
complexity$delay_label <- factor(
  paste0(complexity$delay_ms, " ms"),
  levels = paste0(sort(unique(complexity$delay_ms), decreasing = TRUE), " ms")
)
cache_colours <- c("gnomAD AF ≥ 10%" = orange, "gnomAD AF ≥ 1%" = green)

p3a <- ggplot(complexity, aes(y = delay_label, colour = cache_label)) +
  geom_segment(
    aes(x = direct_hours, xend = cached_hours, yend = delay_label),
    linewidth = 0.9,
    arrow = grid::arrow(length = grid::unit(0.075, "inches"), type = "closed")
  ) +
  geom_point(aes(x = direct_hours), shape = 21, fill = "white", size = 2.5, stroke = 0.75) +
  geom_point(aes(x = cached_hours), shape = 17, size = 2.8) +
  geom_label(
    aes(x = (direct_hours + cached_hours) / 2, label = sprintf("%.1f h", saved_hours)),
    fill = "white", colour = ink, linewidth = 0, size = 2.3,
    show.legend = FALSE
  ) +
  facet_grid(. ~ cache_label) +
  scale_colour_manual(values = cache_colours, guide = "none") +
  scale_x_continuous(breaks = c(0, 5, 10, 15, 20), limits = c(0, 24)) +
  labs(tag = "A", x = "Observed runtime per genome (hours)", y = "Added work per consequence") +
  journal_theme(9.3)

ceiling <- aggregate(speedup ~ cache_label, complexity, max)
ceiling$hit_rate <- vapply(
  ceiling$cache_label,
  function(label) median(complexity$cache_hit_rate[complexity$cache_label == label]),
  numeric(1)
)
ceiling$theoretical <- 1 / (1 - ceiling$hit_rate)

p3b <- ggplot(
  complexity,
  aes(x = direct_hours, y = speedup, colour = cache_label, group = cache_label)
) +
  geom_hline(
    data = ceiling, aes(yintercept = theoretical, colour = cache_label),
    linetype = "dashed", linewidth = 0.6, show.legend = FALSE
  ) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.3) +
  geom_text(
    data = ceiling,
    aes(
      x = 23.2, y = theoretical,
      label = sprintf("%.1f%% hits: %.1f× ceiling", 100 * hit_rate, theoretical),
      colour = cache_label
    ),
    inherit.aes = FALSE, hjust = 1, vjust = -0.55, size = 2.45,
    show.legend = FALSE
  ) +
  scale_colour_manual(values = cache_colours) +
  scale_x_continuous(breaks = c(5, 10, 15, 20), limits = c(3.5, 24)) +
  scale_y_continuous(breaks = c(2, 4, 6, 8, 10), limits = c(1, 11.2)) +
  labs(tag = "B", x = "Direct pipeline runtime (hours)", y = "Speedup") +
  journal_theme(9.3)

figure3 <- p3a / p3b + plot_layout(heights = c(1.1, 1), guides = "collect") &
  theme(legend.position = "bottom")
save_plot(
  figure3,
  file.path(output_dir, "manuscript", "figure_pipeline_time_returned"),
  7.2, 6.7
)
save_plot(p3a, file.path(output_dir, "manuscript", "panels", "figure_pipeline_A_hours"), 7.2, 3.4)
save_plot(p3b, file.path(output_dir, "manuscript", "panels", "figure_pipeline_B_speedup"), 5.1, 3.5)

# Repository headliner: normalized WGS time, intentionally simple and within-tool.
wgs <- three_summary[three_summary$assay_label == "WGS", ]
wgs$configuration <- factor(
  c(
    VEP = "VEP\n--everything",
    fastVEP = "fastVEP\nclinical",
    SnpEff = "SnpEff\nclinical"
  )[as.character(wgs$tool_label)],
  levels = c("VEP\n--everything", "fastVEP\nclinical", "SnpEff\nclinical")
)
headliner <- rbind(
  data.frame(
    configuration = wgs$configuration, condition = "Direct", percent = 100,
    tool_label = wgs$tool_label, stringsAsFactors = FALSE
  ),
  data.frame(
    configuration = wgs$configuration, condition = "With VCFcache",
    percent = 100 / wgs$speedup, tool_label = wgs$tool_label,
    stringsAsFactors = FALSE
  )
)
headliner$condition <- factor(headliner$condition, levels = c("Direct", "With VCFcache"))
headliner$tool_label <- factor(headliner$tool_label, levels = tool_levels)
headliner$label <- sprintf("%.0f%%", headliner$percent)

headline_plot <- ggplot(headliner, aes(x = condition, y = percent, fill = tool_label)) +
  geom_col(width = 0.58, colour = ink, linewidth = 0.45) +
  geom_text(aes(label = label), vjust = -0.45, fontface = "bold", size = 4.0, colour = ink) +
  geom_label(
    data = wgs,
    aes(
      x = 1.5, y = 73,
      label = sprintf("%.1f× faster\n%s returned", speedup, format_seconds(time_saved_seconds)),
      group = configuration
    ),
    inherit.aes = FALSE, fill = "white", colour = ink,
    linewidth = 0.6, size = 3.2, lineheight = 0.95
  ) +
  facet_grid(. ~ configuration) +
  scale_fill_manual(values = tool_colours, guide = "none") +
  scale_y_continuous(limits = c(0, 112), expand = expansion(mult = c(0, 0))) +
  labs(x = NULL, y = "Relative annotation time") +
  theme_minimal(base_size = 12, base_family = "sans") +
  theme(
    text = element_text(colour = ink),
    panel.grid = element_blank(),
    axis.text.y = element_blank(), axis.ticks = element_blank(),
    axis.text.x = element_text(face = "bold", colour = ink),
    axis.title.y = element_text(face = "bold"),
    strip.text = element_text(face = "bold", size = 13, colour = ink),
    plot.margin = margin(14, 14, 10, 14)
  )

save_plot(
  headline_plot,
  file.path(output_dir, "repository", "headliner_native_wgs"),
  9.5, 4.8
)

# Alternative compact visual: each tool as a direct-to-cached clock line.
clock_plot <- ggplot(wgs, aes(y = tool_label, colour = tool_label)) +
  geom_segment(
    aes(x = direct_seconds / 60, xend = cached_seconds / 60, yend = tool_label),
    linewidth = 5.5, alpha = 0.22, lineend = "round"
  ) +
  geom_point(aes(x = direct_seconds / 60), shape = 21, fill = "white", size = 5, stroke = 1.2) +
  geom_point(aes(x = cached_seconds / 60), shape = 19, size = 5) +
  geom_text(
    aes(
      x = cached_seconds / 60,
      label = sprintf("  %.1f×", speedup)
    ),
    hjust = 0, fontface = "bold", size = 4.2, show.legend = FALSE
  ) +
  scale_colour_manual(values = tool_colours, guide = "none") +
  scale_x_log10(
    breaks = c(2, 5, 15, 30, 60, 120, 180),
    labels = function(x) paste0(x, " min"), limits = c(1.5, 220)
  ) +
  labs(x = "End-to-end time for one held-out genome", y = NULL) +
  theme_minimal(base_size = 12, base_family = "sans") +
  theme(
    text = element_text(colour = ink),
    axis.text.y = element_text(face = "bold", size = 12, colour = ink),
    axis.title.x = element_text(face = "bold"),
    panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
    plot.margin = margin(16, 22, 16, 16)
  )

save_plot(
  clock_plot,
  file.path(output_dir, "alternatives", "native_wgs_clock_lines"),
  8.5, 4.2
)

session <- capture.output(sessionInfo())
writeLines(session, file.path(output_dir, "R_SESSION_INFO.txt"))

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg[[1]]), mustWork = TRUE)
provenance_paths <- c(campaign_path, assay_path, complexity_path, script_path)
provenance <- data.frame(
  role = c(
    "final clinical campaign", "matched VEP assay campaign",
    "controlled VEP pipeline campaign", "R/ggplot2 renderer"
  ),
  path = provenance_paths,
  md5 = unname(tools::md5sum(provenance_paths)),
  stringsAsFactors = FALSE
)
write_tsv(provenance, file.path(output_dir, "SOURCE_PROVENANCE.tsv"))

index <- c(
  "# Final native-annotation figure set",
  "",
  "This snapshot is additive and does not overwrite any prior figure set.",
  "All plots were rendered in R with ggplot2; composite panels use patchwork.",
  "PNG files are 600 dpi, and every manuscript panel is also exported as SVG and PDF.",
  "",
  "## Proposed manuscript figures",
  "",
  "**Proposed Figure 2: Assay scale and recoverable work determine VCFcache benefit across native annotators.** Twelve held-out GRCh38 genomes (six KPGP and six SGDP) were evaluated as matched ACMG SF v3.3 Panel, Twist Human Core Exome WES and WGS inputs with the bundled gnomAD AF ≥1% strategy. VEP used `--everything`; fastVEP and SnpEff used their clinical recipes. Each recipe is compared only with its own direct baseline because annotation content differs across tools. (A) Median direct runtime (open circle), cached runtime (arrowhead) and median paired speedup. (B) Sample-level speedups; colour identifies annotator and shape identifies cohort. The dashed line marks equal runtime. (C) Median paired speedup and wall-time change. All 108 pairs passed their annotator-specific semantic gate. AF, allele frequency; KPGP, Korean Personal Genome Project; SGDP, Simons Genome Diversity Project; VEP, Variant Effect Predictor; WES, whole-exome sequencing; WGS, whole-genome sequencing.",
  "",
  "**Proposed Supplementary Figure: Additional annotation work increases the benefit of reuse within SnpEff.** The same 12 GRCh37 PGP genomes and the same gnomAD AF ≥1% cache membership were processed with core SnpEff or a clinical recipe adding ClinVar and gnomAD exome lookups. (A) Sample-level paired change in speedup. Lines connect the same sample across recipes; the dashed line marks equal runtime. The warning-only pair described in the manuscript was excluded before inference, leaving 11 paired genomes for WES and WGS. (B) Median direct runtime (open circle), cached runtime (arrowhead) and median paired speedup. Annotation load, rather than cache coverage, moved median WES speedup from 1.16- to 3.62-fold and WGS speedup from 4.46- to 12.31-fold. PGP, Personal Genome Project; WES, whole-exome sequencing; WGS, whole-genome sequencing.",
  "",
  "**Proposed Figure 3: Absolute time returned grows while relative speedup approaches the cache-coverage ceiling.** One held-out 4.80-million-variant KPGP genome was processed with six controlled VEP `--everything` loads and bundled gnomAD AF ≥10% and AF ≥1% caches. (A) Direct runtime (open circle), cached runtime (arrowhead), and hours returned per genome. (B) Speedup as direct annotation becomes slower. Dashed lines are the theoretical `1/(1-f)` ceilings for observed hit rate `f`. The synthetic plugin emitted no annotation field, and all 12 cached outputs passed semantic comparison. AF, allele frequency; VEP, Variant Effect Predictor; WGS, whole-genome sequencing.",
  "",
  "## Repository and alternative figures",
  "",
  "- `headliner_native_wgs`: deliberately simple normalized-time view for the repository front page; bars compare each native recipe only with its own direct baseline.",
  "- `native_wgs_clock_lines`: compact direct-to-cached runtime alternative; open circles are direct annotation and filled circles are VCFcache.",
  "",
  "## Interpretation boundary",
  "",
  "The combined annotator panels show within-tool speedups under native recipes. They must not be used to rank annotation completeness or claim cross-tool semantic equivalence. The 12 genomes are independent KPGP/SGDP evaluation samples and are not 1000 Genomes samples.",
  ""
)
writeLines(index, file.path(output_dir, "FIGURE_INDEX.md"))

output_files <- sort(list.files(output_dir, recursive = TRUE, full.names = TRUE))
output_files <- output_files[!file.info(output_files)$isdir]
output_manifest <- data.frame(
  path = sub(paste0("^", normalizePath(output_dir), "/"), "", normalizePath(output_files)),
  size_bytes = file.info(output_files)$size,
  md5 = unname(tools::md5sum(output_files)),
  stringsAsFactors = FALSE
)
write_tsv(output_manifest, file.path(output_dir, "OUTPUT_MANIFEST.tsv"))

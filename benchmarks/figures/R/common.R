vcfcache_r_library <- Sys.getenv(
  "VCFCACHE_R_LIBRARY",
  unset = file.path(
    path.expand("~"), ".local", "lib", "R",
    paste(R.version$major, strsplit(R.version$minor, "\\.")[[1]][1], sep = ".")
  )
)
.libPaths(c(vcfcache_r_library, .libPaths()))

required_packages <- c("ggplot2", "patchwork", "scales", "jsonlite")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop(
    "Missing R packages: ", paste(missing_packages, collapse = ", "),
    ". Run Rscript --vanilla benchmarks/figures/R/setup_packages.R"
  )
}

library(ggplot2)
library(patchwork)

vcf_colors <- c(
  ink = "#18323F",
  blue = "#28536B",
  green = "#43A047",
  light_green = "#A7D8B7",
  orange = "#E98B2A",
  purple = "#7B61A8",
  red = "#C94C4C",
  grey = "#687980",
  light_grey = "#E8EEF0",
  paper = "#FCFDFB"
)

theme_vcfcache <- function(base_size = 11) {
  theme_minimal(base_size = base_size, base_family = "sans") +
    theme(
      text = element_text(color = vcf_colors[["ink"]]),
      plot.title = element_text(face = "bold", size = rel(1.2), margin = margin(b = 5)),
      plot.subtitle = element_text(color = vcf_colors[["grey"]], margin = margin(b = 10)),
      plot.caption = element_text(color = vcf_colors[["grey"]], hjust = 0, size = rel(0.78)),
      axis.title = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      legend.position = "bottom",
      legend.title = element_text(face = "bold"),
      strip.text = element_text(face = "bold"),
      plot.margin = margin(12, 14, 12, 14)
    )
}

read_tsv <- function(path) {
  read.delim(
    path,
    sep = "\t",
    header = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    na.strings = c("", "NA")
  )
}

write_tsv <- function(value, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.table(value, path, sep = "\t", row.names = FALSE, quote = FALSE, na = "")
}

save_plot <- function(plot, prefix, width, height, dpi = 300) {
  dir.create(dirname(prefix), recursive = TRUE, showWarnings = FALSE)

  grDevices::svg(
    paste0(prefix, ".svg"), width = width, height = height,
    pointsize = 11, onefile = TRUE, family = "sans", bg = "white"
  )
  print(plot)
  grDevices::dev.off()

  grDevices::cairo_pdf(
    paste0(prefix, ".pdf"), width = width, height = height,
    pointsize = 11, family = "sans", bg = "white"
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

bootstrap_median_interval <- function(values, seed, replicates = 10000) {
  set.seed(seed)
  estimates <- replicate(
    replicates,
    median(sample(values, length(values), replace = TRUE))
  )
  as.numeric(stats::quantile(estimates, c(0.025, 0.975), names = FALSE))
}

format_duration <- function(seconds) {
  ifelse(
    seconds < 60,
    sprintf("%.0f sec", seconds),
    ifelse(
      seconds < 3600,
      sprintf("%.0f min", seconds / 60),
      ifelse(seconds < 172800, sprintf("%.1f h", seconds / 3600), sprintf("%.1f d", seconds / 86400))
    )
  )
}

preliminary_caption <- function(snapshot) {
  counts <- snapshot$sample_counts
  sprintf(
    "PRELIMINARY SNAPSHOT %s | external WGS %d/%d complete, %d formally revalidated",
    sub("T", " ", sub("\\+00:00$", " UTC", snapshot$created_at)),
    counts$external_completed,
    counts$external_expected,
    counts$external_semantically_validated
  )
}

as_logical_semantic <- function(value) {
  tolower(as.character(value)) %in% c("true", "t", "1")
}

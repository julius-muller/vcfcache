arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) != 2) {
  stop(
    "Usage: Rscript --vanilla benchmarks/figures/R/render_all.R ",
    "<snapshot-directory> <output-directory>"
  )
}

full_arguments <- commandArgs(trailingOnly = FALSE)
file_argument <- full_arguments[grepl("^--file=", full_arguments)]
if (length(file_argument) != 1) {
  stop("Could not determine the render_all.R location")
}
script_dir <- dirname(normalizePath(sub("^--file=", "", file_argument)))

source(file.path(script_dir, "common.R"))
source(file.path(script_dir, "render_assay.R"))
source(file.path(script_dir, "render_external.R"))
source(file.path(script_dir, "render_user_impact.R"))
source(file.path(script_dir, "render_alternatives.R"))

input_dir <- normalizePath(arguments[[1]], mustWork = TRUE)
output_dir <- arguments[[2]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

snapshot <- jsonlite::fromJSON(
  file.path(input_dir, "SNAPSHOT.json"),
  simplifyVector = FALSE
)
if (!snapshot$status %in% c("PRELIMINARY", "FINAL")) {
  stop("Snapshot status must be PRELIMINARY or FINAL")
}

render_assay_figure(input_dir, output_dir, snapshot)
render_external_figure(input_dir, output_dir, snapshot)
render_user_impact_figure(input_dir, output_dir, snapshot)
render_alternative_figures(input_dir, output_dir, snapshot)

session <- capture.output(sessionInfo())
writeLines(session, file.path(output_dir, "R_SESSION_INFO.txt"))
writeLines(
  c(
    paste("snapshot", input_dir),
    paste("rendered_at", format(Sys.time(), tz = "UTC", usetz = TRUE)),
    paste("renderer", normalizePath(file.path(script_dir, "render_all.R"), mustWork = FALSE))
  ),
  file.path(output_dir, "RENDERED_FROM.txt")
)

cat("Rendered ggplot2 figures under", output_dir, "\n")

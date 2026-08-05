local_library <- Sys.getenv(
  "VCFCACHE_R_LIBRARY",
  unset = file.path(
    path.expand("~"), ".local", "lib", "R",
    paste(R.version$major, strsplit(R.version$minor, "\\.")[[1]][1], sep = ".")
  )
)
dir.create(local_library, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(local_library, .libPaths()))

packages <- c("ggplot2", "patchwork", "scales", "jsonlite")
missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) {
  install.packages(missing, repos = "https://cloud.r-project.org", Ncpus = 1)
}

for (package in packages) {
  cat(package, as.character(packageVersion(package)), "\n")
}

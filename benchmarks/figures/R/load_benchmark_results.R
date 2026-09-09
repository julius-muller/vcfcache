# Read one benchmark TSV as a tibble.
load_benchmark_table <- function(path, required = NULL) {
  if (!file.exists(path)) {
    stop("Missing benchmark result file: ", path)
  }

  result <- tibble::as_tibble(
    read.delim(
      path,
      sep = "\t",
      header = TRUE,
      stringsAsFactors = FALSE,
      check.names = FALSE,
      na.strings = c("", "NA"),
      comment.char = "",
      quote = ""
    )
  )

  if (!is.null(required)) {
    missing <- setdiff(required, names(result))
    if (length(missing) > 0) {
      stop(
        "Missing columns in ", path, ": ",
        paste(missing, collapse = ", ")
      )
    }
  }

  result
}

# Read every TSV below output_dir into a named list of tibbles.
load_benchmark_results <- function(output_dir = "benchmarks/figures/output") {
  output_dir <- normalizePath(output_dir, mustWork = TRUE)
  paths <- sort(
    list.files(
      output_dir,
      pattern = "\\.tsv$",
      recursive = TRUE,
      full.names = TRUE,
      ignore.case = TRUE
    )
  )

  if (length(paths) == 0) {
    stop("No TSV benchmark results found under: ", output_dir)
  }

  relative_paths <- substring(paths, nchar(output_dir) + 2L)
  result_names <- tools::file_path_sans_ext(relative_paths)
  result_names <- gsub("[/\\\\]", "__", result_names)
  result_names <- gsub("[^A-Za-z0-9_]+", "_", result_names)
  result_names <- make.unique(result_names, sep = "_")

  results <- stats::setNames(
    lapply(paths, load_benchmark_table),
    result_names
  )
  attr(results, "output_dir") <- output_dir
  attr(results, "paths") <- stats::setNames(paths, result_names)
  results
}

# Return the file path associated with one loaded result table.
benchmark_result_path <- function(results, name) {
  paths <- attr(results, "paths")
  if (is.null(paths) || !name %in% names(paths)) {
    stop("Unknown benchmark result name: ", name)
  }
  unname(paths[[name]])
}

# Print the available result-table names without printing the tables themselves.
print_benchmark_results <- function(results) {
  paths <- attr(results, "paths")
  if (is.null(paths)) {
    stop("results is not a list returned by load_benchmark_results()")
  }
  print(data.frame(name = names(results), path = unname(paths), row.names = NULL))
  invisible(results)
}

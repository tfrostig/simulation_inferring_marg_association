compare_directories <- function(dir1, dir2, extension = "csv") {
  dir1 <- normalizePath(dir1, mustWork = TRUE)
  dir2 <- normalizePath(dir2, mustWork = TRUE)

  # Match only files with given extension
  pattern <- paste0("\\.", extension, "$")
  files1 <- list.files(dir1, pattern = pattern, recursive = TRUE, full.names = TRUE)
  files2 <- list.files(dir2, pattern = pattern, recursive = TRUE, full.names = TRUE)

  # Use relative paths with file.path-safe separator
  rel1 <- substring(files1, nchar(dir1) + 2)  # +2 to skip separator
  rel2 <- substring(files2, nchar(dir2) + 2)

  common_files <- intersect(rel1, rel2)

  results <- data.frame(
    file = common_files,
    identical = logical(length(common_files)),
    stringsAsFactors = FALSE
  )

  for (i in seq_along(common_files)) {
    cat("Comparing:", common_files[i], "\n")
    f1 <- file.path(dir1, common_files[i])
    f2 <- file.path(dir2, common_files[i])

    df1 <- tryCatch(read.csv(f1, stringsAsFactors = FALSE), error = function(e) NULL)
    df2 <- tryCatch(read.csv(f2, stringsAsFactors = FALSE), error = function(e) NULL)

    results$identical[i] <- identical(df1, df2)
  }

  cat("Compared", nrow(results), "files (", sum(results$identical), "identical)\n")
  return(results)
}

# Example
compare_directories("output/results2025-07-22", "output/results2025-07-19")

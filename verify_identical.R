compare_directories <- function(dir1, dir2) {
  # List files in both directories (non-recursive)
  files1 <- list.files(dir1, full.names = TRUE)
  files2 <- list.files(dir2, full.names = TRUE)

  # Get filenames without paths for matching
  names1 <- basename(files1)
  names2 <- basename(files2)

  # Find common filenames
  common_files <- intersect(names1, names2)

  # Initialize result list
  results <- data.frame(
    file = common_files,
    identical = logical(length(common_files)),
    stringsAsFactors = FALSE
  )

  # Compare each common file
  for (i in seq_along(common_files)) {
    f1 <- file.path(dir1, common_files[i])
    f2 <- file.path(dir2, common_files[i])

    # Check if the files are identical
    results$identical[i] <- identical(readBin(f1, "raw", file.info(f1)$size),
                                      readBin(f2, "raw", file.info(f2)$size))
  }

  return(results)
}

# Example usage
# compare_directories("output/results2025-07-19", "output/results2025-07-22")

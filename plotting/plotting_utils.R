# Global constants
FIGURE_HEIGHT <- 6
FIGURE_WIDTH <- 9
TEXT_SIZE <- 12
ITER_COUNT <- 1000

load_results <- function(results_file, has_mean_sd_col = TRUE) {
  if (!file.exists(results_file)) {
    warning("Results file not found: ", results_file)
    return(NULL)
  }
  cat("Loading results from:", basename(results_file), "\n")
  df <- read.csv(results_file)

  if (has_mean_sd_col && "X" %in% names(df)) {
    df <- df %>%
      mutate(mean_sd = case_when(
        grepl('mean', X, ignore.case = TRUE) ~ 'mean',
        grepl('sd', X, ignore.case = TRUE) ~ 'sd',
        TRUE ~ NA_character_
      )) %>%
      select(-X)
    if(any(is.na(df$mean_sd))) {
      warning("Could not determine mean/sd for all rows in ", basename(results_file))
    }
  } else if ("X" %in% names(df)) {
    df <- select(df, -X)
  }
  return(df)
}

rename_methods <- function(df) {
  if (!"Method" %in% names(df)) return(df)

  df %>%
    mutate(Method = as.character(Method),
           Method = case_when(
             Method == 'Oracle.Corrected' ~ 'Variability corrected (Oracle)',
             Method == 'MLE.Corrected'    ~ 'Variability corrected (MLE)',
             Method == 'Corrected'        ~ 'Variability corrected (Empirical)',
             Method == 'Naive'            ~ 'Naive',
             Method == 'Oracle'           ~ 'Full',
             Method == 'empirical'        ~ 'Variability corrected (Empirical)',
             Method == 'gaussian'         ~ 'Variability corrected (Gaussian)',
             Method == 'naive'            ~ 'Naive',
             Method == 'oracle'           ~ 'Full',
             Method == 'naive_forward_selection' ~ 'Standard FS',
             Method == 'adjust_forward_selection' ~ 'Adjusted FS',
             TRUE ~ Method
           ))
}

calculate_se <- function(df, n_iter = ITER_COUNT, is_proportion = TRUE) {
  if (!all(c("value", "sd") %in% names(df))) {
    warning("Missing 'value' or 'sd' column for SE calculation.")
    return(mutate(df, se = NA_real_))
  }

  if (is_proportion) {
    df %>% mutate(se = sqrt(pmax(value * (1 - value), 0) / n_iter))
  } else {
    df %>% mutate(se = sd / sqrt(n_iter))
  }
}

preprocess_psat_data <- function(raw_data) {
  if(is.null(raw_data)) return(NULL)

  n_iter_val <- raw_data$n.iter[1]

  number_to_pass_df <- raw_data %>%
    filter(mean_sd == 'mean') %>%
    distinct(n.org, n.ref, h, rho, p, eff.dim, g.sparse, number.to.pass)

  long_data <- raw_data %>%
    pivot_longer(cols = matches('FDR|power'), names_to = "metric_name") %>%
    mutate(metric_name = gsub('\\.symmetric', '_symmetric', metric_name))

  type_dat <- tryCatch({
    bind_rows(unglue(long_data$metric_name, '{Method}_{PSAT_method}.{Value_type}'))
  }, error = function(e) {
    warning("unglue failed to parse metric_name: ", e$message)
    tibble(Method = NA_character_, PSAT_method = NA_character_, Value_type = NA_character_, .rows = nrow(long_data))
  })

  long_data <- long_data %>%
    bind_cols(., type_dat) %>%
    select(-metric_name) %>%
    pivot_wider(id_cols = c("n.org", "n.ref", "h", "rho", "p",
                            "eff.dim", "g.sparse", "n.iter", "cas.snp",
                            "Method", "PSAT_method", "Value_type"),
                values_from = 'value',
                names_from  = 'mean_sd') %>%
    rename('value' = 'mean') %>%
    left_join(number_to_pass_df, by = c("n.org", "n.ref", "h", "rho", "p", "eff.dim", "g.sparse")) %>%
    rename_methods() %>%
    calculate_se(n_iter = n_iter_val, is_proportion = TRUE)

  return(long_data)
}

preprocess_normal_fs_data <- function(raw_data, metrics = c("FDR", "Power")) {
  if(is.null(raw_data)) return(NULL)

  n_iter_val <- raw_data$n.iter[1]

  mean_data <- raw_data %>% filter(mean_sd == 'mean') %>% select(-mean_sd)
  sd_data <- raw_data %>% filter(mean_sd == 'sd') %>% select(-mean_sd, -any_of(c("value", "sd")))

  processed_list <- list()

  for (metric in metrics) {
    metric_pattern <- paste0(".", metric)

    mean_long <- mean_data %>%
      select(-contains(setdiff(metrics, metric))) %>%
      pivot_longer(cols = contains(metric),
                   names_to = 'Method',
                   values_to = 'value') %>%
      mutate(Method = gsub(metric_pattern, '', Method, fixed = TRUE),
             MetricType = metric)

    sd_long <- sd_data %>%
      select(matches("n.org|n.ref|h|rho|p|eff.dim|g.sparse|cas.snp"), contains(metric)) %>%
      pivot_longer(cols = contains(metric),
                   names_to = 'Method',
                   values_to = 'sd') %>%
      mutate(Method = gsub(metric_pattern, '', Method, fixed = TRUE))

    key_cols <- intersect(names(mean_long), names(sd_long))
    key_cols <- setdiff(key_cols, c("value", "sd", "MetricType"))

    joined_data <- left_join(mean_long, select(sd_long, all_of(key_cols), sd), by = key_cols)
    processed_list[[metric]] <- joined_data
  }

  bind_rows(processed_list) %>%
    rename_methods() %>%
    calculate_se(n_iter = n_iter_val, is_proportion = TRUE)
}

get_plot_colors <- function(methods) {
  full_palette <- c(
    'Full'                                = '#20854EFF',
    'Naive'                               = '#E18727FF',
    'Variability corrected (Gaussian)'    = '#0072B5FF',
    'Variability corrected (Empirical)'   = '#BC3C29FF',
    'Variability corrected (MLE)'         = '#7876B1FF',
    'Variability corrected (Oracle)'      = 'grey50',
    'Standard FS'                         = 'darkred',
    'Adjusted FS'                         = 'darkblue'
  )

  present_colors <- full_palette[names(full_palette) %in% methods]

  missing_methods <- setdiff(methods, names(present_colors))
  if (length(missing_methods) > 0) {
    warning("No color defined for methods: ", paste(missing_methods, collapse=", "))
  }

  return(present_colors)
}

get_axis_label <- function(label_type) {
  switch(label_type,
         "n_org" = TeX('Original sample size, $n_o$'),
         "n_ref" = TeX('Reference sample size, $n_r$'),
         "h" = TeX('Explained variance, $h^2$'),
         "rho" = TeX('LD correlation, $\\rho$'),
         "causal_snps" = TeX('Causal SNPs, $s^*$'),
         stop("Unknown label type requested.")
  )
}

save_plot <- function(plot_obj, figure_dir, file_prefix, plot_type,
                      ext = "pdf", height = FIGURE_HEIGHT, width = FIGURE_WIDTH) {
  ensure_directory(figure_dir)

  filename <- file.path(figure_dir, paste0(file_prefix, "_", plot_type, ".", ext))
  cat("Saving plot to:", basename(filename), "\n")

  ggsave(filename = filename, plot = plot_obj, height = height, width = width)
  invisible(filename)
}


format_string <- function(string_vec, prefix = 'n[r]') {
  new_names <- paste0(prefix, ' = ', string_vec)
  names(new_names) <- string_vec
  return(new_names)
}

ensure_directory <- function(dir_path) {
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE)
    message("Created directory: ", dir_path)
  }
}

save_plot <- function(plot, figure_dir, filename, height = PLOT_HEIGHT, width = PLOT_WIDTH) {
  output_file <- file.path(figure_dir, filename)
  ggsave(output_file, plot, height = height, width = width)
  invisible(output_file)
}

read_data_file <- function(results_file) {
  if (!file.exists(results_file)) {
    warning("Results file not found: ", results_file)
    return(NULL)
  }

  data <- read.csv(results_file, row.names = NULL)
  data <- data %>% mutate(is_mean = grepl("mean", X))
  return(data)
}

calculate_target_fdr <- function(cas_snp) {
  return(0.05 * (MAX_CAUSAL_SNPS - cas_snp) / MAX_CAUSAL_SNPS)
}

create_legend <- function(plot) {
  return(cowplot::get_legend(plot + theme(legend.position = "bottom")))
}

# Common preprocessing steps
prepare_data <- function(res_freq) {
  dat_fd <- res_freq %>%
    select(-contains("Power"), -X) %>%
    pivot_longer(cols = contains("FDR"), names_to = "Method", values_to = "FDR") %>%
    mutate(Method = gsub(".FDR", "", Method))

  dat_power <- res_freq %>%
    select(-contains("FDR"), -X) %>%
    pivot_longer(cols = contains("power"), names_to = "Method", values_to = "Power") %>%
    mutate(Method = gsub(".power", "", Method))

  return(list(dat_fd = dat_fd, dat_power = dat_power))
}

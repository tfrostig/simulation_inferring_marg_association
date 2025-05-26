# create_plots_new.R
# Script defining functions to create plots from simulation results

library(tidyverse)
library(cowplot)
library(latex2exp)
library(gridExtra)
library(unglue)
library(ggsci)

source("plotting/plotting_utils.R")

# Global constants
FIGURE_HEIGHT <- 6
FIGURE_WIDTH <- 9
TEXT_SIZE <- 14
ITER_COUNT <- 1000
PLOT_HEIGHT <- 6.5  # For individual plots
PLOT_WIDTH <- 8     # For individual plots
COMBINED_HEIGHT <- 8.5  # For combined plots
COMBINED_WIDTH <- 10    # For combined plots
MAX_CAUSAL_SNPS <- 20   # Total possible causal SNPs
POINT_SIZE <- 0.9       # Standard point size
LINE_SIZE <- 0.7        # Standard line size
LINE_SIZE_THICK <- 1.0  # Thicker line size for emphasis
ERROR_MULTIPLIER <- 2   # For confidence intervals


# 1. Standard Simulation Plots -----------------------------------------------------
createStandardSimPlots <- function(results_file, figure_dir) {
  data <- read_data_file(results_file)
  if (is.null(data)) {
    warning("Standard simulation results file not found: ", results_file)
    return(NULL)
  }

  ensure_directory(figure_dir)

  # Preprocess data
  data$g.sparse[data$g.sparse == 1] <- "rare"
  data$g.sparse[data$g.sparse == 2] <- "common"

  processed_data <- prepare_data(data)
  dat_fd <- processed_data$dat_fd
  dat_power <- processed_data$dat_power

  # Create color vector
  col_vec <- ggsci::pal_nejm()(4)
  names(col_vec) <- unique(dat_fd$Type)

  # Create FDR plot
  dat_fd$FDR[is.na(dat_fd$FDR)] <- 0
  fd_plot <- dat_fd %>%
    filter(is_mean, eff.dim != "9,10,11", g.sparse == "rare", rho == 0.95) %>%
    mutate(sd = sqrt(FDR * (1 - FDR) / ITER_COUNT)) %>%
    ggplot(aes(x = h, y = FDR, group = Method, color = Method)) +
    geom_point(size = POINT_SIZE, show.legend = FALSE) +
    geom_linerange(aes(ymin = pmax(FDR - ERROR_MULTIPLIER * sd, 0),
                       ymax = FDR + ERROR_MULTIPLIER * sd),
                   show.legend = FALSE) +
    geom_hline(aes(yintercept = calculate_target_fdr(cas.snp)), col = "grey") +
    geom_line(size = LINE_SIZE) +
    facet_grid(eff.dim ~ n.ref, scales = "free") +
    theme_bw(base_size = TEXT_SIZE) +
    scale_color_manual(values = col_vec, labels = c("Variability corrected (Empirical)",
                                                    "Variability corrected (Guassian)", "Naive", "Full")) +
    labs(tag = "A") +
    theme(axis.text.x = element_text(angle = 30, hjust = 0.75), legend.position = "bottom")

  # Create Power plot
  power_plot <- dat_power %>%
    filter(Method != "emperical.gauss", Method != "naive",
           eff.dim != "9,10,11", is_mean, g.sparse == "rare", rho == 0.95) %>%
    mutate(sd = sqrt(Power * (1 - Power) / ITER_COUNT)) %>%
    ggplot(aes(x = h, y = Power, group = Method, color = Method)) +
    geom_point(size = POINT_SIZE, show.legend = FALSE) +
    geom_linerange(aes(ymin = pmax(Power - ERROR_MULTIPLIER * sd, 0),
                       ymax = pmin(Power + ERROR_MULTIPLIER * sd, 1)),
                   show.legend = FALSE) +
    geom_line(size = LINE_SIZE) +
    theme_bw() +
    xlab(TeX("Explained variance, $h$")) +
    facet_grid(eff.dim ~ n.ref, scales = "free") +
    scale_color_manual(values = col_vec[-c(2:3)], labels = c("Variability corrected (Empirical)", "Full")) +
    labs(tag = "B") +
    theme(axis.text.x = element_text(angle = 30), legend.position = "none")

  # Create combined plot
  legend <- create_legend(fd_plot)
  combined_plot <- gridExtra::grid.arrange(
    fd_plot + theme(legend.position = "none"),
    power_plot,
    legend,
    ncol = 1,
    nrow = 3,
    layout_matrix = rbind(c(1,1), c(2,2), c(3,3)),
    heights = c(2.5, 2.5, 0.2)
  )

  # Save plots
  save_plot(combined_plot, figure_dir, "FDR_POWER_Standard_Sim_fig_2.pdf",
            height = COMBINED_HEIGHT, width = COMBINED_WIDTH)

  message("Standard plots saved to: ", figure_dir)
  return(list(fd_plot = fd_plot, power_plot = power_plot, combined = combined_plot))
}

# 2. PSAT Simulation Plots --------------------------------------------------------
createPSATPlots <- function(results_file, figure_dir) {
  if (!file.exists(results_file)) {
    warning("PSAT simulation results file not found: ", results_file, ". Skipping plots.")
    return(NULL)
  }

  ensure_directory(figure_dir)

  full_data <- read.csv(results_file)
  full_data <- full_data %>%
    mutate(mean_sd = case_when(grepl('mean', X) ~ 'mean',
                               grepl('sd', X) ~ 'sd')) %>%
    select(-X)

  number_to_pass_df <- full_data %>%
    filter(mean_sd == 'mean') %>%
    select(n.org, n.ref, h, rho, p, eff.dim, g.sparse, number.to.pass)

  long_data <- prepare_psat_data(full_data, number_to_pass_df)

  # Define colors and filter methods for most plots
  col_vec <- c('Full' = '#20854EFF',
               'Naive'  = '#E18727FF',
               'Variability corrected (Empirical)' = '#BC3C29FF',
               'Variability corrected (MLE)' = '#7876B1FF',
               'Variability corrected (Oracle)' = 'grey50')

  methods_to_keep <- c('Full', 'Variability corrected (Empirical)', 'Variability corrected (MLE)')
  long_data_filtered <- long_data %>% filter(Method %in% methods_to_keep)
  col_vec_filtered <- col_vec[methods_to_keep]

  # Generate all the PSAT plots
  plots <- generate_psat_plots(long_data, long_data_filtered, col_vec, col_vec_filtered, figure_dir)

  message("PSAT plots saved to: ", figure_dir)
  return(plots)
}

# Helper function to prepare PSAT data
prepare_psat_data <- function(full_data, number_to_pass_df) {
  long_data <- full_data %>% pivot_longer(cols = matches('FDR|power'))
  long_data <- long_data %>% mutate(name = gsub('\\.symmetric', '_symmetric', name))
  type_dat  <- bind_rows(unglue::unglue(long_data$name, '{Method}_{PSAT_method}.{Value_type}'))
  long_data <- long_data %>% bind_cols(., type_dat)

  long_data <- long_data %>%
    pivot_wider(id_cols = c("n.org", "n.ref", "h", "rho", "p",
                            "eff.dim", "g.sparse", "n.iter", "cas.snp",
                            "name", "Method", "PSAT_method", "Value_type"),
                values_from = 'value',
                names_from  = 'mean_sd') %>%
    rename('value' = 'mean') %>%
    left_join(number_to_pass_df) %>%
    mutate(
      se = sd / sqrt(ITER_COUNT),
      se_prop = sqrt(pmax(value * (1 - value), 0) / ITER_COUNT),
      Method = recode(Method,
                      'Oracle.Corrected' = 'Variability corrected (Oracle)',
                      'MLE.Corrected'    = 'Variability corrected (MLE)',
                      'Corrected'        = 'Variability corrected (Empirical)',
                      'Naive'            = 'Naive',
                      'Oracle'           = 'Full')
    )

  return(long_data)
}

# Helper function to generate all PSAT plots
generate_psat_plots <- function(long_data, long_data_filtered, col_vec, col_vec_filtered, figure_dir) {
  # Plot A: Compare FDR (Naive vs. Corrected)
  plotA_data <- long_data %>%
    filter(n.ref %in% c(500, 2000), rho %in% c(0.7, 0.95),
           eff.dim == '1,20', PSAT_method %in% c('naive', 'polyhedral_symmetric'),
           Value_type == 'FDR',
           Method %in% c('Variability corrected (Empirical)', 'Variability corrected (MLE)', 'Naive', 'Full'))
  plotA_cols <- col_vec[intersect(names(col_vec), unique(plotA_data$Method))]

  pA <- create_psat_plot_a(plotA_data, plotA_cols)
  save_plot(pA, figure_dir, "A_PSAT_ECCM_FDR_Compare.pdf",
            height = FIGURE_HEIGHT, width = FIGURE_WIDTH)

  # Plot B: FDR (Selected Methods, PSAT corrected)
  pB <- create_psat_plot_b(long_data_filtered, col_vec_filtered)
  save_plot(pB, figure_dir, "B_PSAT_ECCCM_FDR_Sim.pdf",
            height = FIGURE_HEIGHT, width = FIGURE_WIDTH)

  # Plot C: Conditional Power
  pC <- create_psat_plot_c(long_data_filtered, col_vec_filtered)
  save_plot(pC, figure_dir, "C_PSAT_ECCCM_Power_Cond.pdf",
            height = FIGURE_HEIGHT, width = FIGURE_WIDTH)

  # Plot D: Unconditional Power
  pD_data <- calculate_unconditional_power(long_data_filtered)
  pD <- create_psat_plot_d(pD_data, col_vec_filtered)
  save_plot(pD, figure_dir, "D_PSAT_ECCCM_Power_Uncond.pdf",
            height = FIGURE_HEIGHT, width = FIGURE_WIDTH)

  # Combined Plot
  pComb <- create_psat_combined_plot(long_data, col_vec)
  save_plot(pComb, figure_dir, "E_PSAT_ECCCM_Combined_fig_3.pdf",
            height = FIGURE_HEIGHT, width = FIGURE_WIDTH)

  return(list(
    plot_a = pA,
    plot_b = pB,
    plot_c = pC,
    plot_d = pD,
    plot_combined = pComb
  ))
}

create_psat_plot_a <- function(plot_data, colors) {
  ggplot(plot_data, aes(x = h, y = value, color = Method, linetype = PSAT_method)) +
    geom_hline(aes(yintercept = calculate_target_fdr(cas.snp)), col = 'grey', linetype="dashed") +
    geom_line(size = LINE_SIZE) +
    geom_point(size = POINT_SIZE, show.legend = FALSE) +
    geom_linerange(aes(ymin = pmax(value - ERROR_MULTIPLIER * se_prop, 0),
                       ymax = pmin(value + ERROR_MULTIPLIER * se_prop, 1)),
                   show.legend = FALSE) +
    facet_grid(eff.dim + rho ~ n.ref, scales = 'free') +
    ylab('FDR') + theme_bw() + ylim(c(0, 0.3)) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1),
          legend.position = 'bottom', legend.box = "vertical") +
    scale_color_manual(values = colors) +
    scale_linetype_manual(name = 'PSAT',
                          labels = c('naive' = 'No selection correction',
                                     'polyhedral_symmetric' = 'Selection correction'),
                          values = c('solid', 'dashed'))
}

create_psat_plot_b <- function(data_filtered, colors_filtered) {
  data_filtered %>%
    filter(n.ref %in% c(500, 2000), rho %in% c(0.7, 0.95),
           PSAT_method == 'polyhedral_symmetric', eff.dim == '1,20', Value_type == 'FDR') %>%
    ggplot(aes(x = h, y = value, color = Method)) +
    geom_hline(aes(yintercept = calculate_target_fdr(cas.snp)), color = 'grey', linetype="dashed") +
    geom_line(size = LINE_SIZE) +
    geom_point(size = POINT_SIZE, show.legend = FALSE) +
    geom_linerange(aes(ymin = pmax(value - ERROR_MULTIPLIER * se_prop, 0),
                       ymax = pmin(value + ERROR_MULTIPLIER * se_prop, 1)),
                   show.legend = FALSE) +
    facet_grid(eff.dim + rho ~ n.ref, scales = 'free') +
    ylab('FDR') + theme_bw() + ylim(c(0, 0.15)) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = 'bottom') +
    scale_color_manual(values = colors_filtered)
}

create_psat_plot_c <- function(data_filtered, colors_filtered) {
  data_filtered %>%
    filter(n.ref %in% c(500, 2000), rho %in% c(0.7, 0.95),
           PSAT_method == 'polyhedral_symmetric', eff.dim == '1,20', Value_type == 'power') %>%
    ggplot(aes(x = h, y = value, color = Method)) +
    geom_line(size = LINE_SIZE) +
    geom_point(size = POINT_SIZE, show.legend = FALSE) +
    geom_linerange(aes(ymin = pmax(value - ERROR_MULTIPLIER * se_prop, 0),
                       ymax = pmin(value + ERROR_MULTIPLIER * se_prop, 1)),
                   show.legend = FALSE) +
    ylim(c(0.5, 1)) + facet_grid(eff.dim + rho ~ n.ref, scales = 'free') +
    ylab('Conditional power') + theme_bw() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = 'bottom') +
    scale_color_manual(values = colors_filtered)
}

calculate_unconditional_power <- function(data_filtered) {
  data_filtered %>%
    filter(n.ref %in% c(500, 2000), rho %in% c(0.7, 0.85, 0.95),
           PSAT_method == 'polyhedral_symmetric', eff.dim == '1,20', Value_type == 'power') %>%
    mutate(uncond_power = value / (number.to.pass + 1),
           se_uncond = se_prop / (number.to.pass + 1))
}

create_psat_plot_d <- function(pd_data, colors_filtered) {
  ggplot(pd_data, aes(x = h, y = uncond_power, color = Method)) +
    geom_line(size = LINE_SIZE) +
    geom_point(size = POINT_SIZE, show.legend = FALSE) +
    geom_linerange(aes(ymin = pmax(uncond_power - ERROR_MULTIPLIER * se_uncond, 0),
                       ymax = pmin(uncond_power + ERROR_MULTIPLIER * se_uncond, 1)),
                   show.legend = FALSE) +
    ylim(c(0, 1)) + facet_grid(eff.dim + rho ~ n.ref, scales = 'free') +
    ylab('Unconditional power') + theme_bw() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = 'bottom') +
    scale_color_manual(values = colors_filtered)
}

create_psat_combined_plot <- function(long_data, col_vec) {
  plot_dat <- bind_rows(
    long_data %>% filter(n.ref == 1000, rho %in% c(0.7, 0.8, 0.95), eff.dim == '1,20',
                         PSAT_method %in% c('naive', 'polyhedral_symmetric'), Value_type == 'FDR',
                         Method %in% names(col_vec)),
    long_data %>% filter(n.ref == 1000, rho %in% c(0.7, 0.8, 0.95),
                         PSAT_method == 'polyhedral_symmetric', eff.dim == '1,20', Value_type == 'power',
                         Method %in% names(col_vec)) %>%
      mutate(Value_type = 'Conditional power'),
    long_data %>% filter(n.ref == 1000, rho %in% c(0.7, 0.8, 0.95),
                         PSAT_method == 'polyhedral_symmetric', eff.dim == '1,20', Value_type == 'power',
                         Method %in% names(col_vec)) %>%
      mutate(value = value / (number.to.pass + 1), Value_type = 'Unconditional power')
  ) %>%
    mutate(Value_type = recode(Value_type, 'FDR' = 'Conditional FDR'),
           PSAT_method = factor(PSAT_method),
           Value_type = factor(Value_type, levels = c('Conditional FDR', 'Conditional power', 'Unconditional power')),
           se_combined = case_when(
             Value_type == 'Conditional FDR' ~ sqrt(pmax(value * (1 - value), 0) / ITER_COUNT),
             Value_type == 'Conditional power' ~ sqrt(pmax(value * (1 - value), 0) / ITER_COUNT),
             Value_type == 'Unconditional power' ~ (sqrt(pmax(value * (number.to.pass + 1) *
                                                                (1 - value * (number.to.pass + 1)), 0) / ITER_COUNT)) / (number.to.pass + 1),
             TRUE ~ NA_real_
           )
    )

  combined_methods <- intersect(names(col_vec), unique(plot_dat$Method))
  combined_cols <- col_vec[combined_methods]

  ggplot(plot_dat, aes(x = h, y = value, color = Method, linetype = PSAT_method)) +
    geom_hline(data = . %>% filter(Value_type == 'Conditional FDR'), aes(yintercept = 0.05),
               color = 'grey', linetype="dashed") +
    geom_line(size = LINE_SIZE) +
    geom_point(size = POINT_SIZE, show.legend = FALSE) +
    geom_linerange(aes(ymin = pmax(value - ERROR_MULTIPLIER * se_combined, 0),
                       ymax = pmin(value + ERROR_MULTIPLIER * se_combined, 1)),
                   show.legend = FALSE, na.rm = TRUE) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1),
          legend.position = 'bottom', legend.box = "vertical") +
    scale_color_manual(values = combined_cols) +
    scale_linetype_manual(name = 'Adjustment for selection',
                          labels = c('naive' = 'No selection adjustment',
                                     'polyhedral_symmetric' = 'Selection adjustment'),
                          values = c('dashed', 'solid')) +
    facet_grid(Value_type ~ rho, scales = 'free') +
    scale_x_continuous(sec.axis = sec_axis(~ . , name = TeX("$\\rho$"), breaks = NULL, labels = NULL)) +
    ylab('')
}

# 3. Normal Simulation Plots ------------------------------------------------------
createNormalPlots <- function(results_file, figure_dir) {
  data <- read_data_file(results_file)
  if (is.null(data)) {
    warning("Normal simulation results file not found: ", results_file, ". Skipping plots.")
    return(NULL)
  }

  ensure_directory(figure_dir)

  processed_data <- prepare_data(data)
  dat_fd <- processed_data$dat_fd
  dat_power <- processed_data$dat_power

  # Create color vector and names
  n_ref_names <- c(`500` = 'n_r = 500', `1000` = 'n_r = 1000', `5000` = 'n_r = 5000')
  n_org_names <- c(`1000` = 'n_o = 1000', `5000` = 'n_o = 5000',
                   `10000` = 'n_o = 10000', `15000` = 'n_o = 15000')

  col_vec <- ggsci::pal_nejm()(4)
  names(col_vec) <- unique(dat_fd$Method)

  dat_fd$FDR[is.na(dat_fd$FDR)] <- 0

  # Create and save plots
  fd_plot <- create_normal_fd_plot(dat_fd, col_vec)
  save_plot(fd_plot, figure_dir, "FDR_GAUSSIAN_fig_1.pdf")

  fd_no_naive_plot <- create_normal_fd_no_naive_plot(dat_fd, col_vec)
  save_plot(fd_no_naive_plot, figure_dir, "FDR_GAUSS_NO_NAIVE_fig_s1.pdf")

  power_plot <- create_normal_power_plot(dat_power, col_vec)
  save_plot(power_plot, figure_dir, "Power_GAUSS_fig_s2.pdf")

  combined_plot <- create_normal_combined_plot(fd_plot, power_plot, figure_dir)

  message("Normal plots saved to: ", figure_dir)
  return(list(
    fd_plot = fd_plot,
    fd_no_naive_plot = fd_no_naive_plot,
    power_plot = power_plot,
    combined = combined_plot
  ))
}

create_normal_fd_plot <- function(dat_fd, col_vec) {
  dat_fd %>%
    filter(is_mean) %>%
    filter(h %in% c(0.005, 0.01, 0.05), rho == 0.8) %>%
    mutate(sd = sqrt(FDR * (1 - FDR) / ITER_COUNT)) %>%
    ggplot(aes(x = n.org, y = FDR, group = Method, color = Method)) +
    geom_point(size = POINT_SIZE, show.legend = FALSE) +
    geom_linerange(aes(ymin = pmax(FDR - ERROR_MULTIPLIER * sd, 0),
                       ymax = FDR + ERROR_MULTIPLIER * sd),
                   show.legend = FALSE) +
    geom_hline(aes(yintercept = calculate_target_fdr(cas.snp)), col = 'grey') +
    geom_line(size = LINE_SIZE) +
    xlab(TeX('Original sample size, $n_o$')) +
    theme_bw() +
    facet_grid(h ~ n.ref, scales = 'free') +
    theme_bw(base_size = TEXT_SIZE) +
    scale_color_manual(labels = c('Variability corrected (Empirical)', 'Variability corrected (Guassian)',
                                  'Naive', 'Full'),
                       values = col_vec) +
    theme(axis.text.x = element_text(angle = 30, hjust=0.75), legend.position = 'bottom') +
    scale_x_continuous(sec.axis = sec_axis(~ . , name = TeX("$n_r$"), breaks = NULL, labels = NULL)) +
    scale_y_continuous(sec.axis = sec_axis(~ . , name = TeX("Explained variance, $h$"),
                                           breaks = NULL, labels = NULL))
}

create_normal_fd_no_naive_plot <- function(dat_fd, col_vec) {
  dat_fd %>%
    filter(is_mean) %>%
    filter(h %in% c(0.005, 0.01, 0.05), rho == 0.8, Method != 'naive') %>%
    mutate(sd = sqrt(FDR * (1 - FDR) / ITER_COUNT)) %>%
    ggplot(aes(x = n.org, y = FDR, group = Method, color = Method)) +
    geom_point(size = POINT_SIZE, show.legend = FALSE) +
    geom_linerange(aes(ymin = pmax(FDR - ERROR_MULTIPLIER * sd, 0),
                       ymax = FDR + ERROR_MULTIPLIER * sd),
                   show.legend = FALSE) +
    geom_hline(aes(yintercept = calculate_target_fdr(cas.snp)), col = 'grey') +
    geom_line(size = LINE_SIZE) +
    xlab(TeX('Original sample size, $n_o$')) +
    theme_bw() +
    facet_grid(h ~ n.ref, scales = 'free') +
    theme_bw(base_size = TEXT_SIZE) +
    scale_color_manual(labels = c('Variability corrected (Empirical)',
                                  'Variability corrected (Guassian)', 'Full'),
                       values = col_vec) +
    theme(axis.text.x = element_text(angle = 30, hjust=0.75), legend.position = 'bottom') +
    scale_x_continuous(sec.axis = sec_axis(~ . , name = TeX("$n_r$"), breaks = NULL, labels = NULL)) +
    scale_y_continuous(sec.axis = sec_axis(~ . , name = TeX("Explained variance, $h$"),
                                           breaks = NULL, labels = NULL))
}

create_normal_power_plot <- function(dat_power, col_vec) {
  dat_power %>%
    filter(is_mean) %>%
    filter(Method != 'naive') %>%
    filter(h <= 0.005, rho == 0.8, n.org <= 50000) %>%
    mutate(sd = sqrt(Power * (1 - Power) / ITER_COUNT)) %>%
    ggplot(aes(x = n.org, y = Power, group = Method, color = Method)) +
    geom_point(size = POINT_SIZE, show.legend = FALSE) +
    geom_linerange(aes(ymin = pmax(Power - ERROR_MULTIPLIER * sd, 0),
                       ymax = pmin(Power + ERROR_MULTIPLIER * sd, 1)),
                   show.legend = FALSE) +
    geom_line(size = LINE_SIZE_THICK) +
    theme_bw() +
    xlab('Original sample size') +
    facet_grid(h ~ n.ref, scales = 'free') +
    scale_color_manual(labels = c('Variability corrected (Empirical)',
                                  'Variability corrected (Guassian)', 'Full'),
                       values = col_vec) +
    theme(axis.text.x = element_text(angle = 30), legend.position = 'bottom') +
    scale_x_continuous(sec.axis = sec_axis(~ . , name = TeX("$n_r$"), breaks = NULL, labels = NULL)) +
    scale_y_continuous(sec.axis = sec_axis(~ . , name = TeX("Explained variance, $h$"),
                                           breaks = NULL, labels = NULL))
}

create_normal_combined_plot <- function(fd_plot, power_plot, figure_dir) {
  legend_a <- create_legend(fd_plot)
  full_plot <- gridExtra::grid.arrange(
    fd_plot + theme(legend.position = 'none'),
    power_plot + theme(legend.position = 'none'),
    legend_a,
    ncol = 1,
    nrow = 3,
    layout_matrix = rbind(c(1,1), c(2,2), c(3,3)),
    heights = c(2.5, 2.5, 0.2)
  )

  save_plot(full_plot, figure_dir, "FDR_POWER_GAUSS_ALL.pdf",
            height = COMBINED_HEIGHT, width = COMBINED_WIDTH)

  # Secondary save with different parameters for compatibility
  cowplot::save_plot(file.path(figure_dir, "FDR_POWER_GAUSS_ALL.pdf"),
                     full_plot, base_width = 7, base_height = 9)

  return(full_plot)
}

# 4. Forward Selection Plots ------------------------------------------------------
createForwardSelectionPlots <- function(results_file, figure_dir) {
  data <- read_data_file(results_file)
  if (is.null(data)) {
    warning("Forward selection results file not found: ", results_file, ". Skipping plots.")
    return(NULL)
  }

  ensure_directory(figure_dir)

  processed_data <- prepare_data(data)
  dat_fd <- processed_data$dat_fd
  dat_power <- processed_data$dat_power

  col_vec <- ggsci::pal_nejm()(4)
  names(col_vec) <- unique(dat_fd$Type)

  dat_fd$FDR[is.na(dat_fd$FDR)] <- 0

  fd_plot <- create_fs_fd_plot(dat_fd)
  power_plot <- create_fs_power_plot(dat_power)
  combined_plot <- create_fs_combined_plot(fd_plot, power_plot, figure_dir)

  message("Forward selection plots saved to: ", figure_dir)
  return(list(fd_plot = fd_plot, power_plot = power_plot, combined = combined_plot))
}

create_fs_fd_plot <- function(dat_fd) {
  dat_fd %>%
    filter(is_mean) %>%
    filter(eff.dim != '9,10,11') %>%
    filter(g.sparse == 'rare', rho == 0.95) %>%
    mutate(sd = sqrt(FDR * (1 - FDR) / ITER_COUNT)) %>%
    ggplot(aes(x = h, y = FDR, group = Method, color = Method)) +
    geom_point(size = POINT_SIZE, show.legend = FALSE) +
    geom_linerange(aes(ymin = pmax(FDR - ERROR_MULTIPLIER * sd, 0),
                       ymax = FDR + ERROR_MULTIPLIER * sd),
                   show.legend = FALSE) +
    geom_hline(aes(yintercept = calculate_target_fdr(cas.snp)), col = 'grey') +
    geom_line(size = LINE_SIZE_THICK) +
    xlab('') +
    theme_bw() +
    facet_grid(eff.dim ~ n.ref, scales = 'free') +
    theme_bw(base_size = TEXT_SIZE) +
    scale_color_manual(labels = c('Standard FS', 'Adjusted FS'),
                       values = c('darkred', 'darkblue')) +
    theme(axis.text.x = element_text(angle = 30, hjust=0.75), legend.position = 'bottom') +
    scale_x_continuous(sec.axis = sec_axis(~ . , name = TeX("$n_r$"), breaks = NULL, labels = NULL)) +
    scale_y_continuous(sec.axis = sec_axis(~ . , name = TeX("Causal SNPs, $s^*$"),
                                           breaks = NULL, labels = NULL))
}

create_fs_power_plot <- function(dat_power) {
  dat_power %>%
    filter(is_mean) %>%
    filter(eff.dim != '9,10,11') %>%
    filter(g.sparse == 'rare', rho == 0.95) %>%
    mutate(sd = sqrt(Power * (1 - Power) / ITER_COUNT)) %>%
    ggplot(aes(x = h, y = Power, group = Method, color = Method)) +
    geom_point(size = POINT_SIZE, show.legend = FALSE) +
    geom_linerange(aes(ymin = pmax(Power - ERROR_MULTIPLIER * sd, 0),
                       ymax = pmin(Power + ERROR_MULTIPLIER * sd, 1)),
                   show.legend = FALSE) +
    geom_hline(aes(yintercept = calculate_target_fdr(cas.snp)), col = 'grey') +
    geom_line(size = LINE_SIZE_THICK) +
    xlab('') +
    theme_bw() +
    facet_grid(eff.dim ~ n.ref, scales = 'free') +
    theme_bw(base_size = TEXT_SIZE) +
    scale_color_manual(labels = c('Adjusted FS', 'Standard FS'),
                       values = c('darkblue', 'darkred')) +
    theme(axis.text.x = element_text(angle = 30, hjust=0.75), legend.position = 'bottom') +
    scale_x_continuous(sec.axis = sec_axis(~ . , name = TeX("$n_r$"), breaks = NULL, labels = NULL)) +
    scale_y_continuous(sec.axis = sec_axis(~ . , name = TeX("Causal SNPs, $s^*$"),
                                           breaks = NULL, labels = NULL))
}

create_fs_combined_plot <- function(fd_plot, power_plot, figure_dir) {
  legend_a <- create_legend(fd_plot)
  full_plot <- cowplot::plot_grid(
    fd_plot + theme(legend.position = 'none'),
    power_plot + theme(legend.position = 'none'),
    legend_a,
    nrow = 3, labels = c('A', 'B', ''),
    rel_heights = c(2.5, 2.5, 0.2)
  )

  save_plot(full_plot, figure_dir, "FS_FDR_power_fig_s6.pdf",
            height = COMBINED_HEIGHT, width = COMBINED_WIDTH)

  return(full_plot)
}

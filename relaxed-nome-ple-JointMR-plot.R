library(ggplot2)
library(RColorBrewer)
library(patchwork)
library(dplyr)

source('F://wsj-doctor//2021-跨种族MR//ago//20230927--no pleo//QQ.plot_20230926.R')
source('F://wsj-doctor//2021-跨种族MR//ago//20230927--no pleo//histplot_20230926.R')

if (basename(getwd()) != "weak" && dir.exists("weak")) {
  setwd("weak")
}
if (dir.exists("F://wsj-doctor//2024-MR-META//weak")) {
  setwd("F://wsj-doctor//2024-MR-META//genemo一修//simulation addition//NOME")
}

#### relaxed-NOME + pleiotropy plugin + delta-method simulation plots ####

gb <- c(200)
beta0b <- c(0, 0.05)
taub <- c(0.2, 0)
rho_1b <- c(0.9, 0.5)
rho_2b <- rho_1b
exposure_rhob <- c(0.9, 0.5)
xy_rhob <- c(0, 0.05)
siga_grid <- data.frame(siga1 = c(0, 0), siga2 = c(0, 0.05))
nxb <- c(3)
nyb <- c(5)
weakp <- 0
F_THRESHOLD <- 10
P_THRESHOLD <- 5e-08

stop_if_missing <- TRUE

build_ccname <- function(nx, ny) {
  methods <- c("dIVW", "IVW_fixed", "IVW_random", "wme")
  pair_count <- nx * ny
  all_beta_names <- unlist(lapply(methods, function(method) {
    paste0("MR_", method, "_", rep("beta", pair_count), seq_len(pair_count))
  }))
  all_se_names <- unlist(lapply(methods, function(method) {
    paste0("MR_", method, "_", rep("se", pair_count), seq_len(pair_count))
  }))
  all_pval_names <- unlist(lapply(methods, function(method) {
    paste0("MR_", method, "_", rep("pval", pair_count), seq_len(pair_count))
  }))

  c(
    "nx", "ny", "g", "beta0", "rho_mean", "exposure_rho", "xy_rho",
    "tau", "siga1", "siga2", "SNP_gaws_meta", "SNPMR", "SNPnew",
    paste0("meta_res_dIVW_", c("beta", "se", "pval")),
    paste0("meta_res_IVW_fixed_", c("beta", "se", "pval")),
    paste0("meta_res_IVW_random_", c("beta", "se", "pval")),
    paste0("meta_res_wme_", c("beta", "se", "pval")),
    all_beta_names, all_se_names, all_pval_names,
    paste0("MR_meta_old_dIVW_", c("beta", "se", "pval")),
    paste0("MR_meta_old_ivw_fixed_", c("beta", "se", "pval")),
    paste0("MR_meta_old_ivw_random_", c("beta", "se", "pval")),
    paste0("MR_meta_old_wme_", c("beta", "se", "pval")),
    paste0("JointMR_RN_", c("beta", "se", "pval")),
    paste0("JointMR_RN_ple_plugin_", c("beta", "se", "pval")),
    "F_tat", "weakp", "alpha0_hat_mean", "alpha0_se_mean",
    "alpha0_true_mean", paste0("alpha0_hat_Y", seq_len(ny)),
    paste0("alpha0_se_Y", seq_len(ny)),
    paste0("alpha0_true_Y", seq_len(ny))
  )
}

keep_names <- function() {
  c(
    "nx", "ny", "g", "beta0", "rho_mean", "exposure_rho", "xy_rho",
    "tau", "siga1", "siga2", "SNP_gaws_meta", "SNPMR", "SNPnew",
    paste0("meta_res_dIVW_", c("beta", "se", "pval")),
    paste0("meta_res_IVW_fixed_", c("beta", "se", "pval")),
    paste0("meta_res_IVW_random_", c("beta", "se", "pval")),
    paste0("meta_res_wme_", c("beta", "se", "pval")),
    paste0("MR_meta_old_dIVW_", c("beta", "se", "pval")),
    paste0("MR_meta_old_ivw_fixed_", c("beta", "se", "pval")),
    paste0("MR_meta_old_ivw_random_", c("beta", "se", "pval")),
    paste0("MR_meta_old_wme_", c("beta", "se", "pval")),
    paste0("JointMR_RN_", c("beta", "se", "pval")),
    paste0("JointMR_RN_ple_plugin_", c("beta", "se", "pval")),
    "F_tat", "weakp", "alpha0_hat_mean", "alpha0_se_mean",
    "alpha0_true_mean"
  )
}

result_filename <- function(nx, ny, g, beta0, rho_1, exposure_rho, xy_rho,
                            tau, siga1, siga2, weakp,
                            F_THRESHOLD, P_THRESHOLD) {
  paste0("Sum_nx_", nx, "_ny_", ny, "_g_", g,
         "_beta0_", beta0, "_rho1_", rho_1,
         "_XXrho_", exposure_rho, "_XYrho_", xy_rho,
         "_tau_", tau, "_siga1_", siga1, "_siga2_", siga2,
         "_weakp_", weakp, "_F_THRESHOLD_", F_THRESHOLD,
         "_P_THRESHOLD_", P_THRESHOLD,
         "_0502_relaxed_nome_ple_plugin_delta.csv")
}

find_result_file <- function(nx, ny, g, beta0, rho_1, exposure_rho, xy_rho,
                             tau, siga1, siga2, weakp,
                             F_THRESHOLD, P_THRESHOLD) {
  filename <- result_filename(nx, ny, g, beta0, rho_1, exposure_rho, xy_rho,
                              tau, siga1, siga2, weakp,
                              F_THRESHOLD, P_THRESHOLD)
  if (file.exists(filename)) {
    return(filename)
  }

  pattern <- paste0("Sum_nx_", nx, "_ny_", ny, "_g_", g,
                    "_beta0_", beta0, "_rho1_", rho_1,
                    "_XXrho_", exposure_rho, "_XYrho_", xy_rho,
                    "_tau_", tau, "_siga1_", siga1, "_siga2_", siga2,
                    "_weakp_", weakp, "_F_THRESHOLD_", F_THRESHOLD,
                    "_P_THRESHOLD_", P_THRESHOLD,
                    "_0428_relaxed_nome_ple_plugin_delta*.csv")
  matches <- Sys.glob(pattern)
  if (length(matches) == 1) {
    return(matches)
  }

  NA_character_
}

read_result_csv <- function(filename, expected_cols) {
  field_counts <- count.fields(filename, sep = ",", quote = "\"",
                               blank.lines.skip = TRUE, comment.char = "")
  keep_rows <- !is.na(field_counts) & field_counts == expected_cols

  if (any(!keep_rows)) {
    message(sprintf("Dropped %s malformed row(s) from %s: expected %s columns.",
                    sum(!keep_rows), basename(filename), expected_cols))
  }

  if (!any(keep_rows)) {
    return(data.frame(matrix(nrow = 0, ncol = expected_cols)))
  }

  lines <- readLines(filename, warn = FALSE)
  lines <- lines[nzchar(trimws(lines))]
  once <- read.csv(text = paste(lines[keep_rows], collapse = "\n"),
                   header = FALSE, quote = "\"", comment.char = "")
  if (nrow(once) > 0 && identical(as.character(once[1, 1]), "nx")) {
    once <- once[-1, , drop = FALSE]
  }
  once
}

read_one_result <- function(nx, ny, g, beta0, rho_1, exposure_rho, xy_rho,
                            tau, siga1, siga2, weakp,
                            F_THRESHOLD, P_THRESHOLD) {
  ccname <- build_ccname(nx, ny)
  filename <- find_result_file(nx, ny, g, beta0, rho_1, exposure_rho,
                               xy_rho, tau, siga1, siga2, weakp,
                               F_THRESHOLD, P_THRESHOLD)
  if (is.na(filename)) {
    return(NULL)
  }

  once <- read_result_csv(filename, length(ccname))
  if (ncol(once) != length(ccname)) {
    stop(sprintf("Column mismatch in %s: expected %s, got %s",
                 filename, length(ccname), ncol(once)))
  }

  colnames(once) <- ccname
  once$datanum <- paste0(once$nx, "-", once$ny)
  once[intersect(c(keep_names(), "datanum"), names(once))]
}

result_all <- NULL
missing_files <- character()

for (tau in taub) {
  for (nx in nxb) {
    for (ny in nyb) {
      for (siga_i in seq_len(nrow(siga_grid))) {
        siga1 <- siga_grid$siga1[siga_i]
        siga2 <- siga_grid$siga2[siga_i]
        for (weakp_i in seq_along(weakp)) {
          weakp_now <- weakp[weakp_i]
          for (g in gb) {
            for (beta0 in beta0b) {
              for (rho_i in seq_along(rho_1b)) {
                rho_1 <- rho_1b[rho_i]
                for (exposure_rho in exposure_rhob) {
                  for (xy_rho in xy_rhob) {
                    once <- read_one_result(
                      nx, ny, g, beta0, rho_1, exposure_rho, xy_rho,
                      tau, siga1, siga2, weakp_now,
                      F_THRESHOLD, P_THRESHOLD
                    )
                    if (is.null(once)) {
                      missing_files <- c(
                        missing_files,
                        result_filename(nx, ny, g, beta0, rho_1,
                                        exposure_rho, xy_rho, tau,
                                        siga1, siga2, weakp_now,
                                        F_THRESHOLD, P_THRESHOLD)
                      )
                    } else {
                      result_all <- rbind(result_all, once)
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}

if (length(missing_files) > 0) {
  message("Missing result files:\n", paste(missing_files, collapse = "\n"))
  if (stop_if_missing) {
    stop("Some target result files are missing. Run the simulation first, or set stop_if_missing <- FALSE for a preview.")
  }
}

if (is.null(result_all) || nrow(result_all) == 0) {
  stop("No result files were read.")
}

method_levels <- c(
  "meta_res_dIVW", "meta_res_FIVW", "meta_res_RIVW", "meta_res_WME",
  "dIVW_meta", "FIVW_meta", "RIVW_meta", "WME_meta",
  "JointMR"
)

joint_beta <- ifelse(as.numeric(result_all$siga2) > 0,
                     result_all$JointMR_RN_ple_plugin_beta,
                     result_all$JointMR_RN_beta)
joint_pval <- ifelse(as.numeric(result_all$siga2) > 0,
                     result_all$JointMR_RN_ple_plugin_pval,
                     result_all$JointMR_RN_pval)

ress <- data.frame(
  beta = c(result_all$meta_res_dIVW_beta,
           result_all$meta_res_IVW_fixed_beta,
           result_all$meta_res_IVW_random_beta,
           result_all$meta_res_wme_beta,
           result_all$MR_meta_old_dIVW_beta,
           result_all$MR_meta_old_ivw_fixed_beta,
           result_all$MR_meta_old_ivw_random_beta,
           result_all$MR_meta_old_wme_beta,
           joint_beta),
  pval = c(result_all$meta_res_dIVW_pval,
           result_all$meta_res_IVW_fixed_pval,
           result_all$meta_res_IVW_random_pval,
           result_all$meta_res_wme_pval,
           result_all$MR_meta_old_dIVW_pval,
           result_all$MR_meta_old_ivw_fixed_pval,
           result_all$MR_meta_old_ivw_random_pval,
           result_all$MR_meta_old_wme_pval,
           joint_pval),
  beta1 = rep(result_all$beta0, length(method_levels)),
  rho_y = rep(result_all$rho_mean, length(method_levels)),
  rho_x = rep(result_all$exposure_rho, length(method_levels)),
  rho_xy = rep(result_all$xy_rho, length(method_levels)),
  tau = rep(result_all$tau, length(method_levels)),
  siga1 = rep(result_all$siga1, length(method_levels)),
  siga2 = rep(result_all$siga2, length(method_levels)),
  datanum = rep(result_all$datanum, length(method_levels)),
  weakp = rep(result_all$weakp, length(method_levels)),
  methods = rep(method_levels, each = nrow(result_all))
)

num_cols <- c("beta", "pval", "beta1", "rho_y", "rho_x", "rho_xy",
              "tau", "siga1", "siga2", "weakp")
ress[num_cols] <- lapply(ress[num_cols], as.numeric)

rho_y_levels <- sort(unique(ress$rho_y))
rho_x_levels <- sort(unique(ress$rho_x))
rho_xy_levels <- sort(unique(ress$rho_xy))
tau_levels <- sort(unique(ress$tau))
pleio_levels <- sort(unique(ress$siga2))

ress$beta0 <- factor(ress$beta1,
                     levels = sort(unique(beta0b)),
                     ordered = TRUE,
                     labels = paste0("causal effect = ",
                                     sprintf("%.2f", sort(unique(beta0b)))))
ress$rho_y_f <- factor(ress$rho_y,
                       levels = rho_y_levels,
                       ordered = TRUE,
                       labels = paste0("YY rho = ", rho_y_levels))
ress$rho_x_f <- factor(ress$rho_x,
                       levels = rho_x_levels,
                       ordered = TRUE,
                       labels = paste0("XX rho = ", rho_x_levels))
ress$rho_xy_f <- factor(ress$rho_xy,
                        levels = rho_xy_levels,
                        ordered = TRUE,
                        labels = paste0("XY rho = ", rho_xy_levels))
ress$tau_f <- factor(ress$tau,
                     levels = tau_levels,
                     ordered = TRUE,
                     labels = paste0("tau = ", tau_levels))
ress$pleio_f <- factor(ress$siga2,
                       levels = pleio_levels,
                       ordered = TRUE,
                       labels = ifelse(pleio_levels == 0,
                                       "pleiotropy = 0",
                                       paste0("pleiotropy = ", pleio_levels)))
ress$methods <- factor(ress$methods, levels = method_levels)

coloo <- c(brewer.pal(9, "Set3"))[c(1:3, 5:9, 4)]

theme_compare <- theme_bw() +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5),
        axis.text = element_text(size = 11),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 12),
        plot.title = element_text(size = 16, hjust = 0.5),
        legend.position = "bottom",
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 11),
        legend.key.size = unit(1.2, "lines"))

scenario_title <- function(tau_value, pleio_value) {
  paste0("tau = ", tau_value, ", pleiotropy = ", pleio_value)
}

scenario_suffix <- function(tau_value, pleio_value) {
  paste0("_tau_", tau_value, "_siga2_", pleio_value)
}

####################### 1. Boxplot ##########################################

for (tau_value in tau_levels) {
  for (pleio_value in pleio_levels) {
    pdf(paste0("weak_0428_relaxed_nome_boxplot_by_XX_XY_rho",
               scenario_suffix(tau_value, pleio_value), ".pdf"),
        width = 22, height = 14)
    box_data <- ress[ress$tau == tau_value & ress$siga2 == pleio_value &
                       is.finite(ress$beta), ]
    box_hline <- unique(box_data[c("beta0", "beta1", "rho_xy_f", "rho_x_f")])

    p_box <- ggplot(box_data, aes(x = rho_y_f, y = beta)) +
      geom_boxplot(aes(fill = methods),
                   position = position_dodge2(width = 0.85,
                                              preserve = "single"),
                   outlier.alpha = 0.25, linewidth = 0.25) +
      geom_hline(data = box_hline,
                 aes(yintercept = beta1),
                 linetype = "dashed", color = "darkred") +
      scale_fill_manual(values = coloo) +
      facet_grid(beta0 + rho_xy_f ~ rho_x_f, scales = "free_y") +
      labs(title = scenario_title(tau_value, pleio_value),
           x = "outcome/outcome correlation",
           y = "Estimation", fill = "Method") +
      guides(fill = guide_legend(nrow = 2)) +
      theme_compare
    print(p_box)
    dev.off()
  }
}

####################### 2. Q-Q plot: beta0 = 0 ##############################

method_names <- as.character(method_levels)

for (tau_value in tau_levels) {
  for (pleio_value in pleio_levels) {
    pdf(paste0("weak_0428_relaxed_nome_qq_by_XX_XY_rho",
               scenario_suffix(tau_value, pleio_value), ".pdf"),
        width = 24, height = 12)
    plot_list <- list()
    for (xy_value in rho_xy_levels) {
      for (x_value in rho_x_levels) {
        for (y_value in rho_y_levels) {
          once <- ress[ress$beta1 == 0 &
                         ress$tau == tau_value &
                         ress$siga2 == pleio_value &
                         ress$rho_y == y_value &
                         ress$rho_x == x_value &
                         ress$rho_xy == xy_value, ]

          p_value <- lapply(method_names, function(method) {
            as.numeric(once$pval[once$methods == method])
          })
          complete_p <- Reduce(
            `&`,
            lapply(p_value, function(p) is.finite(p) & p > 0 & p <= 1)
          )
          p_value <- lapply(p_value, function(p) p[complete_p])

          legend_pos <- if (xy_value == tail(rho_xy_levels, 1) &&
                              x_value == tail(rho_x_levels, 1) &&
                              y_value == tail(rho_y_levels, 1)) {
            "right"
          } else {
            "none"
          }

          plot_list[[length(plot_list) + 1]] <- suppressMessages(
            qq.plot(p_value, legend_label = method_names) +
              scale_y_continuous() +
              coord_cartesian(ylim = c(0, 4))
          ) +
            labs(title = paste0("YY rho = ", y_value,
                                "\nXX rho = ", x_value,
                                ", XY rho = ", xy_value)) +
            theme(plot.title = element_text(size = 16),
                  axis.title.x = element_text(size = 14),
                  axis.title.y = element_text(size = 14),
                  legend.position = legend_pos,
                  legend.title = element_text(size = 20),
                  legend.text = element_text(size = 18),
                  legend.key.size = unit(2.0, "lines"))
        }
      }
    }

    print(
      wrap_plots(plot_list, ncol = length(rho_y_levels) * length(rho_x_levels)) +
        plot_annotation(title = scenario_title(tau_value, pleio_value)) &
        theme(plot.title = element_text(hjust = 0.5))
    )
    dev.off()
  }
}

####################### 3. Power plot: beta0 = 0.05 #########################

calc_empirical_threshold <- function(dat_h0) {
  split_dat <- split(
    dat_h0,
    list(dat_h0$tau_f, dat_h0$pleio_f, dat_h0$rho_y_f, dat_h0$rho_x_f,
         dat_h0$rho_xy_f, dat_h0$methods),
    drop = TRUE
  )
  out <- lapply(split_dat, function(x) {
    p <- sort(x$pval[is.finite(x$pval) & x$pval > 0 & x$pval <= 1])
    if (length(p) == 0) {
      return(NULL)
    }
    index <- floor(0.05 * length(p)) + 1
    index <- min(index, length(p))
    data.frame(
      tau_f = x$tau_f[1],
      pleio_f = x$pleio_f[1],
      rho_y_f = x$rho_y_f[1],
      rho_x_f = x$rho_x_f[1],
      rho_xy_f = x$rho_xy_f[1],
      methods = x$methods[1],
      thres = p[index]
    )
  })
  do.call(rbind, out)
}

thres_df <- calc_empirical_threshold(ress[ress$beta1 == 0, ])
power_input <- ress[ress$beta1 == 0.05, ]

power_df <- power_input %>%
  left_join(thres_df,
            by = c("tau_f", "pleio_f", "rho_y_f", "rho_x_f",
                   "rho_xy_f", "methods")) %>%
  group_by(tau_f, pleio_f, rho_y_f, rho_x_f, rho_xy_f, methods) %>%
  summarise(power = mean(pval < thres, na.rm = TRUE), .groups = "drop")

for (tau_value in tau_levels) {
  for (pleio_value in pleio_levels) {
    pdf(paste0("weak_0428_relaxed_nome_power_by_XX_XY_rho",
               scenario_suffix(tau_value, pleio_value), ".pdf"),
        width = 24, height = 12)
    plot_list <- list()
    tau_label <- paste0("tau = ", tau_value)
    pleio_label <- ifelse(pleio_value == 0,
                          "pleiotropy = 0",
                          paste0("pleiotropy = ", pleio_value))
    for (xy_value in rho_xy_levels) {
      for (x_value in rho_x_levels) {
        for (y_value in rho_y_levels) {
          y_label <- paste0("YY rho = ", y_value)
          x_label <- paste0("XX rho = ", x_value)
          xy_label <- paste0("XY rho = ", xy_value)
          once_power <- power_df[power_df$tau_f == tau_label &
                                   power_df$pleio_f == pleio_label &
                                   power_df$rho_y_f == y_label &
                                   power_df$rho_x_f == x_label &
                                   power_df$rho_xy_f == xy_label, ]

          powt <- rep(NA_real_, length(method_names))
          for (i in seq_along(method_names)) {
            powt[i] <- once_power$power[once_power$methods == method_names[i]]
          }

          power_plot <- hist.plot(power = round(powt, 2),
                                  legend_label = method_names)
          power_plot$layers[[length(power_plot$layers)]] <- geom_text(
            mapping = aes(label = sprintf("%.2f", pow))
          )

          legend_pos <- if (xy_value == tail(rho_xy_levels, 1) &&
                              x_value == tail(rho_x_levels, 1) &&
                              y_value == tail(rho_y_levels, 1)) {
            "right"
          } else {
            "none"
          }

          plot_list[[length(plot_list) + 1]] <- power_plot +
            labs(title = paste0("YY rho = ", y_value,
                                "\nXX rho = ", x_value,
                                ", XY rho = ", xy_value)) +
            theme(plot.title = element_text(size = 16),
                  axis.text.x = element_blank(),
                  axis.title.x = element_text(size = 14),
                  axis.title.y = element_text(size = 14),
                  legend.position = legend_pos,
                  legend.title = element_text(size = 20),
                  legend.text = element_text(size = 18),
                  legend.key.size = unit(2.0, "lines"))
        }
      }
    }

    print(
      wrap_plots(plot_list, ncol = length(rho_y_levels) * length(rho_x_levels)) +
        plot_annotation(title = scenario_title(tau_value, pleio_value)) &
        theme(plot.title = element_text(hjust = 0.5))
    )
    dev.off()
  }
}

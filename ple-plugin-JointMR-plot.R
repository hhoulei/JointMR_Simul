source("F://wsj-doctor//2021-跨种族MR//ago//20230927--no pleo//QQ.plot_20230926.R")
source("F://wsj-doctor//2021-跨种族MR//ago//20230927--no pleo//histplot_20230926.R")

library(ggplot2)
library(RColorBrewer)
library(gridExtra)
library(cowplot)
library(patchwork)
library(outliers)
library(dplyr)

######### Pleio-adjusted JointMR: compare tau=0 and tau=0.2 under nx=3, ny=5 #######
setwd("F://wsj-doctor//2024-MR-META//genemo一修//simulation addition//ple-new")

gb <- c(200)
beta0b <- c(0, 0.05)
taub <- c(0, 0.2)
rho_1b <- c(0.5, 0.8, 0.85, 0.9, 0.95)
rho_2b <- c(0.5, 0.8, 0.85, 0.9, 0.95)
siga1b <- c(0)
siga2b <- c(0.03)
nxb <- c(3)
nyb <- c(5)
weakp <- 0
F_THRESHOLD <- 10
P_THRESHOLD <- 5e-8

make_colnames <- function(nx, ny) {
  methods <- c("dIVW", "IVW_fixed", "IVW_random", "wme")
  all_beta_names <- unlist(lapply(methods, function(method) {
    paste0("MR_", method, "_", rep("beta", nx * ny), seq_len(nx * ny))
  }))
  all_se_names <- unlist(lapply(methods, function(method) {
    paste0("MR_", method, "_", rep("se", nx * ny), seq_len(nx * ny))
  }))
  all_pval_names <- unlist(lapply(methods, function(method) {
    paste0("MR_", method, "_", rep("pval", nx * ny), seq_len(nx * ny))
  }))

  c(
    "nx", "ny", "g", "beta0", "rho_mean", "tau", "siga1",
    "siga2", "SNP_gaws_meta", "SNPMR", "SNPnew",
    paste0("meta_res_dIVW_", c("beta", "se", "pval")),
    paste0("meta_res_IVW_fixed_", c("beta", "se", "pval")),
    paste0("meta_res_IVW_random_", c("beta", "se", "pval")),
    paste0("meta_res_wme_", c("beta", "se", "pval")),
    all_beta_names, all_se_names, all_pval_names,
    paste0("MR_meta_old_dIVW_", c("beta", "se", "pval")),
    paste0("MR_meta_old_ivw_fixed_", c("beta", "se", "pval")),
    paste0("MR_meta_old_ivw_random_", c("beta", "se", "pval")),
    paste0("MR_meta_old_wme_", c("beta", "se", "pval")),
    paste0("MRmeta_", c("beta", "se", "pval")),
    "F_tat", "weakp",
    "alpha0_hat_mean", "alpha0_se_mean", "alpha0_true_mean",
    paste0("alpha0_hat_Y", seq_len(ny)),
    paste0("alpha0_se_Y", seq_len(ny)),
    paste0("alpha0_true_Y", seq_len(ny))
  )
}

tiquname <- c(
  "nx", "ny", "g", "beta0", "rho_mean", "tau", "siga1",
  "siga2", "SNP_gaws_meta", "SNPMR", "SNPnew",
  paste0("meta_res_dIVW_", c("beta", "se", "pval")),
  paste0("meta_res_IVW_fixed_", c("beta", "se", "pval")),
  paste0("meta_res_IVW_random_", c("beta", "se", "pval")),
  paste0("meta_res_wme_", c("beta", "se", "pval")),
  paste0("MR_meta_old_dIVW_", c("beta", "se", "pval")),
  paste0("MR_meta_old_ivw_fixed_", c("beta", "se", "pval")),
  paste0("MR_meta_old_ivw_random_", c("beta", "se", "pval")),
  paste0("MR_meta_old_wme_", c("beta", "se", "pval")),
  paste0("MRmeta_", c("beta", "se", "pval")),
  "alpha0_hat_mean", "alpha0_se_mean", "alpha0_true_mean",
  "datanum", "F_tat", "weakp"
)

read_one_result <- function(nx, ny, g, beta0, rho_1, tau, siga1, siga2, weakp) {
  ccname <- make_colnames(nx, ny)
  filename <- paste0(
    "Sum_nx_", nx, "_ny_", ny, "_g_", g, "_beta0_", beta0,
    "_rho1_", rho_1, "_tau_", tau, "_siga1_", siga1,
    "_siga2_", siga2, "_weakp_", weakp,
    "_F_THRESHOLD_", F_THRESHOLD,
    "_P_THRESHOLD_", P_THRESHOLD,
    "_0501_bootstrap_ple_plugin_delta.csv"
  )

  if (!file.exists(filename)) {
    warning("Missing file: ", filename)
    return(NULL)
  }

  once <- read.csv(filename, header = TRUE, check.names = FALSE)
  if (ncol(once) != length(ccname)) {
    once <- read.csv(filename, header = FALSE, check.names = FALSE)
  }
  if (nrow(once) > 0 && identical(as.character(once[1, 1]), ccname[1])) {
    once <- once[-1, , drop = FALSE]
  }
  if (ncol(once) != length(ccname)) {
    stop(
      "Column number mismatch in ", filename,
      ": observed ", ncol(once), ", expected ", length(ccname)
    )
  }

  colnames(once) <- ccname
  once$datanum <- paste0(once$nx, "-", once$ny)
  once[tiquname]
}

result_all <- NULL
for (yh in seq_along(taub)) {
  tau <- taub[[yh]]
  for (th in seq_along(gb)) {
    g <- gb[[th]]
    for (xh in seq_along(siga1b)) {
      siga1 <- siga1b[[xh]]
      siga2 <- siga2b[[xh]]
      for (uh in seq_along(beta0b)) {
        beta0 <- beta0b[[uh]]
        for (rh in seq_along(rho_1b)) {
          rho_1 <- rho_1b[[rh]]
          rho_2 <- rho_2b[[rh]]
          for (zh in seq_along(nxb)) {
            nx <- nxb[[zh]]
            for (qh in seq_along(nyb)) {
              ny <- nyb[[qh]]
              once <- read_one_result(nx, ny, g, beta0, rho_1, tau,
                                      siga1, siga2, weakp)
              if (!is.null(once)) {
                result_all <- rbind(result_all, once)
              }
            }
          }
        }
      }
    }
  }
}

if (is.null(result_all) || nrow(result_all) == 0) {
  stop("No simulation result files were loaded.")
}

ress <- data.frame(
  beta = c(
    result_all$meta_res_dIVW_beta,
    result_all$meta_res_IVW_fixed_beta,
    result_all$meta_res_IVW_random_beta,
    result_all$meta_res_wme_beta,
    result_all$MR_meta_old_dIVW_beta,
    result_all$MR_meta_old_ivw_fixed_beta,
    result_all$MR_meta_old_ivw_random_beta,
    result_all$MR_meta_old_wme_beta,
    result_all$MRmeta_beta
  ),
  pval = c(
    result_all$meta_res_dIVW_pval,
    result_all$meta_res_IVW_fixed_pval,
    result_all$meta_res_IVW_random_pval,
    result_all$meta_res_wme_pval,
    result_all$MR_meta_old_dIVW_pval,
    result_all$MR_meta_old_ivw_fixed_pval,
    result_all$MR_meta_old_ivw_random_pval,
    result_all$MR_meta_old_wme_pval,
    result_all$MRmeta_pval
  ),
  beta1 = rep(result_all$beta0, 9),
  rho = rep(result_all$rho_mean, 9),
  tau = rep(result_all$tau, 9),
  datanum = rep(result_all$datanum, 9),
  weakp = rep(result_all$weakp, 9),
  alpha0_hat_mean = rep(result_all$alpha0_hat_mean, 9),
  alpha0_true_mean = rep(result_all$alpha0_true_mean, 9),
  methods = c(rep(c(
    "meta_dIVW", "meta_FIVW", "meta_RIVW",
    "meta_WME", "dIVW_meta", "FIVW_meta",
    "RIVW_meta", "WME_meta", "JointMR"
  ), each = nrow(result_all)))
)

ress$beta <- as.numeric(ress$beta)
ress$pval <- as.numeric(ress$pval)
ress$beta1 <- as.numeric(ress$beta1)
ress$rho <- as.numeric(ress$rho)
ress$tau <- as.numeric(ress$tau)
ress$alpha0_hat_mean <- as.numeric(ress$alpha0_hat_mean)
ress$alpha0_true_mean <- as.numeric(ress$alpha0_true_mean)
ress <- ress[is.finite(ress$beta1) & is.finite(ress$rho) &
               is.finite(ress$tau), ]

ress$beta0 <- factor(
  ress$beta1,
  levels = beta0b,
  ordered = TRUE,
  labels = c("causal effect = 0.00",
             "causal effect = 0.05")
)
ress$rho_f <- factor(
  ress$rho,
  levels = rho_1b,
  ordered = TRUE,
  labels = c("0.5", "0.8", "0.85", "0.9", "0.95")
)
ress$tau_f <- factor(
  ress$tau,
  levels = taub,
  ordered = TRUE,
  labels = c("tau = 0 (fixed effect)",
             "tau = 0.2 (random effect)")
)
ress$datanum <- factor(
  ress$datanum,
  levels = c("3-5"),
  ordered = TRUE,
  labels = c("number of X-Y datasets = 3-5")
)
ress$methods <- factor(
  ress$methods,
  levels = c(
    "meta_dIVW", "meta_FIVW", "meta_RIVW",
    "meta_WME", "dIVW_meta", "FIVW_meta",
    "RIVW_meta", "WME_meta", "JointMR"
  )
)

method_names <- levels(ress$methods)
coloo <- c(brewer.pal(9, "Set3"))[c(1:3, 5:9, 4)]

####################### boxplot #################################
box_data <- ress[is.finite(ress$beta), ]

for (t in seq_along(taub)) {
  tau_value <- taub[t]
  tau_data <- box_data[box_data$tau == tau_value, ]
  beta_lines <- unique(tau_data[, c("beta0", "beta1")])

  p_box <- ggplot(data = tau_data, aes(x = rho_f, y = beta)) +
    geom_boxplot(aes(fill = methods), outlier.alpha = 0.5) +
    scale_fill_manual(values = coloo) +
    geom_hline(
      data = beta_lines,
      aes(yintercept = beta1),
      linetype = "dashed",
      color = "darkred"
    ) +
    theme_bw() +
    coord_cartesian(ylim = range(beta0b) + c(-0.25, 0.25)) +
    xlab("dataset correlation") +
    ylab("Estimation") +
    labs(tag = "(A)", title = levels(ress$tau_f)[t]) +
    facet_grid(. ~ beta0, scales = "free_x") +
    theme(
      axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5),
      strip.text.x = element_text(size = 35),
      strip.text.y = element_text(size = 35),
      plot.title = element_text(size = 35, hjust = 0.5, vjust = 0.5),
      axis.text = element_text(size = 25),
      axis.title = element_text(size = 30),
      legend.text = element_text(size = 25),
      legend.title = element_text(size = 25),
      plot.tag = element_text(size = 40),
      legend.spacing = unit(2, "cm"),
      legend.key.size = unit(3, "lines")
    ) +
    guides(fill = guide_legend(title = "Method"))

  ggsave(paste0("ple_plugin_tau_compare_3-5_Estimation_boxplot_tau_", tau_value, ".pdf"),
         plot = p_box, width = 30, height = 8)
}

####################### Q-Q plot #################################
rho_meanb <- rho_1b
rho_meanname <- c(
  "data correlation is 0.5",
  "data correlation is 0.8",
  "data correlation is 0.85",
  "data correlation is 0.9",
  "data correlation is 0.95"
)
tauname <- c("tau = 0 (fixed effect)",
             "tau = 0.2 (random effect)")

for (m in 1:1) {
  for (t in seq_along(taub)) {
    plot_list <- list()
    for (n in seq_along(rho_meanb)) {
      once <- ress[(ress$beta1 == beta0b[m] &
                      ress$rho == rho_meanb[n] &
                      ress$tau == taub[t]), ]

      p_value <- list()
      for (i in seq_along(method_names)) {
        p_once <- as.numeric(once$pval[once$methods == method_names[i]])
        p_value[[i]] <- p_once
      }
      complete_p <- Reduce(
        `&`,
        lapply(p_value, function(p) is.finite(p) & p > 0 & p <= 1)
      )
      p_value <- lapply(p_value, function(p) p[complete_p])

      plot <- suppressMessages(
        qq.plot(p_value, legend_label = method_names) +
          scale_y_continuous() +
          coord_cartesian(ylim = c(0, 4))
      ) +
        labs(title = paste0(rho_meanname[n])) +
        theme(
          plot.title = element_text(size = 25),
          plot.tag = element_text(size = 23),
          axis.title.x = element_text(size = 20),
          axis.title.y = element_text(size = 20)
        )

      plot_list[[length(plot_list) + 1]] <- plot
    }

    big_plot <- wrap_plots(plot_list, ncol = length(rho_meanb)) +
      plot_layout(guides = "collect") &
      theme(
        legend.position = "right",
        legend.spacing = unit(1.5, "cm"),
        legend.key.size = unit(3, "lines"),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20)
      )

    title_plot <- ggplot() +
      theme_void() +
      annotate("text", x = 0, y = 0, label = tauname[t], size = 15) +
      theme(plot.margin = margin(b = -15, unit = "pt")) +
      labs(tag = "(B)")

    final_plot <- title_plot + big_plot +
      plot_layout(ncol = 1, heights = c(0.3, 1)) &
      theme(plot.tag = element_text(size = 40))

    pdf(paste0("ple_plugin_tau_compare_3-5_qq_tau_", taub[t], ".pdf"),
        width = 30, height = 8)
    print(final_plot)
    dev.off()
  }
}

################# power plot ########################
thres_df <- expand.grid(methods = method_names, rho = rho_meanb, tau = taub)
thres_df$thres <- NA

for (r in rho_meanb) {
  for (tt in taub) {
    once_H0 <- ress[ress$beta1 == beta0b[1] &
                      ress$rho == r &
                      ress$tau == tt, ]
    for (mth in method_names) {
      pvals <- as.numeric(once_H0$pval[once_H0$methods == mth])
      pvals <- pvals[is.finite(pvals) & pvals > 0 & pvals <= 1]
      if (length(pvals) > 0) {
        pvals_sorted <- sort(pvals)
        index <- floor(0.05 * length(pvals_sorted)) + 1
        thres_value <- if (index > length(pvals_sorted)) {
          pvals_sorted[length(pvals_sorted)]
        } else {
          pvals_sorted[index]
        }
        thres_df$thres[thres_df$methods == mth &
                         thres_df$rho == r &
                         thres_df$tau == tt] <- thres_value
      } else {
        stop(sprintf("No H0 data for method %s, rho=%s, tau=%s", mth, r, tt))
      }
    }
  }
}

for (m in 2:length(beta0b)) {
  for (t in seq_along(taub)) {
    plot_list <- list()

    for (n in seq_along(rho_meanb)) {
      once <- ress[(ress$beta1 == beta0b[m] &
                      ress$rho == rho_meanb[n] &
                      ress$tau == taub[t]), ]

      powt <- NULL
      for (i in seq_along(method_names)) {
        once1 <- once[once$methods == method_names[i], ]
        thres_value <- thres_df$thres[thres_df$methods == method_names[i] &
                                        thres_df$rho == rho_meanb[n] &
                                        thres_df$tau == taub[t]]
        p_once <- as.numeric(once1$pval)
        p_once <- p_once[is.finite(p_once) & p_once > 0 & p_once <= 1]
        pow <- if (length(p_once) > 0 && is.finite(thres_value)) {
          round(mean(p_once < thres_value, na.rm = TRUE), digits = 2)
        } else {
          NA_real_
        }
        powt <- c(powt, pow)
      }

      plot <- hist.plot(power = powt, legend_label = method_names) +
        labs(title = paste0(rho_meanname[n])) +
        theme(
          plot.title = element_text(size = 25),
          plot.tag = element_text(size = 23),
          axis.text.x = element_blank(),
          axis.title.x = element_text(size = 20),
          axis.title.y = element_text(size = 20),
          legend.position = "none"
        )

      plot_list[[length(plot_list) + 1]] <- plot
    }

    big_plot <- wrap_plots(plot_list, ncol = length(rho_meanb)) +
      plot_layout(guides = "collect") &
      theme(
        legend.position = "right",
        legend.spacing = unit(1.5, "cm"),
        legend.key.size = unit(3, "lines"),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20)
      )

    title_plot <- ggplot() +
      theme_void() +
      annotate("text", x = 0, y = 0, label = tauname[t], size = 15) +
      theme(plot.margin = margin(b = -15, unit = "pt")) +
      labs(tag = "(C)")

    final_plot <- title_plot + big_plot +
      plot_layout(ncol = 1, heights = c(0.3, 1)) &
      theme(plot.tag = element_text(size = 40))

    pdf(paste0("ple_plugin_tau_compare_3-5_power_tau_", taub[t], ".pdf"),
        width = 30, height = 8)
    print(final_plot)
    dev.off()
  }
}

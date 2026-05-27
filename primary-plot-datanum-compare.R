library(ggplot2)
library(RColorBrewer)
library(patchwork)
library(dplyr)

source('F://wsj-doctor//2021-跨种族MR//ago//20230927--no pleo//QQ.plot_20230926.R')
source('F://wsj-doctor//2021-跨种族MR//ago//20230927--no pleo//histplot_20230926.R')

if (basename(getwd()) != "weak" && dir.exists("weak")) {
  setwd("weak")
}
setwd("F://wsj-doctor//2024-MR-META//weak")

#### Five X-Y database-pair comparison: 1-5, 2-2, 2-5, 3-3, 3-5 ####

gb <- c(200)
beta0b <- c(0, 0.05, 0.1)
taub <- c(0)
rho_1b <- c(0.5, 0.8, 0.85, 0.9, 0.95)
rho_2b <- rho_1b
siga1b <- c(0)
siga2b <- c(0)
weakp <- 0
F_THRESHOLD <- 10
P_THRESHOLD <- 5e-08

target_pairs <- data.frame(
  nx = c(1, 2, 2, 3, 3),
  ny = c(5, 2, 5, 3, 5)
)

# TRUE is recommended for the final five-pair figure. Set to FALSE while
# only part of the simulations have finished and you want a quick preview.
stop_if_missing <- TRUE

tiquname <- c('nx', 'ny', 'g', 'beta0', "rho_mean", 'tau', "siga1",
              "siga2", "SNP_gaws_meta", "SNPMR", "SNPnew",
              paste0('meta_dIVW_', c('beta', 'se', 'pval')),
              paste0('meta_IVW_fixed_', c('beta', 'se', 'pval')),
              paste0('meta_IVW_random_', c('beta', 'se', 'pval')),
              paste0('meta_WME_', c('beta', 'se', 'pval')),
              paste0('dIVW_meta_', c('beta', 'se', 'pval')),
              paste0('IVW_fixed_meta_', c('beta', 'se', 'pval')),
              paste0('IVW_random_meta_', c('beta', 'se', 'pval')),
              paste0('WME_meta_', c('beta', 'se', 'pval')),
              paste0('JointMR_', c('beta', 'se', 'pval')),
              "datanum", "F_tat", "weakp")

make_colnames <- function(nx, ny) {
  methods <- c("dIVW", "IVW_fixed", "IVW_random", "WME")
  all_beta_names <- unlist(lapply(methods, function(method) {
    paste0("MR_", method, "_", rep("beta", nx * ny), seq_len(nx * ny))
  }))
  all_se_names <- unlist(lapply(methods, function(method) {
    paste0("MR_", method, "_", rep("se", nx * ny), seq_len(nx * ny))
  }))
  all_pval_names <- unlist(lapply(methods, function(method) {
    paste0("MR_", method, "_", rep("pval", nx * ny), seq_len(nx * ny))
  }))

  c('nx', 'ny', 'g', 'beta0', "rho_mean", 'tau', "siga1",
    "siga2", "SNP_gaws_meta", "SNPMR", "SNPnew",
    paste0('meta_dIVW_', c('beta', 'se', 'pval')),
    paste0('meta_IVW_fixed_', c('beta', 'se', 'pval')),
    paste0('meta_IVW_random_', c('beta', 'se', 'pval')),
    paste0('meta_WME_', c('beta', 'se', 'pval')),
    all_beta_names, all_se_names, all_pval_names,
    paste0('dIVW_meta_', c('beta', 'se', 'pval')),
    paste0('IVW_fixed_meta_', c('beta', 'se', 'pval')),
    paste0('IVW_random_meta_', c('beta', 'se', 'pval')),
    paste0('WME_meta_', c('beta', 'se', 'pval')),
    paste0('JointMR_', c('beta', 'se', 'pval')), "F_tat", "weakp")
}

result_filename <- function(nx, ny, g, beta0, rho_1, tau,
                            siga1, siga2, weakp,
                            F_THRESHOLD, P_THRESHOLD) {
  paste0("Sum_nx_", nx, "_ny_", ny, "_g_", g, "_beta0_", beta0,
         "_rho1_", rho_1, "_tau_", tau, "_siga1_", siga1,
         "_siga2_", siga2, "_weakp_", weakp,
         "_F_THRESHOLD_", F_THRESHOLD,
         "_P_THRESHOLD_", P_THRESHOLD, "_1128.csv")
}

find_result_file <- function(nx, ny, g, beta0, rho_1, tau,
                             siga1, siga2, weakp,
                             F_THRESHOLD, P_THRESHOLD) {
  filename <- result_filename(nx, ny, g, beta0, rho_1, tau,
                              siga1, siga2, weakp,
                              F_THRESHOLD, P_THRESHOLD)
  if (file.exists(filename)) {
    return(filename)
  }

  filename_nx2 <- sub("\\.csv$", "_nx2.csv", filename)
  if (file.exists(filename_nx2)) {
    return(filename_nx2)
  }

  # If you changed siga2 in weak-1128.R, or added a suffix such as _nx2,
  # this fallback finds the matching file.
  pattern <- paste0("Sum_nx_", nx, "_ny_", ny, "_g_", g,
                    "_beta0_", beta0, "_rho1_", rho_1,
                    "_tau_", tau, "_siga1_", siga1,
                    "_siga2_*_weakp_", weakp,
                    "_F_THRESHOLD_", F_THRESHOLD,
                    "_P_THRESHOLD_", P_THRESHOLD, "_1128*.csv")
  matches <- Sys.glob(pattern)
  if (length(matches) == 1) {
    return(matches)
  }

  NA_character_
}

read_one_result <- function(nx, ny, g, beta0, rho_1, tau,
                            siga1, siga2, weakp,
                            F_THRESHOLD, P_THRESHOLD) {
  ccname <- make_colnames(nx, ny)
  filename <- find_result_file(nx, ny, g, beta0, rho_1, tau,
                               siga1, siga2, weakp,
                               F_THRESHOLD, P_THRESHOLD)
  if (is.na(filename)) {
    return(NULL)
  }

  once <- read.csv(filename, header = FALSE)
  if (ncol(once) != length(ccname)) {
    stop(sprintf("Column mismatch in %s: expected %s, got %s",
                 filename, length(ccname), ncol(once)))
  }

  colnames(once) <- ccname
  once <- na.omit(once)
  once$datanum <- paste0(once$nx, "-", once$ny)
  once[tiquname]
}

result_all <- NULL
missing_files <- character()

for (pair_i in seq_len(nrow(target_pairs))) {
  nx <- target_pairs$nx[pair_i]
  ny <- target_pairs$ny[pair_i]
  for (tau in taub) {
    for (g in gb) {
      for (siga_i in seq_along(siga1b)) {
        siga1 <- siga1b[siga_i]
        siga2 <- siga2b[siga_i]
        for (beta0 in beta0b) {
          for (rho_i in seq_along(rho_1b)) {
            rho_1 <- rho_1b[rho_i]
            rho_2 <- rho_2b[rho_i]
            once <- read_one_result(nx, ny, g, beta0, rho_1, tau,
                                    siga1, siga2, weakp,
                                    F_THRESHOLD, P_THRESHOLD)
            if (is.null(once)) {
              missing_files <- c(
                missing_files,
                result_filename(nx, ny, g, beta0, rho_1, tau,
                                siga1, siga2, weakp,
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

if (length(missing_files) > 0) {
  message("Missing result files:\n", paste(missing_files, collapse = "\n"))
  if (stop_if_missing) {
    stop("Some target result files are missing. Run weak-1128.R first, or set stop_if_missing <- FALSE for a preview.")
  }
}

if (is.null(result_all) || nrow(result_all) == 0) {
  stop("No result files were read.")
}

method_levels <- c('meta_dIVW', 'meta_FIVW', 'meta_RIVW',
                   'meta_WME', "dIVW_meta", 'FIVW_meta',
                   'RIVW_meta', 'WME_meta', 'JointMR')

ress <- data.frame(
  beta = c(result_all$meta_dIVW_beta,
           result_all$meta_IVW_fixed_beta,
           result_all$meta_IVW_random_beta,
           result_all$meta_WME_beta,
           result_all$dIVW_meta_beta,
           result_all$IVW_fixed_meta_beta,
           result_all$IVW_random_meta_beta,
           result_all$WME_meta_beta,
           result_all$JointMR_beta),
  pval = c(result_all$meta_dIVW_pval,
           result_all$meta_IVW_fixed_pval,
           result_all$meta_IVW_random_pval,
           result_all$meta_WME_pval,
           result_all$dIVW_meta_pval,
           result_all$IVW_fixed_meta_pval,
           result_all$IVW_random_meta_pval,
           result_all$WME_meta_pval,
           result_all$JointMR_pval),
  beta1 = rep(result_all$beta0, length(method_levels)),
  rho = rep(result_all$rho_mean, length(method_levels)),
  tau = rep(result_all$tau, length(method_levels)),
  datanum = rep(result_all$datanum, length(method_levels)),
  weakp = rep(result_all$weakp, length(method_levels)),
  methods = rep(method_levels, each = nrow(result_all))
)

ress$beta <- as.numeric(ress$beta)
ress$pval <- as.numeric(ress$pval)
ress$beta1 <- as.numeric(ress$beta1)
ress$rho <- as.numeric(ress$rho)
ress$tau <- as.numeric(ress$tau)

datanum_levels <- paste0(target_pairs$nx, "-", target_pairs$ny)

ress$beta0 <- factor(ress$beta1,
                     levels = beta0b,
                     ordered = TRUE,
                     labels = paste0("causal effect = ",
                                     sprintf("%.2f", beta0b)))
ress$rho_f <- factor(ress$rho,
                     levels = rho_1b,
                     ordered = TRUE,
                     labels = as.character(rho_1b))
ress$tau_f <- factor(ress$tau,
                     levels = taub,
                     ordered = TRUE,
                     labels = unname(c("0" = "tau = 0 (fixed effect)",
                                       "0.2" = "tau = 0.2 (random effect)")[
                                         as.character(taub)
                                       ]))
ress$datanum <- factor(ress$datanum,
                       levels = datanum_levels,
                       ordered = TRUE,
                       labels = paste0("X-Y datasets = ", datanum_levels))
ress$methods <- factor(ress$methods, levels = method_levels)

coloo <- c(brewer.pal(9, "Set3"))[c(1:3, 5:9, 4)]

theme_compare <- theme_bw() +
  theme(axis.text.x = element_text(angle = 35, vjust = 1, hjust = 1),
        strip.text = element_text(size = 13),
        plot.title = element_text(size = 18, hjust = 0.5),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 13),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 11),
        legend.key.size = unit(1.1, "lines"))

####################### 1. Boxplot: beta0 = 0 and 0.05 #######################
box_data <- ress[ress$beta1 %in% c(0, 0.05), ]
box_hline <- unique(box_data[c("beta0", "beta1", "datanum")])

p_box <- ggplot(box_data, aes(x = rho_f, y = beta)) +
  geom_boxplot(aes(fill = methods),
               position = position_dodge2(width = 0.85, preserve = "single"),
               outlier.alpha = 0.25, linewidth = 0.25) +
  geom_hline(data = box_hline,
             aes(yintercept = beta1),
             linetype = "dashed", color = "darkred") +
  scale_fill_manual(values = coloo) +
  facet_grid(beta0 ~ datanum, scales = "free_y") +
  xlab("dataset correlation") +
  ylab("Estimation") +
  guides(fill = guide_legend(title = "Method", nrow = 3)) +
  theme_compare +
  theme(axis.text.x = element_text(size = 22, angle = 0,
                                   vjust = 0.5, hjust = 0.5),
        axis.text.y = element_text(size = 22),
        axis.title.x = element_text(size = 30),
        axis.title.y = element_text(size = 30),
        strip.text.x = element_text(size = 24),
        strip.text.y = element_text(size = 24),
        legend.position = "bottom",
        legend.text = element_text(size = 22),
        legend.title = element_text(size = 24),
        legend.key.size = unit(2.2, "lines"))

ggsave("weak_five_datanum_boxplot_beta0_0_05.pdf",
       plot = p_box, width = 40, height = 18)

####################### 2. Q-Q plot: beta0 = 0 ###############################
method_names <- as.character(method_levels)
rho_meanname <- paste0("data correlation is ", rho_1b)
tau_target <- 0

qq_rows <- list()
for (d in seq_along(datanum_levels)) {
  rho_plots <- list()
  for (n in seq_along(rho_1b)) {
    once <- ress[ress$beta1 == 0 &
                   ress$rho == rho_1b[n] &
                   ress$tau == tau_target &
                   ress$datanum == levels(ress$datanum)[d], ]

    p_value <- list()
    for (i in seq_along(method_names)) {
      p_value[[i]] <- as.numeric(once$pval[once$methods == method_names[i]])
    }

    legend_pos <- if (d == 3 && n == length(rho_1b)) "right" else "none"

    rho_plots[[n]] <- qq.plot(p_value, legend_label = method_names) +
      labs(title = rho_meanname[n]) +
      theme(plot.title = element_text(size = 24),
            plot.tag = element_text(size = 22),
            axis.title.x = element_text(size = 20),
            axis.title.y = element_text(size = 20),
            legend.position = legend_pos,
            legend.title = element_text(size = 28),
            legend.text = element_text(size = 28),
            legend.key.size = unit(3.0, "lines"))
  }

  qq_rows[[d]] <- wrap_plots(rho_plots, ncol = length(rho_1b)) +
    plot_annotation(title = levels(ress$datanum)[d]) &
    theme(plot.title = element_text(size = 24, hjust = 0.5))
}

pdf("weak_five_datanum_qq.pdf",
    width = 40, height = 8 * length(datanum_levels))
print(wrap_plots(qq_rows, ncol = 1))
dev.off()

####################### 3. Power plot: beta0 = 0.05 ##########################
calc_empirical_threshold <- function(dat_h0) {
  split_dat <- split(dat_h0, list(dat_h0$datanum, dat_h0$rho_f,
                                  dat_h0$methods),
                     drop = TRUE)
  out <- lapply(split_dat, function(x) {
    p <- sort(x$pval[is.finite(x$pval) & x$pval > 0 & x$pval <= 1])
    if (length(p) == 0) {
      return(NULL)
    }
    index <- floor(0.05 * length(p)) + 1
    index <- min(index, length(p))
    data.frame(
      datanum = x$datanum[1],
      rho_f = x$rho_f[1],
      methods = x$methods[1],
      thres = p[index]
    )
  })
  do.call(rbind, out)
}

thres_df <- calc_empirical_threshold(ress[ress$beta1 == 0, ])
power_input <- ress[ress$beta1 == 0.05, ]

power_df <- power_input %>%
  left_join(thres_df, by = c("datanum", "rho_f", "methods")) %>%
  group_by(datanum, rho_f, methods) %>%
  summarise(power = mean(pval < thres, na.rm = TRUE), .groups = "drop")

power_rows <- list()
for (d in seq_along(datanum_levels)) {
  rho_plots <- list()
  for (n in seq_along(rho_1b)) {
    once_power <- power_df[power_df$rho_f == as.character(rho_1b[n]) &
                             power_df$datanum == levels(ress$datanum)[d], ]

    powt <- rep(NA_real_, length(method_names))
    for (i in seq_along(method_names)) {
      powt[i] <- once_power$power[once_power$methods == method_names[i]]
    }

    power_plot <- hist.plot(power = round(powt, 2),
                            legend_label = method_names)
    power_plot$layers[[length(power_plot$layers)]] <- geom_text(
      mapping = aes(label = sprintf("%.2f", pow))
    )

    legend_pos <- if (d == 3 && n == length(rho_1b)) "right" else "none"

    rho_plots[[n]] <- power_plot +
      labs(title = rho_meanname[n]) +
      theme(plot.title = element_text(size = 18),
            plot.tag = element_text(size = 16),
            axis.text.x = element_blank(),
            axis.title.x = element_text(size = 14),
            axis.title.y = element_text(size = 14),
            legend.position = legend_pos,
            legend.title = element_text(size = 28),
            legend.text = element_text(size = 28),
            legend.key.size = unit(3.0, "lines"))
  }

  power_rows[[d]] <- wrap_plots(rho_plots, ncol = length(rho_1b)) +
    plot_annotation(title = levels(ress$datanum)[d]) &
    theme(plot.title = element_text(size = 20, hjust = 0.5))
}

pdf("weak_five_datanum_power_beta0_0.05.pdf",
    width = 40, height = 8 * length(datanum_levels))
print(wrap_plots(power_rows, ncol = 1))
dev.off()
